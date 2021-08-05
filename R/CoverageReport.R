# Author: Walter Xie
# Accessed on 27 Jul 2021

#' @name Coverage
#' @title Coverage Report from Model Validation
#'
#' @description
#' The pipeline to generate a coverage report of the model validation.
#'
#' @details
#' \code{summariseTracesAndTrees} creates the statistical summary of
#' traces from a BEAST log file and tree statistics from BEAST tree log file,
#' and save them into tsv files.
#'
#' @param log.file  The BEAST log file containing posterior samples.
#' @param tree.file The BEAST tree log file containing posterior trees.
#'                  Default to NA, which will skip to report the tree statistics.
#' @param burn.in   The proportion of samples treated as the burn-in stage
#'                  during MCMC. Default to 0.1.
#' @param stats.fn.fun,tree.stats.fn.fun The function to create the result file names
#'                  from the given log file names. Set to NA, if you do not want to
#'                  create the files.
#' @keywords Coverage
#' @export
#' @examples
#' WD = file.path("~/WorkSpace/MyValidations/")
#' setwd(WD)
#' log.files = list.files(pattern = ".log")
#' for(lg in log.files) {
#'   # assume same file stem
#'   tree.file=paste0(sub('\\.log$', '', lg), ".trees")
#'   res <- summariseTracesAndTrees(lg, tree.file)
#' }
#'
#' # skip tree logs
#' summariseTracesAndTrees(lg, tree.file=NA)
#'
#' @rdname Coverage
summariseTracesAndTrees <- function(log.file, tree.file=NA, burn.in=0.1,
                                    stats.fn.fun=function(x){ paste0(sub('\\.log$', '', x), ".tsv") },
                                    tree.stats.fn.fun=function(x){ paste0(x, ".tsv") }) {
  require(tidyverse)

  cat("\nProcess ", log.file, "...\n")
  if (!file.exists(log.file)) stop("Cannot find file : ", log.file)

  # read MCMC log
  mcmc.log <- readMCMCLog(log.file)
  # get traces and remove burn in
  traces <- getTraces(mcmc.log, burn.in=burn.in)
  # get stats
  stats <- analyseTraces(traces)

  if (!is.function(stats.fn.fun)) write_tsv(stats, stats.fn.fun(log.file))

  # add tree stats
  tre.sta <- NULL
  if (!is.na(tree.file)) {
    if (!file.exists(tree.file)) stop("Cannot find file : ", tree.file)

    tre.sta.df <- readTrees(tree.file)
    tre.sta <- analyseTreeStats(tre.sta.df)

    # ? HPD95.lower.STATE_14150000
    #tre.stats$trace <- sub("\\.STATE.*$","",tre.stats$trace, ignore.case = T)
    if (!is.function(tree.stats.fn.fun))
      write_tsv(tre.sta, tree.stats.fn.fun(tree.file))
  }
  list(stats=stats, tree.stats=tre.sta)
}


selectResultByESS <- function(i.sta=0, i.end=99, prefix="sim",
                              file.steam.fun=function(prefix, i){ paste0(prefix, "_", i) }
                              file.postfix="tsv", tree.file.postfix="trees.tsv" ) {
  require(tidyverse)

  # record 100 files whose ESS >= 200
  tracesDF <- list()
  # if any of 100 has <200 ESS, then it will be replaced by one of extra 10
  lowESS <- tibble()
  etr <- 1
  for(i in i.sta:i.end) {
    fn <- file.steam.fun(prefix, i)
    fi <- paste0(fn, ".", file.postfix)
    stopifnot(file.exists(fi))
    tre.sta.fi <- paste0(fn, ".", tree.file.postfix)

    # "mean", "HPD95.lower", "HPD95.upper", "ESS"
    traces <- readTraces(fi) %>% select(trace, params)
    ESS <- traces %>% filter(trace == "ESS") %>% select(!trace)
    minESS <- min(ESS %>% as.numeric)
    cat(fi, ", min ESS = ", minESS, "\n")

    if (file.exists(tre.sta.fi)) {
      # add tree stats
      tre.sta <- try(read_tsv(tre.sta.fi)) %>%
        filter(trace %in% stats.name) %>% select(!trace)
      traces <- cbind(traces, tre.sta)
    }

    if (minESS < 200) {
      selected <- pullAValidResultFromExtraPool(fi, etr, ext.pool=extraStats, params=c("trace", params))
      # selected extra
      file.selected <- selected$file.selected
      fn.sel <- sub('\\.tsv$', '', file.selected)
      tracesDF[[fn.sel]] <- selected$traces

      # sync current index in the pool
      etr <- selected$ext.idx

      # record low ESS
      tmp.low <- ESS %>% as_tibble %>% add_column(fn, .before=1) %>% add_column(fn.sel, .before=2)
      lowESS <- lowESS %>% rbind(tmp.low)

    } else {
      tracesDF[[fn]] <- traces
    }
  }

  list()
}



#' @details
#' \code{summariseTrueValues} creates one file to contain all of the true values
#' from all LPhy simulations. The generated "true.log" can be used to BEAST 2
#' model validation pipeline \url{https://github.com/rbouckaert/DeveloperManual}.
#'
#'
#'
#' @param true.logs The list of one-line log files containing true values
#'                  from LPhy simulations, where one LPhy log file per simulation.
#' @param true.trees  The list of one-tree log files containing true trees
#'                    from LPhy simulations, where one LPhy tree file per simulation.
#' @param lphy.params  The list of
#' @param beast.params The list of
#' @keywords Coverage

#' @examples
#' WD = file.path("~/WorkSpace/MyValidations/")
#' setwd(WD)
#' true.logs = list.files(pattern = "_true.log")
#' true.trees = list.files(pattern = "_true_ψ.trees")
#' summariseTrueValues()
#'
#' @rdname Coverage
summariseTrueValues <- function(true.logs, true.trees=NA,
                                lphy.params=c("μ","Θ"),
                                beast.params=c("mu","Theta") ) {
  require(tidyverse)

  params = c("mu","Theta", "r_0", "r_1", "r_2",
             "kappa.1", "kappa.2", "kappa.3",
             "pi_0.A", "pi_0.C", "pi_0.G", "pi_0.T",
             "pi_1.A", "pi_1.C", "pi_1.G", "pi_1.T",
             "pi_2.A", "pi_2.C", "pi_2.G", "pi_2.T" )#,"psi.height")
  tre.params = c("total.br.len","tree.height")
  stats.name = c("mean", "HPD95.lower", "HPD95.upper", "ESS")
  # must have the same order of param
  params2 = c("μ","Θ","r_0","r_1","r_2","κ_0","κ_1","κ_2","π_0_0","π_0_1","π_0_2","π_0_3",
              "π_1_0","π_1_1","π_1_2","π_1_3","π_2_0","π_2_1","π_2_2","π_2_3")

}


