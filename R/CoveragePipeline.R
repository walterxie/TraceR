# Author: Walter Xie
# Accessed on 11 Aug 2021

#' @name CovgPip
#' @title Coverage Pipeline from Model Validation
#'
#' @description
#' The pipeline to generate a coverage report of the model validation.
#' The input files for this pipeline are:
#' BEAST logs (tree logs) containing samples from the posterior,
#' and LPhy logs (tree logs) containing true values.
#' They suppose to be in the same working directory.
#'
#' The steps of this pipeline are listed in the details:
#'
#'
#' @details
#' Step 1: \code{pipCreateSimulationSummaries} summarises traces statistics for
#' every simulations (normally starting 110, where 10 extra simulations
#' are used for replacing any low ESS results in the next step),
#' and save each result into a tsv file.
#'
#' It looks like:
#' trace           posterior  likelihood  ...
#' mean            -16848.58  -16269.03   ...
#' stderr.of.mean	0.249      0.200       ...
#' ...
#'
#' Note: if tree stats are logged into BEAST log file,
#' you can set \code{tree.file=NA} to skip creating tree stats tsv file
#' from BEAST tree logs.
#'
#' It includes \code{\link{summariseTracesAndTrees}}.
#'
#' @param log.files The vector of BEAST log file names.
#' @param burn.in   The proportion of samples treated as the burn-in stage
#'                  during MCMC. Default to 0.1.
#' @param process.tree.file Default to FALSE, to skip creating tree stats tsv file
#'                  from BEAST tree logs.
#' @keywords Coverage
#' @export
#' @examples
#' WD = file.path("~/WorkSpace/MyValidations/")
#' setwd(WD)
#' # 110 simulations
#' log.files = list.files(pattern = "_([0-9]+).log")
#' # default to skip trees, because their stats are logged in the log files
#' beast.params <- pipCreateSimulationSummaries(log.files, burn.in=0.1)
#'
#' @rdname CovgPip
pipCreateSimulationSummaries <- function(log.files, burn.in=0.1,
                                         process.tree.file=FALSE) {
  require(tidyverse)
  cat("\nPipeline step 1: summarising traces statistics ...\n")

  # output x is BEAST log file name
  stats.fn.fun=function(x){ paste0(sub('\\.log$', '', x), ".tsv") }
  tree.stats.fn.fun=function(x){ paste0(x, ".tsv") }

  for(lg in log.files) {
    # assume same file stem
    if (process.tree.file)
      tree.file=paste0(sub('\\.log$', '', lg), ".trees")
    else
      tree.file=NA

    res <- summariseTracesAndTrees(lg, tree.file=tree.file, burn.in=burn.in,
                                   stats.fn.fun=stats.fn.fun,
                                   tree.stats.fn.fun=tree.stats.fn.fun)
  }
  # parameters in BEAST log
  names(res$stats)
}


#' @details
#' Step 2: \code{pipSelectValidResults} selects 100 simulations
#' (run 110 in total), where the ESS of every parameters are guaranteed >= 200.
#' If not, then replace it to the one from the extra 10 sequentially,
#' and check ESS.
#'
#' If all extra 10 are used but there exists any low-ESS simulations,
#' then the pipeline stops and inform to re-run all simulations with longer
#' MCMC chain length.
#'
#' It includes \code{\link{selectResultByESS}}.
#'
#' @param i.sta,i.end  The start index or end index of batch runs,
#'                     also used in the file names.
#' @param prefix  The prefix of statistics summary files.
#' @param extra.file.fun,extra.tree.file.fun The function to list
#'                the statistics summary file names of (10) extra runs.
#'                Set extra.tree.file.fun = NA, if you do not need to
#'                handle BEAST tree log files.
#' @keywords Coverage
#' @export
#' @examples
#' all.stats = list.files(pattern = "_([0-9]+).tsv")
#' sele.list <- pipSelectValidResults(i.sta=0, i.end=99, prefix="al2",
#'                                    extra.tree.file.fun=NA)
#'
#' @rdname CovgPip
pipSelectValidResults <- function(i.sta=0, i.end=99, prefix="sim",
                                  extra.file.fun=function(ext){ list.files(pattern=paste0("_10([0-9])", ".", ext)) },
                                  extra.tree.file.fun=function(ext){ list.files(pattern=paste0("_10([0-9])", ".", ext)) } ) {
  require(tidyverse)
  n.sim <- i.end-i.sta+1
  cat("\nPipeline step 2: selecting ", n.sim, " simulation results (ESS >= 200) ...\n")

  file.steam.fun=function(prefix, i){ paste0(prefix, "_", i) }
  extra.file.steam.fun=function(f){ sub('\\.tsv$', '', f) }

  file.postfix="tsv"
  tree.file.postfix="trees.tsv"
  # extra 10 simulations: *_100.tsv ~ *_109.tsv
  extra.files = extra.file.fun(file.postfix)
  cat("Prepare ", length(extra.files), " extra stats files : ",
      paste(extra.files, collapse = T), ".\n")

  extra.tree.files = c()
  if (is.function(extra.tree.file.fun)) {
    extra.tree.files = extra.tree.file.fun(tree.file.postfix)
    cat("Prepare ", length(extra.tree.files), " extra tree stats files : ",
        paste(extra.tree.files, collapse = T), ".\n")
  } else
    tree.file.postfix=NA # Set tree.file.postfix=NA and keep extra.tree.files=c()

  sele.res <- selectResultByESS(i.sta=i.sta, i.end=i.end, prefix=prefix,
                                file.steam.fun=file.steam.fun,
                                file.postfix=file.postfix, tree.file.postfix=tree.file.postfix,
                                extra.file.steam.fun=extra.file.steam.fun,
                                extra.files = extra.files, extra.tree.files=extra.tree.files)

  if (nrow(sele.res$lowESS)>0)
    write_tsv(lowESS, file.path(paste0("low-ESS.tsv")) )

  stopifnot(length(sele.res$selected) == n.sim)
  return(sele.res$selected)
}

#' @details
#' Step 3: \code{pipCreateParameterSummaries} summarise BEAST results for
#' each of parameters.It includes \code{\link{summariseParameters}}.
#'
#' @param selected   The list of data frames containing traces summary,
#'                   produced by \code{\link{pipSelectValidResults}}.
#' @param params     The vector of parameter names in the BEAST log files.
#' @param stats.name The vector of names of statistics. They have to be one of
#'                   names from the 1st column 'trace' in *.tsv file.
#' @keywords Coverage
#' @export
#' @examples
#' beast.params = c("mu","Theta", "r_0", "r_1", "r_2", "psi.treeLength", "psi.height")
#' summ <- summariseParameters(sele.list, params = beast.params)
#' summ$param.summaries[["mu"]]
#'
#' @rdname CovgPip
pipCreateParameterSummaries <- function(selected=list(),
                                        params = c("mu","Theta", "psi.treeLength", "psi.height"),
                                        stats.name = c("mean", "HPD95.lower", "HPD95.upper", "stdev", "ESS") ) {
  require(tidyverse)
  cat("\nPipeline step 3:\n")

  summ <- summariseParameters(selected, params = params, stats.name = stats.name,
                              output.file.fun=function(x){ paste0(x, ".tsv") } )

  for (pa in params) {
    df <- summ$param.summaries[[pa]]
    write_tsv(df, paste0(pa, ".tsv"))
    cat("Write ", paste0(pa, ".tsv"), "\n")
  }
  cat("min ESS of selected results for each parameter = ", paste(summ$minESS, collapse = ", "), "\n")
  cat("min ESS of all = ", min(summ$minESS), "\n")
  return(summ)
}

#' @details
#' Step 4: \code{pipCreateTrueValueSummaries} creates one final summary file
#' containing true values for every parameters.
#' It includes \code{\link{summariseTrueValues}}.
#'
#' @param selected.fn.steam The vector of file name steams selected by
#'                          \code{pipSelectValidResults} to summarise
#'                          the coverage.
#' @param params  The vector of parameter names in LPhy.
#' @param add.tree.stats Default to TRUE, to summarise true tree stats
#'                       which are fixed to "total.br.len" (total branch length)
#'                       and "tree.height" at the moment.
#'                       Set to False, if tree stats is not required
#'                       in the data frame.
#' @param log.file.fun  The function to get one-line log file name,
#'                      containing true values from LPhy simulations,
#'                      where one LPhy log file per simulation.
#' @param tree.file.fun The function to get one-tree log file name,
#'                      containing true trees from LPhy simulations,
#'                      where one LPhy tree file per simulation.
#' @keywords Coverage
#' @export
#' @examples
#' # list.files(pattern = "_true.log")
#' df.tru <- pipCreateTrueValueSummaries(names(sele.list), params=c("μ","Θ"), add.tree.stats=TRUE)
#'
#' @rdname CovgPip
pipCreateTrueValueSummaries <- function(selected.fn.steam=c(),
                                        params=c("μ","Θ", "r_0", "r_1", "r_2"),
                                        add.tree.stats=TRUE,
                                        log.file.fun=function(x){ paste0(x,"_true.log") },
                                        tree.file.fun=function(x){ paste0(x,"_true_ψ.trees") } ) {
  require(tidyverse)
  cat("\nPipeline step 4:\n")

  # LPhy parameters in the true-value log
  tru.log <- log.file.fun(selected.fn.steam[[1]])
  pa.nm <- read_tsv(tru.log) %>% names
  cat("Find ", length(pa.nm), " LPhy parameters in ", tru.log, " : ",
      paste(pa.nm, collapse = T), "\n")
  # list.files(pattern = "_true_ψ.trees")

  # the order of parameters same to BEAST parameters
  df.tru <- summariseTrueValues(selected.fn.steam, params=params,
                                add.tree.stats=TRUE,
                                log.file.fun = log.file.fun, tree.file.fun = tree.file.fun )
  write_tsv(df.tru, "trueValue.tsv")
  return(df.tru)
}

#' @details
#' Step 5: \code{reportCoverages} marks how many true values are falling
#' into or outside the 95% HPD interval of posteriors for each parameter,
#' and report the overall coverage.
#' It includes \code{\link{markInOut}}.
#'
#' Note: the same parameter may be given different names between LPhy script
#' and BEAST XML/log, please ensure that you match them correctly.
#'
#' @param beast.params,lphy.params The vector of BEAST or LPhy parameters,
#'                                 and they have to exactly match each other.
#' @param beast.summ.file.fun The function to create beast posterior summary file name.
#' @param true.val.file The file containing true values for every parameters.
#' @keywords Coverage
#' @export
#' @examples
#' covg <- reportCoverages(beast.params = c("mu","Theta", "psi.treeLength", "psi.height"),
#'                         lphy.params = c("μ","Θ", "total.br.len","tree.height"))
#'
#' @rdname CovgPip
reportCoverages <- function(beast.params = c("mu","Theta", "psi.treeLength", "psi.height"),
                            lphy.params = c("μ","Θ", "total.br.len","tree.height"),
                            beast.summ.file.fun=function(x){ paste0(x,".tsv") },
                            true.val.file="trueValue.tsv") {
  require(tidyverse)
  cat("\nPipeline step 5: report the overall coverage ...\n")
  stopifnot(length(beast.params)==length(lphy.params) && length(beast.params)>0)

  # true
  stopifnot(file.exists(true.val.file))
  df.tru <- read_tsv(true.val.file, col_types = cols())

  covg <- c()

  for (i in 1:length(beast.params)) {
    pos.fn <- beast.summ.file.fun(beast.params[i])
    stopifnot(file.exists(pos.fn))

    df.pos <- read_tsv(pos.fn, col_types = cols())
    inOut <- markInOut(df.pos, df.tru, tru.val.par=lphy.params[i])

    write_tsv(inOut, paste0(beast.params[i], "-coverage.tsv"))

    covg <- c(covg, nrow(inOut[inOut$is.in==T,])/nrow(inOut))
  }

  return(tibble(lphy=lphy.params, beast=beast.params, covg=covg))
}



