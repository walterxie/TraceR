# Author: Walter Xie
# Accessed on 6 Aug 2021

#' @name TracerAdv
#' @title The Pipeline Functions to Summarise Bayesian MCMC Results
#'
#' @description
#' The several integrated functions based on \code{\link{Tracer}} and
#' \code{\link{TreeStats}}, which are useful to create a pipeline to
#' summarise groups of batch runs.
#'
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
#' @keywords TracerAdv
#' @export
#' @examples
#' WD = file.path("~/WorkSpace/MyValidations/")
#' setwd(WD)
#' log.files = list.files(pattern = "_([0-9]+).log")
#' log.files # exclude _true.log
#' for(lg in log.files) {
#'   # assume same file stem
#'   tree.file=paste0(sub('\\.log$', '', lg), ".trees")
#'   res <- summariseTracesAndTrees(lg, tree.file)
#' }
#'
#' # skip tree logs
#' summariseTracesAndTrees(log.file="sim_0.log", tree.file=NA)
#'
#' @rdname TracerAdv
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

  if (is.function(stats.fn.fun))
    write_tsv(stats, stats.fn.fun(log.file))

  # add tree stats
  tre.sta <- NULL
  if (!is.na(tree.file)) {
    if (!file.exists(tree.file)) stop("Cannot find file : ", tree.file)

    tre.sta.df <- readTrees(tree.file, burn.in=burn.in)
    tre.sta <- analyseTreeStats(tre.sta.df)

    # ? HPD95.lower.STATE_14150000
    #tre.stats$trace <- sub("\\.STATE.*$","",tre.stats$trace, ignore.case = T)
    if (is.function(tree.stats.fn.fun))
      write_tsv(tre.sta, tree.stats.fn.fun(tree.file))
  }
  list(stats=stats, tree.stats=tre.sta)
}

#' @details
#' \code{readTracesTSV} loads traces tsv format file into a data frame,
#' and report the minimum ESS.
#' @seealso \code{\link{read_tsv}}
#'
#' @param traces.file  The tsv file containing traces statistics.
#' @keywords TracerAdv
#' @export
#' @examples
#' res <- readTracesTSV(traces.file="sim_0.tsv")
#' # res$minESS
#' traces <- res$traces
#'
#' @rdname TracerAdv
readTracesTSV <- function(traces.file="sim_0.tsv", ...) {
  require(tidyverse)
  stopifnot(file.exists(traces.file))
  # "trace" col has "mean", ..., "HPD95.lower", "HPD95.upper", "ESS"
  # suppress messages from read_tsv
  traces <- try(read_tsv(traces.file, col_types = cols(), ...))
  # take ESS row
  ESS <- traces %>% filter(trace == "ESS") %>% select(!trace)
  minESS <- min(ESS %>% as.numeric)
  cat(traces.file, " has min ESS ", minESS, ".\n")
  list(traces=traces, minESS=minESS)
}


