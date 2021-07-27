# Author: Walter Xie
# Accessed on 27 Jul 2021

#' @name Coverage
#' @title Coverage Report from Model Validation
#'
#' @description
#' The coverage report of a model validation result.
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
#' @keywords Coverage
#' @export
#' @examples
#' WD = file.path("~/WorkSpace/MyValidations/")
#' setwd(WD)
#' log.files = list.files(pattern = ".log")
#' for(lg in log.files) {
#'   tree.file=paste0(sub('\\.log$', '', lg), ".trees")
#'   summariseTracesAndTrees(lg, tree.file)
#' }
#'
#' @rdname Coverage
summariseTracesAndTrees <- function(log.file, tree.file=NA, burn.in=0.1) {
  require(tidyverse)

  cat("\nProcess ", log.file, "...\n")
  if (!file.exists(log.file)) stop("\nRequire log file ", log.file)

  # read MCMC log
  mcmc.log <- readMCMCLog(log.file)
  # get traces and remove burn in
  traces <- getTraces(mcmc.log, burn.in=burn.in)
  # get stats
  stats <- analyseTraces(traces)

  write_tsv(stats, paste0(sub('\\.log$', '', log.file), ".tsv"))

  # add tree stats
  if (!is.na(tree.file)) {

    if (!file.exists(tree.file)) stop("\nRequire tree file ", tree.file)

    tre.sta.df <- readTrees(tree.file)
    tre.sta <- analyseTreeStats(tre.sta.df)

    # ? HPD95.lower.STATE_14150000
    #tre.stats$trace <- sub("\\.STATE.*$","",tre.stats$trace, ignore.case = T)

    write_tsv(tre.sta, paste0(tree.file, ".tsv"))
  }

}




