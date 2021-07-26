# Author: Walter Xie
# Accessed on 26 Jul 2021

#' @name TreeStats
#' @title Summarising tree statistics from the posterior
#'
#' @description
#' Summarise tree statistics from the posterior trees logged by
#' BEAST 2 \url{http://www.beast2.org}, such as the total branch length,
#' tree height, etc.
#'
#'
#' @details
#' \code{readTrees} return a data frame of  the posterior trees logged in *.trees file.
#'
#' @param file The file containing BEAST trees in the nexus format.
#'             The tree should look like "tree STATE_7500 = (1:7.260,(3:3.298 ... ))",
#'             where "tree STATE_7500" is the ID and "7500" is the state of MCMC.
#' @param burn.in the proportion of samples treated as burn-in stage from the MCMC.
#'                The calculation of burnin given the proportion is same to Tracer.
#' @param min.trees it will throw error, if the number of trees in the file
#'                  is smaller than this number, default to 10.
#' @keywords TreeStats
#' @export
#' @examples
#' tre.sta.df <- readTrees("data/RSV2long.trees")
#'
#' @rdname TreeStats
readTrees <- function(file, burn.in=0.1, min.trees=10) {
  require(tidyverse)
  require(ape)
  require(phytools)

  if (!file.exists(file)) stop("\nCannot find tree file ", file)

  tre.list <- read.nexus(file)
  n.tre <- length(tre.list)
  if (n.tre < min.trees) stop("Not enough trees in posterior ! ", tre.list)
  cat("Load ", n.tre, "trees from ", file, " having", Ntip(tre.list[[1]]), "tips ...\n")

  # same to tracer burnin method
  start = round(burn.in * n.tre) + 1
  tres <- tre.list[start:n.tre]
  #sum(tres[[1]]$edge.length)

  # 1801 elements
  tot.br.len.list <- mapply(function(x) sum(x$edge.length), tres)
  tre.height.list <- mapply(function(x) max(nodeHeights(x)), tres)

  stopifnot(length(tres) == length(tot.br.len.list) && length(tres) == length(tre.height.list))

  # remove STATE to retrive state numbers
  states <- names(tot.br.len.list) %>% str_remove_all("[A-Z]|[a-z]|\\_") %>% as.numeric
  cat(length(tres), "trees after burn-in", burn.in, ", start from state ",
      formatC(states[1], format="f", big.mark=",", digits=0), ".\n")

  # data frame
  tibble(states=states, total.br.len=(tot.br.len.list %>% unlist),
                    tree.height=(tre.height.list %>% unlist))
}


#' @details
#' \code{analyseTreeStats} get stats from data frame returned by
#' \code{readTrees}.
#' \code{id} is set to analyse the default selected columns.
#'
#' @param tree.stats The data frame.
#' @keywords TreeStats
#' @export
#' @examples
#' tre.sta <- analyseTreeStats(tre.sta.df)
#'
#' @rdname TreeStats
analyseTreeStats <- function(tree.stats, id=c("total.br.len","tree.height")) {
  require(TracerR)
  tre.sta <- TracerR::analyseTraces(tree.stats, id=id)
  # ? HPD95.lower.STATE_14150000
  tre.sta$trace <- sub("\\.STATE.*$","",tre.sta$trace, ignore.case = T)
  return(tre.sta)
}
