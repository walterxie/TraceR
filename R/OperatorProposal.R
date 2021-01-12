# Author: Walter Xie
# Accessed on 29 May 2020

#' @name Operator
#' @title Summarising proposals for the operator
#'
#' @description
#' Extract the summary of operator proposals from BEAST 2 \url{http://www.beast2.org} output.
#'
#' Tuning: The value of the operator's tuning parameter, or '-' if the operator can't be optimized.
#' accept: The total number of times a proposal by this operator has been accepted.
#' reject: The total number of times a proposal by this operator has been rejected.
#' Pr(m): The probability this operator is chosen in a step of the MCMC (i.e. the normalized weight).
#' Pr(acc|m): The acceptance probability (\code{accept} as a fraction of the total proposals for this operator).
#'
#'
#' @details
#' \code{readState} reads MCMC state log file and return a data frame
#' to summarise operators' acceptence and rejection.
#'
#' @param file The file to read/write.
#' @param acc.cols,id.cols The column names which matches the names in json.
#' @keywords Operator
#' @export
#' @examples
#' ops <- readState("data/star.beast.state"")
#'
#' @rdname Operator
readState <- function(file, acc.cols=c("accept","reject"), id.cols=c("id","p")) {
  require(jsonlite)
  require(tidyverse)

  con  <- file(file, open = "r")
  lins = readLines(con)
  close(con)

  # parse after {"operators":[
  json = lins[grepl("\\{\"id\"\\:", lins)] %>% paste(collapse = " ")  %>% paste("[",.,"]")
  ops = fromJSON(json)
  # Operator Tuning #accept #reject Pr(m) Pr(acc|m)
  ops = ops %>% select(!!id.cols,!!acc.cols) %>%
    mutate(acc.rej=rowSums(select(.,!!acc.cols))) %>%
    mutate(Pr.m=acc.rej/sum(select(.,!!acc.cols))) %>%
    mutate(Pr.acc.m=(!!as.name(acc.cols[1]))/acc.rej) %>%
    select(-acc.rej)
  ops
}
