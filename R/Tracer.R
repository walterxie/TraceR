# Author: Rambaut A, Suchard MA, Xie D & Drummond AJ (2014) Tracer v1.6
# Modified by Walter Xie (D Xie)
# Accessed on 18 May 2016


#' @name Tracer
#' @title Tracer in R
#'
#' @description
#' A simplified implementation of Tracer \url{http://beast.bio.ed.ac.uk/tracer} in R.
#' It can read MCMC log files from BEAST 1 \url{http://beast.bio.ed.ac.uk},
#' BEAST 2 \url{http://www.beast2.org}, or MrBayes \url{http://beast.bio.ed.ac.uk}.
#'
#' It aims to provide various visualizations to understand the Bayesian MCMC result,
#' and alos focus on making batch processing of multipe logs easier.
#'
#' @details
#' \code{readMCMCLog} reads MCMC log files and return a data frame
#' whose column names are parameters and row names are the number
#' of states at each sample.
#'
#' @param file The file to read/write.
#' @param rm.na.col If TRUE, then remove all columns with all
#' missing values (NA). Default to TRUE.
#' @param ... Other arguments passed to \code{\link{readFile}}.
#' @keywords Tracer
#' @export
#' @examples
#' mcmc.log <- readMCMCLog("data/star.beast.log")
#'
#' @rdname Tracer
readMCMCLog <- function(file, rm.na.col=TRUE, ...) {
  require(ComMA)
  mcmc.log <- ComMA::readFile(file, comment.char = "#", msg.file="MCMC log",
                              msg.col="parameters", msg.row="samples", ...)

  if (rm.na.col) {
    col.names <- ComMA::trimAll(names(mcmc.log))
    # exclude all columns with all NA, such as BEAST (< 2.4.1) log
    no.empty.col <- col.names[sapply( mcmc.log, function(x) !all(is.na(x)) )]
    if (length(col.names) != length(no.empty.col)) {
      mcmc.log <- mcmc.log[,no.empty.col]
      warning("Remove ", length(col.names) - length(no.empty.col), " empty column(s) !\n",
              "There are ", length(no.empty.col), " parameters for analysis !\n")
    }
  }

  attr(mcmc.log,"file") <- file
  return(mcmc.log)
}

#' @details
#' \code{getBurnIn} returns a list of statistics inlcuding \code{burn.in}
#' and \code{step.size}, given the data frame read from \code{readMCMCLog}.
#' It also calculate \code{burn.in} from \code{burn.in.perc}, if it is NULL.
#'
#' @param burn.in The number of states of burn in. Default to NULL
#' to use \code{burn.in.perc}. Otherwise, use \code{burn.in}.
#' @param burn.in.perc The percentage of total chain length to be burn in.
#' Default to 0.1. Only used if \code{burn.in} is NULL.
#' @keywords Tracer
#' @export
#' @examples
#' mcmc.traces <- getTraces(mcmc.log)
#'
#' @rdname Tracer
getBurnIn <- function(samples, burn.in=NULL, burn.in.perc=0.1,
                      verbose=TRUE) {
  samples <- as.double(samples)
  n.sample <- length(samples)
  last.state <- max(samples)
  step.size <- samples[2] - samples[1]

  if (is.null(burn.in))
    burn.in <- last.state * burn.in.perc

  if (verbose)
    cat("Set burn in =", burn.in, "for", n.sample, "samples.\n",
        "step.size =", step.size, ", last.state =", last.state, "\n")

  list(
    n.sample = n.sample,
    last.state = last.state,
    step.size = step.size,
    burn.in = burn.in
  )
}

#' @details
#' \code{getTraces} preprocesses the data frame returned from
#' \code{readMCMCLog}.
#' If \code{burn.in} is NULL, then use \code{burn.in.perc}, else
#' consider \code{burn.in} first.
#'
#' @param mcmc.log The data frame from MCMC log file whose
#' column names are parameters and row names are the number
#' of states at each sample.
#' @keywords Tracer
#' @export
#' @examples
#' mcmc.traces <- getTraces(mcmc.log)
#'
#' @rdname Tracer
getTraces <- function(mcmc.log, burn.in=NULL, burn.in.perc=0.1,
                      verbose=TRUE) {
  # states in rownames, double from 2e-308 to 2e+308
  samples <- as.double(rownames(mcmc.log))

  burn.in.stats <- getBurnIn(samples, burn.in=burn.in,
                            burn.in.perc=burn.in.perc, verbose=verbose)
  burn.in <- burn.in.stats$burn.in
  n.sample <- burn.in.stats$n.sample
  last.state <- burn.in.stats$last.state
  step.size <- burn.in.stats$step.size

  burn.in.ind <- match(burn.in, samples)
  if (is.na(burn.in.ind))
    stop("Incorrect burn in", burn.in, ", cannot find it from samples !")

  # assume samples are ordered in log
  traces <- mcmc.log[burn.in.ind:n.sample,]

  if (verbose)
    cat("Remaining", nrow(mcmc.log), "samples after burn-in.\n")

  attr(traces,"burn.in") <- burn.in
  attr(traces,"step.size") <- step.size
  attr(traces,"last.state") <- last.state
  return(traces)
}

trace.stats <- c("mean", "stderr of mean", "stdev", "variance", "median",
                 "mode", "geometric mean", "95% HPD Interval",
                 "auto-correlation time (ACT)", "effective sample size (ESS)")


analyseTraces <- function(traces, id=NULL, verbose=TRUE) {
  attr.traces <- attributes(traces)
  step.size <- attr.traces$step.size

  if (! is.null(id)) {
    summary.stats <- apply(traces[,id], 2, function(x) analyse(x, step.size, verbose))
  } else {
    summary.stats <- apply(traces, 2, function(x) analyse(x, step.size, verbose))
  }

  return(summary.stats)
}

analyse <- function(trace, step.size, verbose=TRUE) {
  m <- mean(trace)
  sd <- sd(trace)
  sem <- std.err(trace)
  v <- var(trace)
  md <- median(trace)
  mo <- mode(trace)
  gm <- gm_mean(trace)

  suppressMessages(require(TeachingDemos))
  hpd95 <- emp.hpd(trace, conf=0.95)
  act <- auto.corr.time(trace)
  ess <- effective.sample.size(act, step.size, length(trace))

  list(
    mean = m,
    stdev = sd,
    stderr.of.mean = sem,
    variance = v,
    median = md,
    mode = mo,
    geometric.mean = gm,
    hpd95 = hpd95,
    ACT = act,
    ESS = ess
  )
}

# The standard error is just the standard deviation divided by the square root of the sample size.
# http://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
std.err <- function(x) sd(x)/sqrt(length(x))

# http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm.mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

auto.corr.time <- function(x){

}

effective.sample.size <- function(act, step.size, n.sample){
  return(step.size * n.sample / act)
}




