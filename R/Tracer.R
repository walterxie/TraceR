# Author: Rambaut A, Suchard MA, Xie D & Drummond AJ (2014) Tracer v1.6
# Modified by Walter Xie (D Xie)
# Accessed on 14 Jan 2021


#' @name Tracer
#' @title Tracer in R
#'
#' @description
#' A simplified implementation of Tracer \url{http://beast.bio.ed.ac.uk/tracer} in R.
#' It can read MCMC log files from BEAST 1 \url{http://beast.bio.ed.ac.uk},
#' BEAST 2 \url{http://www.beast2.org}, or MrBayes \url{http://beast.bio.ed.ac.uk}.
#'
#' It aims to provide various visualizations to understand the Bayesian MCMC result,
#' and also focus on making batch processing of multiple logs easier.
#'
#' @details
#' \code{readMCMCLog} reads MCMC log files and return a data frame
#' whose column names are parameters and row names are the number
#' of states at each sample.
#'
#' @param file The file to read/write.
#' @param samp.col The column to contain the state number of MCMC in the log file.
#' Default to the 1st column as the convention of BEAST.
#' @param rm.na.col If TRUE, then remove all columns with all
#' missing values (NA). Default to TRUE.
#' @param ... Other arguments passed to \code{\link{read_delim}}.
#' @keywords Tracer
#' @export
#' @examples
#' mcmc.log <- readMCMCLog("star.beast.log")
#'
#' @rdname Tracer
readMCMCLog <- function(file, delim="\t", comment = "#", samp.col=1, rm.na.col=TRUE, ...) {
  require(tidyverse)
  # rm comments #
  mcmc.log <- read_delim(file, delim, comment = comment, ...)
  #double from 2e-308 to 2e+308
  mcmc.log[[samp.col]] = as.double(mcmc.log[[samp.col]])

  # MCMC summary
  n.samples = nrow(mcmc.log)
  log.every = mcmc.log[[samp.col]][2]-mcmc.log[[samp.col]][1] # "Sample"
  last.state <- mcmc.log[[samp.col]][n.samples]

  cat("Loaded ", prettyNum(n.samples, big.mark=",",scientific=FALSE), " samples, ",
      "chain length ", prettyNum(last.state, big.mark=",",scientific=FALSE),
      ", log every ", prettyNum(log.every, big.mark=",",scientific=FALSE),
      ", file = ", file, ".\n", sep="")

  if (rm.na.col) {
    col.names <- colnames(mcmc.log)
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
#' \code{getTraces} preprocesses the data frame returned from
#' \code{readMCMCLog} to get the traces.
#' If \code{burn.in.state} is NULL, then use \code{burn.in}, else
#' consider \code{burn.in.state} first.
#'
#' @param mcmc.log The data frame from MCMC log file whose
#' column names are parameters and row names are the number
#' of states at each sample.
#' @keywords Tracer
#' @export
#' @examples
#' mcmc.traces <- getTraces(mcmc.log, burn.in=0.1)
#'
#' @rdname Tracer
getTraces <- function(mcmc.log, burn.in.state=NULL, burn.in=0.1, samp.col=1, verbose=TRUE) {
  # states in rownames, double from 2e-308 to 2e+308
  samples <- mcmc.log[[samp.col]]
  #row.names(mcmc.log) <- samples # Setting row names on a tibble is deprecated
  n.sample <- length(samples)
  last.state <- samples[n.sample]
  step.size <- samples[2]-samples[1]

  if (!is.null(burn.in.state)) {
    # e.g. burn.in.state = 10000
    traces <- mcmc.log[mcmc.log[[samp.col]] > burn.in.state, -1]
    burn.in <- burn.in.state / last.state
  } else {
    # exclude the state == burn.in * n.sample, making same as Tracer
    start = round(burn.in * n.sample) + 1
    traces <- mcmc.log[start:n.sample, -1]
    burn.in.state <- burn.in * last.state
  }

  if (verbose)
    cat("Remaining ", nrow(traces), " samples after burn-in ", burn.in,
        ", where burn-in ends at state ", burn.in.state,
        ", step size ", step.size, ", last state ", last.state, ".\n", sep="");

  attr(traces,"burn.in") <- burn.in
  attr(traces,"step.size") <- step.size
  attr(traces,"last.state") <- last.state
  return(traces)
}


#' @details
#' \code{analyseTraces} get stats from data frame returned by
#' \code{getTraces}.
#' If \code{id} is empty, then analyse all columns, otherwise
#' analyse the selected column(s) only.
#'
#' @param traces The data frame.
#' @keywords Tracer
#' @export
#' @examples
#' stats <- analyseTraces(traces)
#'
#' @rdname Tracer
analyseTraces <- function(traces, id=c(),
                          trace.stats = c("mean", "stderr of mean", "stdev", "variance", "median",
                                          "mode", "geometric mean", "95% HPD Interval",
                                          "effective sample size (ESS)"), #"auto-correlation time (ACT)",
                          verbose=TRUE) {
  if (length(id) > 0) {
    traces <- traces %>% select(id)
  }
  attrs <- attributes(traces)
  if ("step.size" %in% names(attrs) && "last.state" %in% names(attrs)) {
    stats <- lapply(traces, function(x)
      analyse(x, chain.length=as.double(attrs$last.state),
              log.every=as.integer(attrs$step.size), verbose=verbose))
  } else {
    stats <- lapply(traces, function(x) analyse(x, verbose))
  }

  # Row names aren't preserved with many dplyr and tidyr operations.
  # It's best to not rely on them but use rownames_to_column() instead.
  stats <- stats %>% lapply(unlist) %>% data.frame %>%
    rownames_to_column("trace") %>% as_tibble
  return(stats)
}

#' @details
#' \code{analyse} calculates statistics of one parameter,
#' such as \code{ESS} and \code{standard error of mean}.
#'
#' @param x The trace of one parameter.
#' @keywords Tracer
#' @export
#' @examples
#' stats.ls <- analyse(traces[["likelihood"]])
#'
#' @rdname Tracer
analyse <- function(trace, chain.length=NA, log.every=NA, verbose=TRUE) {
  trace.corr = analyseCorrelation(trace)

  m <- trace.corr$mean
  sem <- trace.corr$stderr.of.mean
  ess <- trace.corr$ESS
  sd <- sd(trace)
  v <- var(trace)
  md <- median(trace)
  mo <- mode(trace)
  gm <- gmMean(trace)
  n.x <- length(trace)

  suppressMessages(require(TeachingDemos))
  hpd95 <- emp.hpd(trace, conf=0.95)

  # attrs not passed by lapply
  act <- NA
  if (!is.na(log.every)) {
    # ESS = (log.every * n.x) / ACT
    act <- (as.integer(log.every) * n.x) / ess;
  }

  list(
    mean = m,
    stderr.of.mean = sem,
    stdev = sd,
    variance = v,
    median = md,
    mode = mo, geometric.mean = gm,
    HPD95.lower = hpd95[1], HPD95.upper = hpd95[2],
    ACT = act,
    ESS = ess,
    n.samples = n.x,
    chain.length = chain.length,
    log.every = log.every
  )
}

### internal functions ###

# the code is directly imported from Tracer
# log.every (step size) is the size of the step between states, also sampling frequency
analyseCorrelation <- function(x, MAX.LAG=2000, verbose=FALSE) {
  n.x = length(x);
  max.lag = min(n.x-1, MAX.LAG); # Note: lag index starts from 1 NOT 0
  samp.mean = mean(x)

  if (verbose) cat("Input", n.x, "samples, ")

  gamma.stat=c()
  for (lag in 1:max.lag) {
    gamma.stat[lag] = 0
    # gamma.stat[1] = SUM_{i=1}( (x_i - mean)(x_{i+1} - mean) ) / N
    # gamma.stat[lag] = SUM_{i=1}( (x_i - mean)(x_{i+lag-1} - mean) ) / N
    for (s in 1:(n.x - lag + 1)) {
      del1 = x[s] - samp.mean;
      # Note: both samples and lag index start from 1 NOT 0
      del2 = x[s + lag - 1] - samp.mean;
      gamma.stat[lag] = gamma.stat[lag] + (del1 * del2);
    } # end s loop
    # Note: lag index starts from 1 NOT 0
    gamma.stat[lag] = gamma.stat[lag] / (n.x - lag - 1);

    if (lag == 1) {
      var.stats = gamma.stat[1];
    } else if (lag %% 2 != 0) { # Note: lag index starts from 1 NOT 0
      # fancy stopping criterion :)
      if (gamma.stat[lag - 1] + gamma.stat[lag] > 0) {
        var.stats = var.stats + 2.0 * (gamma.stat[lag - 1] + gamma.stat[lag]);
      } else { # stop for loop
        gamma.stat = gamma.stat[1:lag]
        break
      }
    }
  }

  # standard error of mean
  stderr.of.mean = sqrt(var.stats / n.x);

  ### avoid to provide log.every
  # auto correlation time
  # if (gamma.stat[1] == 0) {
  #   ACT = 0;
  # } else {
  #   ACT = log.every * var.stats / gamma.stat[1];
  # }

  ### effective sample size
  # ESS = N * gamma.stat[1] / ( gamma.stat[1] + SUM_{lag=3}(2.0 * (gamma.stat[lag-1] + gamma.stat[lag])) )
  if (gamma.stat[1] == 0) {
    ESS = 1;
  } else {
    #ESS = (log.every * n.x) / ACT;
    ESS = n.x * gamma.stat[1] / var.stats;
  }

  if (verbose) cat("ESS is ", ESS, ".\n", sep="")

  list(ESS=ESS, mean=samp.mean, stderr.of.mean=stderr.of.mean, n.samples = n.x,
       var.stats=var.stats, gamma.stat=gamma.stat) #
}

# http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gmMean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
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


