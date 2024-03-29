% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tracer.R
\name{Tracer}
\alias{Tracer}
\alias{readMCMCLog}
\alias{getTraces}
\alias{analyseTraces}
\alias{analyse}
\title{Tracer in R}
\usage{
readMCMCLog(
  file,
  delim = "\\t",
  comment = "#",
  col_types = cols(),
  samp.col = 1,
  rm.na.col = TRUE,
  ...
)

getTraces(
  mcmc.log,
  burn.in.state = NULL,
  burn.in = 0.1,
  samp.col = 1,
  verbose = TRUE
)

analyseTraces(
  traces,
  id = c(),
  trace.stats = c("mean", "stderr of mean", "stdev", "variance", "median", "mode",
    "geometric mean", "95\% HPD Interval", "effective sample size (ESS)"),
  verbose = TRUE
)

analyse(trace, chain.length = NA, log.every = NA, verbose = TRUE)
}
\arguments{
\item{file}{The file to read/write.}

\item{samp.col}{The column to contain the state number of MCMC in the log file.
Default to the 1st column as the convention of BEAST.}

\item{rm.na.col}{If TRUE, then remove all columns with all
missing values (NA). Default to TRUE.}

\item{...}{Other arguments passed to \code{\link{read_delim}}.}

\item{mcmc.log}{The data frame from MCMC log file whose
column names are parameters and row names are the number
of states at each sample.}

\item{burn.in}{the proportion of samples treated as the burn-in stage
during MCMC. Default to 0.1.}

\item{traces}{The data frame.}

\item{x}{The trace of one parameter.}
}
\description{
A simplified implementation of Tracer \url{http://beast.bio.ed.ac.uk/tracer} in R.
It can read MCMC log files from BEAST 1 \url{http://beast.bio.ed.ac.uk},
BEAST 2 \url{http://www.beast2.org}, or MrBayes \url{http://beast.bio.ed.ac.uk}.

It aims to provide various visualizations to understand the Bayesian MCMC result,
and also focus on making batch processing of multiple logs easier.
}
\details{
\code{readMCMCLog} reads MCMC log files and return a data frame
whose column names are parameters and row names are the number
of states at each sample.
\code{col_types = cols()} is used to suppress messages from
\code{\link{read_delim}}.

\code{getTraces} preprocesses the data frame returned from
\code{readMCMCLog} to get the traces.
If \code{burn.in.state} is NULL, then use \code{burn.in}, else
consider \code{burn.in.state} first.

\code{analyseTraces} get stats from data frame returned by
\code{getTraces}.
If \code{id} is empty, then analyse all columns, otherwise
analyse the selected column(s) only.

\code{analyse} calculates statistics of one parameter,
such as \code{ESS} and \code{standard error of mean}.
}
\examples{
mcmc.log <- readMCMCLog("star.beast.log")

mcmc.traces <- getTraces(mcmc.log, burn.in=0.1)

stats <- analyseTraces(traces)

stats.ls <- analyse(traces[["likelihood"]])

}
\keyword{Tracer}
