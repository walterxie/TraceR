# Author: Rambaut A, Suchard MA, Xie D & Drummond AJ (2014) Tracer v1.6
# Modified by Walter Xie (D Xie)
# Accessed on 18 May 2016

#' @name TracerGraph
#' @title Tracer Graph in R
#'
#' @description
#' The visualization functions to understand the Bayesian MCMC result.
#'
#' @details
#' \code{gtTraces} returns a \code{\link{gtable}} object to plot
#' the trace given MCMC log data frame,
#' whose column names are parameters and row names are the number
#' of states at each sample.
#'
#' @param mcmc.log The data frame from MCMC log file whose
#' column names are parameters and row names are the number
#' of states at each sample.
#' @param trace.id A vector of name(s) of trace(s) to plot.
#' Default to "posterior".
#' @param drop.1st.row If TRUE, also the default, then remove
#' the 1st row of MCMC log.
#' @param line.or.point If 1, the default, then display the line only;
#' if 2, then only points; if 3, then both the line and points.
#' @param point.size The size of points. Default to 1.
#' @param ... Other arguments passed to \code{\link{readFile}}.
#' @keywords Tracer
#' @export
#' @examples
#' mcmc.log <- readMCMCLog("data/star.beast.log")
#' burn.in.stats <- getBurnIn(rownames(mcmc.log))
#' burn.in.stats
#' names(mcmc.log)
#' gt <- gtTraces(mcmc.log, burn.in.stats$burn.in)
#' plot(gt)
#'
#' gt <- gtTraces(mcmc.log, burn.in.stats$burn.in, trace.id=c("TreeHeight.Species", "TreeHeight.t:tree_0_0", "TreeHeight.t:tree_0_1"))
#'
#' @rdname TracerGraph
gtTraces <- function(mcmc.log, burn.in, trace.id=c("posterior"), drop.1st.row=TRUE,
                     line.or.point=1, point.size=1, point.alpha=1, line.alpha=0.8,
                     title="Trace", x.lab=NULL, y.lab=NULL, legend.title="parameters",
                     verbose=TRUE, ...) {
  if (length(trace.id) < 1)
    stop("Give at least 1 parameter to plot trace, where the 1st should be 'state' !\n",
         paste(trace.id, collapse = ","))
  if (!all(is.element(trace.id, colnames(mcmc.log))))
    stop("Cannnot find '", paste(trace.id, collapse = ","), "' columns to plot !")

  if (drop.1st.row)
    mcmc.log <- mcmc.log[-1,]

  suppressMessages(require(reshape2))
  if (verbose)
    cat("Traces to plot: ", paste(trace.id, collapse = ","), "\n")
  mcmc.log$state <- as.double(rownames(mcmc.log))
  melt.df <- melt(mcmc.log[,c("state",trace.id)], id="state")

  if (length(trace.id) > 1) {
    # multipe parameters
    if (!is.null(legend.title) && length(legend.title) > 1) {
      # labs(fill=legend.title) not work for colour=?, rename "variable"
      colnames(melt.df)[colnames(melt.df) == "variable"] <- legend.title
      group.id=legend.title
      colour.id=legend.title
    } else {
      colnames(melt.df)[colnames(melt.df) == "variable"] <- "parameters"
      group.id="parameters"
      colour.id="parameters"
    }
    if (line.alpha==0.8)
      line.alpha=0.6
    if (point.alpha==1)
      point.alpha=0.6
    if (is.null(y.lab))
      y.lab = ""
  } else {
    # 1 parameter
    group.id=NULL
    colour.id=NULL
    if (is.null(y.lab))
      y.lab = trace.id
  }

  require(ComMA)
  gg.plot <- ComMA::ggLineWithPoints(melt.df, x.id="state", y.id="value", group.id=group.id,
                                colour.id=colour.id, point.size=point.size, point.alpha=point.alpha,
                                title=title, x.lab=x.lab, y.lab=y.lab,
                                line.or.point=line.or.point, line.alpha=line.alpha,
                                verbose=verbose, ...)
  # turns off clipping
  gt <- ComMA::unclip.ggplot(gg.plot)
  return(gt)
}

