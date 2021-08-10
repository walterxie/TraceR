# Author: Walter Xie and Fábio K Mendes
# Accessed on 27 Jul 2021

#' @name CovgGraph
#' @title Coverage Test Graph in R
#'
#' @description
#' The visualization functions to summarise the model validation result.
#' The data frame provided to the plot function must contain
#' the following columns, and same column names:
#'
#' analysis  mean   HPD95.lower HPD95.upper   ESS   true.val  is.in
#' RSV2_58   0.438  0.327       0.554         948.  0.474     TRUE
#'
#' @details
#' \code{ggCoverage} returns a \code{\link{ggplot}} object to plot
#' the coverage of a model validation,
#' where the data frame has columns as statistics and rows as analyses.
#' The data frame can be a subset of full analyses, therefore the plot
#' will zoom into a better view to skip some outliers.
#'
#' @param df The data frame contains summary statistics of
#'           the model validation result.
#' @param cov.per The percentage of coverage, which indicates how frequent
#'                the "true" values are falling into the 95% HPD interval
#'                of the posterior.
#' @param transp  The transparency level of 95% HPD bar, default to 0.3.
#' @param x.lab,y.lab   The label of x-axis or y-axis".
#' @param x.max.lim,y.max.lim The maximum value in the axis,
#'                which can be used to adjust the x and y to the same scale.
#'                Default to NA.
#' @param x.txt,y.txt The position of "covg. =" in the figure.
#'                    Default to NA, which will be automatically assigned
#'                    by minimum and maximum values. They are only used
#'                    for log scale.
#' @param x.txt.just Where the text "covg. =" starts, default to 0.
#' @keywords Tracer
#' @export
#' @examples
#' param = "mu"
#' df.pos <- read_tsv(paste0(param, ".tsv"))
#' df.tru <- read_tsv("trueValue.tsv")
#' inOut <- markInOut(df.pos, df.tru, tru.val.par="μ")
#' # coverage
#' cov <- nrow(subset(inOut, is.in==TRUE)) / nrow(inOut)
#' p <- ggCoverage(inOut, x.lab=paste("True",param,"value"))
#' ggsave(paste0(param, ".pdf"), p, width = 4, height = 3)
#'
#' # zoom in
#' df.sub <- inOut %>% filter(mean < 0.045)
#' nrow(df.sub)
#' cov.per <- round(cov * 100)
#' p <- ggCoverage(df.sub, cov.per, x.lab=paste("True",param,"value"))
#' ggsave(paste0(param, "-sub-",nrow(df.sub),".pdf"), p, width = 4, height = 3)
#'
#' @rdname CovgGraph
ggCoverage <- function(df, cov.per=-1, transp = 0.3, x.lab="", y.lab="Mean posterior",
                       x.txt=NA, y.txt=NA, x.max.lim=NA, y.max.lim=NA, x.txt.just=0) {
  require("ggplot2")

  # colnames must have: "analysis mean HPD95.lower HPD95.upper ESS true.val is.in"

  # default to assume 100 runs
  if (cov.per < 0)
    cov.per = round(nrow(subset(df, is.in==TRUE)) / nrow(df) * 100)

  if (is.na(x.txt)) x.txt = min(df$true.val)
  if (is.na(y.txt)) y.txt = max(df$HPD95.upper)
  if (is.na(x.txt.just)) x.txt.just = max(df$HPD95.upper) * 0.1

  p <- ggplot(data=df, aes(x=true.val, y=mean, group = is.in, colour = is.in)) +
    geom_linerange(aes(ymin=HPD95.lower, ymax=HPD95.upper), size=1.2, alpha=transp) +
    geom_point(size=.2) +
    geom_abline(intercept = 0, slope = 1, color="black", linetype="dotted", size=.2) +
    annotate("text", x=x.txt, y=y.txt, label= paste("covg. =", cov.per, "%"),
             hjust = x.txt.just, size = 5) +
    xlab(x.lab) + ylab(y.lab) +
    guides(colour=FALSE) + theme_classic() + theme(text = element_text(size=15))

  # same scale in x and y
  if (!is.na(x.max.lim))
    p <- p + xlim(NA, x.max.lim)
  else
    p <- p + xlim(NA, y.txt)
  if (!is.na(y.max.lim))
    p <- p + ylim(NA, y.max.lim)
  else
    p <- p + ylim(NA, y.txt)

  return(p)
}


