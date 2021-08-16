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
#' @keywords Coverage
#' @export
#' @examples
#' param = "mu"
#' df.pos <- read_tsv(paste0(param, ".tsv"))
#' df.tru <- read_tsv("trueValue.tsv")
#' inOut <- markInOut(df.pos, df.tru, tru.val.par="μ")
#' write_tsv(inOut, paste0(param, "-coverage.tsv"))
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

#' @details
#' The data saturation in the simluated data may cause MCMC defficult to converge.
#' \code{ggDataSaturation} returns a \code{\link{ggplot}} object to visualise
#' the maximum distances of overall partitions \code{max(root height * mu * r)}.
#' It is also coloured by the coverages.
#'
#' @param gene.dist The data frame contains simulation and distances.
#' @param x.col,y.col,colour.col The column names.
#' @param x.lab,y.lab   The label of x-axis or y-axis".
#' @keywords Coverage
#' @export
#' @examples
#'
#' ### create data frame for data saturation visualisation
#'
#' df.mu <- read_tsv("mu-coverage.tsv", col_types = cols()) %>%
#'   select("simulation", "true.val") %>% rename(mu=true.val)
#' df.r0 <- read_tsv("r_0-coverage.tsv", col_types = cols()) %>%
#'   select("simulation", "true.val") %>% rename(r0=true.val)
#' df.r1 <- read_tsv("r_1-coverage.tsv", col_types = cols()) %>%
#'   select("simulation", "true.val") %>% rename(r1=true.val)
#' df.r2 <- read_tsv("r_2-coverage.tsv", col_types = cols()) %>%
#'   select("simulation", "true.val") %>% rename(r2=true.val)
#' df.hei <- read_tsv("psi.height-coverage.tsv", col_types = cols()) %>%
#'   select("simulation", "true.val") %>% rename(height=true.val)
#'
#' df.colour <- read_tsv("Theta-coverage.tsv", col_types = cols()) %>%
#'   select("simulation", "is.in")
#'
#' gene.dist <- df.mu %>% inner_join(df.hei, by = "simulation") %>%
#'   inner_join(df.r0, by = "simulation") %>%
#'   inner_join(df.r1, by = "simulation") %>%
#'   inner_join(df.r2, by = "simulation") %>%
#'   mutate(max.r = pmax(r0, r1, r2)) %>%
#'   mutate(distance = height*mu*max.r) %>%
#'   select(-r0, -r1, -r2) %>%
#'   inner_join(df.colour, by = "simulation")
#'
#' gene.dist
#' #   simulation  mu     height max.r  distance is.in
#' #   <chr>       <dbl>  <dbl>  <dbl>   <dbl>   <lgl>
#' # 1 al2_0      0.0165   72.0  1.36    1.62    TRUE
#' # 2 al2_1      0.00145  158.  1.99    0.458   TRUE
#'
#' p <- ggDataSaturation(gene.dist)
#' @rdname CovgGraph
ggDataSaturation <- function(gene.dist, x.col = "simulation", y.col = "distance",
                             colour.col = "is.in", x.lab="", y.lab="max(root height * mu * r)") {

  require("tidyverse")
  require("ggplot2")

  stopifnot( any(c(x,y,colour) %in% colnames(gene.dist)) )

  p <- ggplot(gene.dist, aes(x = !!sym(x.col), y = !!sym(y.col),
                             group = !!sym(colour.col),
                             colour = !!sym(colour.col)) ) +
    geom_point(aes(shape = !!sym(colour.col)), size = 0.2, alpha = 0.9) +
    scale_shape_manual(values=c(4, 1, 0))+
    geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.2) +
    geom_hline(yintercept = 0.5, linetype = "dotted", size = 0.2) +
    scale_y_log10() +
    xlab(x.lab) + ylab(y.lab) + guides(colour = FALSE, shape = FALSE) +
    theme_classic() +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  return(p)
}

