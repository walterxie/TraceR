# Author: Walter Xie
# Accessed on 27 Jul 2021

#' @name Coverage
#' @title Coverage Report from Model Validation
#'
#' @description
#' The pipeline to generate a coverage report of the model validation.
#'
#' 1. Summarise traces statistics for every simulations (normally starting 110,
#'    where 10 extra simulations are used for replacing any low ESS results
#'    in the next step), and save each result into a tsv file. It looks like:
#'    trace           posterior  likelihood  ...
#'    mean            -16848.58  -16269.03   ...
#'    stderr.of.mean	0.249      0.200       ...
#'    ...
#'    Note: if tree stats are logged into BEAST log file,
#'    you can set \code{tree.file=NA} to skip creating tree stats tsv file.
#'    See \code{\link{summariseTracesAndTrees}}.
#'
#' 2. Select 100 simulations where the ESS of every parameters are guaranteed
#'    >= 200. If not, then replace it to the one from the extra 10 sequentially,
#'    and check ESS. If all extra 10 are used but there exists any low-ESS simulations,
#'    then the pipeline stops and inform to re-run all simulations with longer
#'    MCMC chain length. See \code{\link{selectResultByESS}}.
#'
#' 3. Create a final summary file for visualisations.
#'    See \code{\link{summariseParameters}} and \code{\link{summariseTrueValues}}.
#'
#' 4. Plot coverage.
#'
#'
#' @details
#' \code{selectResultByESS} selects 100 simulation results, where the ESS of
#' every parameters are guaranteed >= 200.
#' The inputs files are produced by \code{\link{summariseTracesAndTrees}}.
#' Set \code{tree.file.postfix=NA} and \code{extr.tree.files=c()},
#' if it does not require tree stats tsv files.
#'
#' @param i.sta,i.end  The start index or end index of batch runs,
#'                     also used in the file names.
#' @param prefix  The prefix of statistics summary files.
#' @param file.postfix,tree.file.postfix
#'                The extension of statistics summary files.
#'                This function uses \code{\link{readTracesTSV}},
#'                The file extension can be different, but the file
#'                has to be tab-delimited text file.
#' @param file.steam.fun,extra.file.steam.fun
#'                The function to extract the file steam,
#'                which is used as a key to identify the simulation.
#' @param extra.files,extr.tree.files
#'                The extra files used to replace the low ESS results.
#' @keywords Coverage
#' @export
#' @examples
#' WD = file.path("~/WorkSpace/MyValidations/")
#' setwd(WD)
#' all.stats = list.files(pattern = "_([0-9]+).tsv")
#' all.stats # a char vector
#' extra.stats = all.stats[grep("_10([0-9]).tsv", all.stats, ignore.case = T)]
#' extra.stats # extra 10 simulations: *_100.tsv ~ *_109.tsv
#'
#' # skip tree stats
#' sele.res <- selectResultByESS(i.sta=0, i.end=99, prefix="sim",
#'                               tree.file.postfix=NA,
#'                               extra.files = extra.stats)
#' low.ess <- sele.res$lowESS
#' sele.list <- sele.res$selected
#' names(sele.list)
#'
#' @rdname Coverage
selectResultByESS <- function(i.sta=0, i.end=99, prefix="sim",
                              file.steam.fun=function(prefix, i){ paste0(prefix, "_", i) },
                              file.postfix="tsv", tree.file.postfix="trees.tsv",
                              extra.file.steam.fun=function(f){ sub('\\.tsv$', '', f) },
                              extra.files = c(), extr.tree.files=c() ) {
  require(tidyverse)
  # record 100 files whose ESS >= 200
  tracesDF <- list()
  # if any of 100 has <200 ESS, then it will be replaced by one of extra 10
  lowESS <- tibble(origin=character(),minESS=numeric(),replace=character())
  extr <- 1
  for(i in i.sta:i.end) {
    fn <- file.steam.fun(prefix, i)
    fi <- paste0(fn, ".", file.postfix)
    cat("Load ", fi, "...\n")

    res <- readTracesTSV(traces.file=fi)
    minESS <- res$minESS
    traces <- res$traces

    if (!is.na(tree.file.postfix)) {
      # if tree stats file exists
      tre.sta.fi <- paste0(fn, ".", tree.file.postfix)
      if (file.exists(tre.sta.fi)) {
        cat("Load ", tre.sta.fi, "...\n")
        # add tree stats (in cols)
        traces <- addTreeStats(tre.sta.fi, traces)
      }
    }

    # if low ESS
    if (minESS < 200) {
      # skip tree stats, if length(extr.tree.files)<=0
      selected <- pullAnthorFromExtra(extr,
                                      extra.file.steam.fun=extra.file.steam.fun,
                                      extra.files=extra.files,
                                      extr.tree.files=extr.tree.files )
      # selected extra
      file.selected <- selected$file.selected
      cat("Replace the low ESS result", fi, "to", file.selected, "\n")

      fn.sel <- extra.file.steam.fun(file.selected)
      tracesDF[[fn.sel]] <- selected$traces

      # sync current index in the pool
      extr <- selected$curr.idx

      # record low ESS
      lowESS <- lowESS %>% add_row(origin=fn,minESS=minESS,replace=fn.sel)
      if (nrow(selected$lowESS)>0) {
        lowESS <- lowESS %>% rbind(selected$lowESS)
      }

    } else {
      tracesDF[[fn]] <- traces
    }
  } # end for loop

  list(selected = tracesDF, lowESS = lowESS)
}

### private functions

# add tree stats DF (cols) into traces DF
addTreeStats <- function(tree.stats.file="sim_0.trees.tsv", traces) {
  # add tree stats
  tre.sta <- try(read_tsv(tre.sta.fi, col_types = cols())) %>% select(!trace)
  # they should have the same rows as stats file
  if (! all(traces$trace == tre.sta$trace) )
    stop("Trace stats ", nrow(traces), " should contain the same stats in 'trace' column ",
         "as the tree stats file ", tre.sta.fi, " !\n")
  return( traces %>% left_join(tre.sta, by = "trace") )
}

#sim_0.tsv
#trace	posterior	likelihood
#mean	-1395.71800741416	-1260.57129157319
#stderr.of.mean	0.948421280126143	0.142003271885284
pullAnthorFromExtra <- function(curr.idx,
                                extra.file.steam.fun=function(f){ sub('\\.tsv$', '', f) },
                                extra.files=c(), extr.tree.files=c()) {
  require(tidyverse)
  if (curr.idx >= length(extra.files)) # length(extra.files) extra
    stop("Not enough extra simulations to replace low-ESS simulations !\n",
         "Please re-run all simulations with longer MCMC chain length.\n",
         curr.idx, " >= ", length(extra.files))

  lowESS <- tibble(origin=character(),minESS=numeric(),replace=character())

  while (curr.idx < length(extra.files)) {
    # pick up "curr.idx" result from extra pool
    et.fi <- extra.files[curr.idx]

    # must incl trace column
    res <- readTracesTSV(traces.file=et.fi)
    minESS <- res$minESS
    traces <- res$traces

    # all ESS >= 200 then break, otherwise continue the loop
    if (minESS >= 200) {
      if (length(extr.tree.files)>0) {
        tre.sta.fi <- extr.tree.files[curr.idx]
        stopifnot(file.exists(tre.sta.fi))
        # add tree stats
        traces <- addTreeStats(tre.sta.fi, traces)
      }
      # point to the next before break loop
      curr.idx <- curr.idx + 1
      break
    } # end if minESS >= 200

    # record low ESS
    fn.extra <- extra.file.steam.fun(et.fi)
    lowESS <- lowESS %>% add_row(origin=fn.extra,minESS=minESS,replace="")

    # must increase after pick up tree stats from extra pool
    curr.idx <- curr.idx + 1
  } # end while

  list(file.selected=et.fi, curr.idx=curr.idx, traces=traces,
       minESS = minESS, lowESS = lowESS )
}
### end of private functions

#' @details
#' \code{summariseParameters} creates a file for each given parameters
#' storing the final stats summary of the posterior.
#' It also returns a list of data frame containing a parameter stats summary in each,
#' where rows are simulations and columns are statistics.
#'
#' @param selected   The list of data frames containing traces summary,
#'                   produced by \code{\link{selectResultByESS}}.
#' @param params     The vector of parameter names in the BEAST log files.
#' @param stats.name The vector of names of statistics. They have to be one of
#'                   names from the 1st column 'trace' in *.tsv file.
#' @param output.file.fun The function to create the summary file names
#'                   from the given parameters in the BEAST log.
#'                   Set to NA, if you do not want to create the files.
#' @keywords Coverage
#' @export
#' @examples
#' summ <- summariseParameters(sele.list,
#'            params = c("mu","Theta", "psi.treeLength", "psi.height"))
#' min(summ$minESS)
#'
#' @rdname Coverage
summariseParameters <- function(selected=list(),
          params = c("mu","Theta", "psi.treeLength", "psi.height"),
          stats.name = c("mean", "HPD95.lower", "HPD95.upper", "stdev", "ESS"),
          output.file.fun=function(x){ paste0(x, ".tsv") }) {
  require(tidyverse)
  cat("\nSelect ", length(selected)," simulations : ", paste(names(selected), collapse = ", "), ".\n")
  cat("Summarise BEAST parameters : ", paste(params, collapse = ", "), ".\n")

  minESS <- c()
  param.summaries <- list()
  for (pa in params) {
    cat("Analyse parameter : ", pa, "...\n")
    df <- tibble(trace=stats.name)

    pa.stats <- NULL
    for(i in 1:length(selected)) {
      fns <- names(selected)[i]
      # "mean", "HPD95.lower", "HPD95.upper", "ESS"
      pa.stats <- selected[[i]] %>% select(trace, !!pa) %>%
        filter(trace %in% stats.name)
      stopifnot(ncol(pa.stats) > 1 && nrow(pa.stats) > 0)

      df <- df %>% left_join(pa.stats, by = "trace") %>%
        rename(!!fns := !!pa)
    }

    ESS <- df %>% filter(trace == "ESS") %>% select(!trace) %>% unlist
    # to numeric
    tmp.minESS <- min(ESS %>% as.numeric)
    cat(pa, " min ESS = ", tmp.minESS, "\n")

    stopifnot("trace" %in% names(df))
    # keep 'trace' in previous order
    df <- df %>% mutate(trace = factor(trace, levels = unique(trace))) %>%
      # rotate df to make rows be simulations, cols are stats
      gather(key = simulation, value = value, -trace) %>%
      group_by(trace) %>% # keep 'trace' in previous order
      spread(key = trace, value = value)
    # save to list
    param.summaries[[pa]] <- df

    if (!is.na(tmp.minESS) && tmp.minESS >= 200) {
      if (is.function(output.file.fun))
        write_tsv(df, paste0(pa, ".tsv"))
    } else
      warning("Summary not generated ! ", pa, " min ESS = ", tmp.minESS, "\n")
    minESS <- c(minESS, tmp.minESS)
  }

  list(param.summaries=param.summaries, minESS=minESS)
}

#' @details
#' \code{summariseTrueValues} return a data frame to contain all of
#' the true values from all LPhy simulations.
#' LPhy does not log tree stats as a parameter,
#' but they can be computed from the true tree.
#'
#' The returned data frame for "true values" uses the same format of BEAST 2
#' model validation pipeline \url{https://github.com/rbouckaert/DeveloperManual}.
#'
#' @param selected.fn.steam The vector of file name steams selected
#'                          by \code{selectResultByESS} to summarise
#'                          the coverage.
#' @param log.file.fun  The function to get one-line log file name,
#'                      containing true values from LPhy simulations,
#'                      where one LPhy log file per simulation.
#' @param tree.file.fun The function to get one-tree log file name,
#'                      containing true trees from LPhy simulations,
#'                      where one LPhy tree file per simulation.
#'                      The tree stats are fixed to "total.br.len" (total branch length)
#'                      and "tree.height" at the moment.
#'                      Set to NA, if tree stats (e.g total branch length,
#'                      root height) is not required.
#' @param params  The vector of parameter names in LPhy.
#' @keywords Coverage
#' @export
#' @examples
#' # list.files(pattern = "_true.log")
#' df.tru <- summariseTrueValues(names(sele.list), params=c("μ","Θ"))
#' getwd()
#' write_tsv(df.tru, "trueValue.tsv")
#'
#' @rdname Coverage
summariseTrueValues <- function(selected.fn.steam=c(),
                                log.file.fun=function(x){ paste0(x,"_true.log") },
                                tree.file.fun=function(x){ paste0(x,"_true_ψ.trees") },
                                params=c("μ","Θ") ) {
  require(tidyverse)
  df.tru <- tibble(parameter = params)

  if (is.function(tree.file.fun)) {
    require(ape)
    require(phytools)
    # add 2 tree stats
    df.tru <- tibble(parameter = c(params, "total.br.len","tree.height"))
  }
  cat("\nSelect ", length(selected.fn.steam)," true-value files : ",
      paste(selected.fn.steam, collapse = ", "), ".\n")
  cat("Summarise LPhy parameters : ", paste(params, collapse = ", "), ".\n")

  # true values from a file
  tru <- NULL
  # selected is a list of file steams
  for(lg in selected.fn.steam) {
    lg.fi <- log.file.fun(lg)
    cat("Load ", lg.fi, "...\n")
    stopifnot(file.exists(lg.fi))

    # must 1 line
    tru <- read_tsv(lg.fi, col_types = cols()) %>% select(params) %>% unlist # need vector here

    if (is.function(tree.file.fun)) {
      # add tree stats
      tre.fi <- tree.file.fun(lg)
      stopifnot(file.exists(tre.fi))

      tru.tre <- read.nexus(tre.fi)
      cat("Load ", tre.fi, "having", Ntip(tru.tre), "tips ...\n")

      # total branch len and tree height
      tru <- c(tru, sum(tru.tre$edge.length), max(nodeHeights(tru.tre)))
    }

    df.tru <- try(df.tru %>% add_column(!!(lg) := tru))
  }
  stopifnot("parameter" %in% names(df.tru))
  # keep 'parameter' in previous order
  df.tru <- df.tru %>% mutate(parameter = factor(parameter, levels = unique(parameter)))
  return( df.tru %>% # start rotating
            gather(key = simulation, value = value, -parameter) %>%
            group_by(parameter) %>% # keep 'parameter' in previous order
            spread(key = parameter, value = value) )
}

#' @details
#' \code{reportCoverage} produces a data frame as the report for
#' the selected simulations to measure how many true values are falling into
#' or outside the 95% HPD interval in the sampled posterior.
#'
#' Note: the same parameter may be given different names between LPhy script
#' and BEAST XML/log, please ensure that you match them correctly.
#'
#' @param df.posterior The data frame containing the final stats summary of
#'                     the posterior of the selected parameter.
#' @param df.tru.val   The data frame containing all of the true values
#'                     from all LPhy simulations. Created by \code{summariseTrueValues}
#' @param tru.val.par  The parameter name in LPhy script, which may be different
#'                     to one in the BEAST log.
#' @param sim.coln     The default column name as the key both in
#'                     \code{df.posterior} and \code{df.tru.val}.
#' @keywords Coverage
#' @export
#' @examples
#' # df.pos <- read_tsv("mu.tsv")
#' # df.tru <- read_tsv("trueValue.tsv")
#' covg <- reportCoverage(df.pos, df.tru, tru.val.par="μ")
#' covg
#' write_tsv(covg, "mu-coverage.tsv")
#'
#' @rdname Coverage
reportCoverage <- function(df.posterior, df.tru.val, tru.val.par="μ",
                           sim.coln="simulation") {
  require(tidyverse)
  stopifnot(tru.val.par %in% colnames(df.tru.val))
  stopifnot(sim.coln %in% colnames(df.tru.val))
  stopifnot(sim.coln %in% colnames(df.posterior))
  cat("Compute coverage for ", tru.val.par, "from", nrow(df.tru.val) ,
      "simulations ...\n")

  # tru.val.par="μ"
  tru.vals <- df.tru.val %>% select(!!sim.coln, !!tru.val.par)

  # merge
  covg <- tru.vals %>% inner_join(df.posterior, by = sim.coln)  %>%
    mutate(is.in = (!!tru.val.par >= HPD95.lower & !!tru.val.par <= HPD95.upper) )

  print(covg, n = 5)
  return(covg)
}







