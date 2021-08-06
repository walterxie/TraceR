# Simplified Tracer in R

In development. 

## Not only the subset, but also more

This package is a simplified implementation of [Tracer v1.7](http://beast.community/tracer) in R.  

This package is not only the subset of Tracer, 
but also aims to provide more visualizations to understand the Bayesian MCMC result, 
and also to make batch processing of multiple logs easier. 

For example, `readState` can extract the summary of operator proposals from 
[BEAST 2](http://www.beast2.org) `*.state` log file.

```
readState("data/star.beast.state"")
```

Summarising tree statistics from the posterior trees logged by BEAST 2, 
such as the total branch length, tree height, etc.

```
tre.sta.df <- readTrees("data/RSV2long.trees")
tre.sta <- analyseTreeStats(tre.sta.df)
```

## Citation

Rambaut A, Drummond AJ, Xie D, Baele G and Suchard MA (2018),  
Posterior summarisation in Bayesian phylogenetics using Tracer 1.7, 
Systematic Biology. syy032. 
[doi:10.1093/sysbio/syy032](https://doi.org/10.1093/sysbio/syy032)

## Installation

You can use the **devtools** *install\_github()* function to install the lastest development version directly from the GitHub.

```R
library("devtools")
remove.packages("TraceR")
devtools::install_github("walterxie/TraceR")
library("TraceR")
```

Or install from `TraceR_*.tar.gz`.

```R
setwd("~/WorkSpace/TraceR")
install.packages("TraceR_?.?.?.tar.gz", repos = NULL, type = "source")
library("TraceR")
```

To see all exported functions:
```R
help(package = "TraceR")
```

Try the following commands:
```
# read MCMC log
mcmc.log <- readMCMCLog("data/star.beast.log")
# get traces and remove burn in
traces <- getTraces(mcmc.log, burn.in=0.1)
# get stats
stats <- analyseTraces(traces)
```

[ComMA](https://github.com/walterxie/ComMA) package here only provides plotting functions.  
It will be replaced by original ggplot2 code, and then removed in future.  

```R
devtools::install_github("walterxie/ComMA")
library("ComMA")
```


## Tutorial (TODO update is required)



