# Simplified Tracer in R

In development. 

## Not only the subset, but also more

This package is a simplified implementation of [Tracer v1.7](http://beast.community/tracer) in R.  

This package is not only the subset of Tracer, 
but also aims to provide more visualizations to understand the Bayesian MCMC result, 
and also to make batch processing of multipe logs easier. 

For example, `readState` can extract the summary of operator proposals from 
[BEAST 2](http://www.beast2.org) `*.state` log file.

```
readState("data/star.beast.state"")
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
devtools::install_github("walterxie/TraceR")
library("TraceR")
```

```
# read MCMC log
mcmc.log <- readMCMCLog("star.beast.log")
# preprocess to remove burn in
mcmc.traces <- getTraces(mcmc.log, burn.in=0.1)
# get stats
stats <- analyseTraces(traces)
```


Or install from `TraceR_*.tar.gz`.

```R
install.packages("TraceR_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library("TraceR")
```

To see all exported functions:
```R
help(package = "TraceR")
```

[ComMA](https://github.com/walterxie/ComMA) package here only provides plotting functions.  
It will be replaced by orignal ggplot2 code, and then removed in future development.  

```R
devtools::install_github("walterxie/ComMA")
library("ComMA")
```


## Tutorial (TODO update is required)



