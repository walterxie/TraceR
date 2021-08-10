# Tracer in R

In development. 

## Not only the subset, but also more

This package is initially developed as a simplified implementation 
of [Tracer v1.7](http://beast.community/tracer) in R.
It can take an advantage of ggplot2 to create high quality images, 
and also make batch processing of multiple logs easier.

Through the development, it is growing to a package not only used as 
a simple version of Tracer, 
but also aiming to provide more post-analysis methods to understand   
the posterior distributions from the Bayesian phylogenetic inference using MCMC. 
In addition, coverage tests for a model validation are recently added. 

1. For example, `readState` can extract the summary of operator proposals from 
   [BEAST 2](http://www.beast2.org) `*.state` log file.

```
readState("data/star.beast.state"")
```

2. Summarising tree statistics from the posterior trees logged by BEAST 2, 
   such as the total branch length, tree height, etc.

```
tre.sta.df <- readTrees("data/RSV2long.trees")
tre.sta <- analyseTreeStats(tre.sta.df)
```

3. The [pipeline](examples/Pipeline.R) to generate a coverage-test report 
for validating Bayesian phylogenetic models.



## Citation

Rambaut A, Drummond AJ, Xie D, Baele G and Suchard MA (2018),  
Posterior summarisation in Bayesian phylogenetics using Tracer 1.7, 
Systematic Biology. syy032. 
[doi:10.1093/sysbio/syy032](https://doi.org/10.1093/sysbio/syy032)


For the coverage-test pipeline, please cite:

Drummond AJ, Xie D, Mendes F  (),
LinguaPhylo: a probabilistic model specification language for reproducible phylogenetic analyses,
in preparation.

## Installation

You can use the **devtools** *install\_github()* function to install 
the latest development version directly from the GitHub.

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

Please note: [ComMA](https://github.com/walterxie/ComMA) package here   
is only used for plotting. You can use ggplot2 instead.
ComMA will be replaced by the original ggplot2 code, 
and then removed from the future dependecy.  

```R
devtools::install_github("walterxie/ComMA")
library("ComMA")
```


## Tutorial (TODO update is required)



