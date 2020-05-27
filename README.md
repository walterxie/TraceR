# Simplified Tracer in R

In development. 

## Not only the subset, but also more

This package is a simplified implementation of [Tracer v1.6](http://beast.bio.ed.ac.uk/tracer) in R.  

It can take MCMC log files from [BEAST 1](http://beast.bio.ed.ac.uk),
[BEAST 2](http://www.beast2.org), or [MrBayes](http://beast.bio.ed.ac.uk).

This package is not only the subset of Tracer v1.6, 
but also aims to provide more visualizations to understand the Bayesian MCMC result, 
and also to make batch processing of multipe logs easier. 

## Citation
Rambaut A, Suchard MA, Xie D & Drummond AJ (2014) Tracer v1.6

## Installation

You can use the **devtools** *install\_github()* function to install the lastest development version directly from the GitHub.

```R
library("devtools")
devtools::install_github("walterxie/TraceR")
library("TraceR")
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

The dependency of [ComMA](https://github.com/walterxie/ComMA) package will be removed in future development.  

```R
devtools::install_github("walterxie/ComMA")
library("ComMA")
```


## Tutorial



