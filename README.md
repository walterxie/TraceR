# Simplified Tracer in R

In development. 

##Not only the subset, but also more

This package is a simplified implementation of Tracer v1.6 \url{http://beast.bio.ed.ac.uk/tracer} in R.  

It can take MCMC log files from BEAST 1 \url{http://beast.bio.ed.ac.uk},
BEAST 2 \url{http://www.beast2.org}, or MrBayes \url{http://beast.bio.ed.ac.uk}.

This package is not only the subset of Tracer v1.6, 
but also aims to provide more visualizations to understand the Bayesian MCMC result, 
and also to make batch processing of multipe logs easier. 

##Citation
Rambaut A, Suchard MA, Xie D & Drummond AJ (2014) Tracer v1.6

##Installation

You can use the **devtools** *install\_github()* function to install the lastest development version directly from the GitHub.

```R
library("devtools")
devtools::install_github("walterxie/TraceR")
library("TraceR")
```

To see all exported functions:
```R
help(package = "TraceR")
```

##Tutorial



