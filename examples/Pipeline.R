
library("devtools")
remove.packages("TraceR")
devtools::install_github("walterxie/TraceR")
library("TraceR")

WD = file.path("~/WorkSpace/linguaPhylo/manuscript/test/")
setwd(WD)

log.files = list.files(pattern = "_([0-9]+).log")
log.files
for(lg in log.files) {
  res <- summariseTracesAndTrees(lg, tree.file=NA)
}

all.stats = list.files(pattern = "_([0-9]+).tsv")
extra.stats = all.stats[grep("_10([0-9]).tsv", all.stats, ignore.case = T)]
extra.stats # extra 10 simulations: *_100.tsv ~ *_109.tsv

sele.res <- selectResultByESS(i.sta=0, i.end=99, prefix="al1",
                               tree.file.postfix=NA, extra.files = extra.stats)

low.ess <- sele.res$lowESS
sele.list <- sele.res$selected
sele.list[[1]]
names(sele.list)

summ <- summariseParameters(sele.list, params = c("mu","Theta", "r_0", "r_1",
                                                  "psi.treeLength", "psi.height"))

cat("min ESS = ", paste(summ$minESS, collapse = ", "), "\n")
cat("min of min ESS = ", min(summ$minESS), "\n")

