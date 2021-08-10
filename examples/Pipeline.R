
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
# parameters in BEAST log
names(res$stats)

all.stats = list.files(pattern = "_([0-9]+).tsv")
extra.stats = all.stats[grep("_10([0-9]).tsv", all.stats, ignore.case = T)]
extra.stats # extra 10 simulations: *_100.tsv ~ *_109.tsv

sele.res <- selectResultByESS(i.sta=0, i.end=99, prefix="al2",
                               tree.file.postfix=NA, extra.files = extra.stats)

sele.res$lowESS
sele.list <- sele.res$selected
sele.list[[1]]
names(sele.list)

summ <- summariseParameters(sele.list,
                            params = c("mu","Theta", "r_0", "r_1", "r_2",
                                       "psi.treeLength", "psi.height"))

cat("min ESS = ", paste(summ$minESS, collapse = ", "), "\n")
cat("min of min ESS = ", min(summ$minESS), "\n")

true.log.files = list.files(pattern = "_true.log")
true.log.files
# LPhy parameters in the true-value log
read_tsv("al2_0_true.log") %>% names
# list.files(pattern = "_true_ψ.trees")

# the order of parameters same to BEAST parameters
df.tru <- summariseTrueValues(names(sele.list),
                              params=c("μ","Θ", "r_0", "r_1", "r_2"),
                              add.tree.stats=TRUE)
df.tru
getwd()
write_tsv(df.tru, "trueValue.tsv")

# df.pos <- summ$param.summaries[["mu"]]
df.pos <- read_tsv("mu.tsv")
inOut <- markInOut(df.pos, df.tru, tru.val.par="μ")
write_tsv(inOut, "mu-coverage.tsv")

nrow(inOut[inOut$is.in==T,])/nrow(inOut)

p <- ggCoverage(inOut, x.lab="True mu value")
p
# ggsave("mu.pdf", p, width = 4, height = 3)



