
library("devtools")
remove.packages("TraceR")
devtools::install_github("walterxie/TraceR")
library("TraceR")

# cp ../alpha2/*.log ../alpha2/*.trees .
WD = file.path("~/WorkSpace/linguaPhylo/manuscript/test/")
setwd(WD)

log.files = list.files(pattern = "_([0-9]+).log")
log.files # 110 simulations
# parameters in BEAST log
all.beast.params <- pipCreateSimulationSummaries(log.files, burn.in=0.1)
all.beast.params

all.stats = list.files(pattern = "_([0-9]+).tsv")
all.stats
# select vaild results
sele.list <- pipSelectValidResults(i.sta=0, i.end=99, prefix="al2",
                                   extra.tree.file.fun=NA)



beast.params = c("mu","Theta", "r_0", "r_1", "r_2", "psi.treeLength", "psi.height")
summ <- summariseParameters(sele.list, params = beast.params)
for (pa in beast.params) {
  df <- summ$param.summaries[[pa]]
  write_tsv(df, paste0(pa, ".tsv"))
  cat("Write ", paste0(pa, ".tsv"), "\n")
}
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



