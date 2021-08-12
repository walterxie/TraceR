
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

# params to report
beast.params = c("mu","Theta", "r_0", "r_1", "r_2",
                 "kappa.1", "kappa.2", "kappa.3",
                 "pi_0.A", "pi_0.C", "pi_0.G", "pi_0.T",
                 "pi_1.A", "pi_1.C", "pi_1.G", "pi_1.T",
                 "pi_2.A", "pi_2.C", "pi_2.G", "pi_2.T",
                 "psi.treeLength", "psi.height")
summ <- pipCreateParameterSummaries(sele.list, params = beast.params)
print(summ$param.summaries[["mu"]], n=3)

lphy.params = c("μ","Θ","r_0","r_1","r_2","κ_0","κ_1","κ_2",
                "π_0_0","π_0_1","π_0_2","π_0_3",
                "π_1_0","π_1_1","π_1_2","π_1_3",
                "π_2_0","π_2_1","π_2_2","π_2_3")
# the order of parameters same to BEAST parameters
df.tru <- pipCreateTrueValueSummaries(names(sele.list), params=lphy.params,
                                                add.tree.stats=TRUE)
print(df.tru, n=3)


covg <- reportCoverages(beast.params = beast.params,
                        lphy.params = c(lphy.params, "total.br.len","tree.height"))

covg

inOut <- read_tsv("mu-coverage.tsv", col_types = cols())
p <- ggCoverage(inOut, x.lab="True mu value")
p
# ggsave("mu.pdf", p, width = 4, height = 3)

# log scale
inOut <- read_tsv("Theta-coverage.tsv", col_types = cols())
p <- ggCoverage(inOut, x.lab=paste0("True log-theta value"),
                y.lab="Log-mean posterior", x.txt=1, y.txt=9000)
# log scale and fix labels and text
p <- p + scale_x_log10(limits = c(1,1e4),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(1,1e4), breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
p


