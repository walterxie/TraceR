# plot the coverage graph for well-calibrated study
# https://github.com/walterxie/TraceR
# install libraires below
library("TraceR")
library(ape)
require(phytools)
require(tidyverse)

# change to your path
setwd("~/WorkSpace/toroidalDiffusion/tmp")

#Sample	posterior	likelihood	prior	drift.1	drift.2	driftCorr	mu.1	mu.2	sigma.1	sigma.2
#lambda	psi2.height	psi2.treeLength	PhyloWrappedBivariateDiffusion	alpha3
# Sample	drift_0	drift_1	driftCorr	mu_0	mu_1	sigma_0	sigma_1	lambda
tree.prefix = "psi2"
# paramter names logged in beast log
params <- c("drift.1", "drift.2", "driftCorr", "mu.1", "mu.2", "sigma.1", "sigma.2",
            "lambda", paste0(tree.prefix, ".height"), paste0(tree.prefix, ".treeLength"))

res <- tibble()
# 100 beast results
for (i in 0:99) {
  cat("Start analysis ", i, " ...\n")

  fil <- file.path("log_files", paste0("test-", i,".log"))
  stopifnot(file.exists(fil))
  # read MCMC log
  mcmc.log <- readMCMCLog(fil)
  traces <- TraceR::getTraces(mcmc.log, burn.in=0.1)
  stats <- TraceR::analyseTraces(traces) %>% column_to_rownames(var = "trace")

  # print column if ESS < 200
  all.ess <- stats %>% filter(row.names(stats) == "ESS") %>%
    # make sure they are numbers
    mutate_if(is.character, as.numeric)
  cat("min ESS = ", min(all.ess), " at ", colnames(all.ess)[which(all.ess==min(all.ess))], "\n")

  # extract HPD intervals
  hpd.lo <- stats %>% select(!!params) %>% filter(row.names(stats) == "HPD95.lower")
  hpd.up <- stats %>% select(!!params) %>% filter(row.names(stats) == "HPD95.upper")

  # beast posterior trees
  tres <- file.path("trees", paste0("test-", i,".trees"))
  tres.b <- readTrees(tres, burn.in = 0.1)
  tre.sta <- analyseTreeStats(tres.b)
  tre.hpd <- tre.sta %>% filter(row.names(stats) == "HPD95.lower" |
                                  row.names(stats) == "HPD95.upper")

  # compute tree stats HPD
  hpd.lo[paste0(tree.prefix, ".height")] <- tre.hpd[tre.hpd$trace == "HPD95.lower", "tree.height"]
  hpd.lo[paste0(tree.prefix, ".treeLength")] <- tre.hpd[tre.hpd$trace == "HPD95.lower", "total.br.len"]
  hpd.up[paste0(tree.prefix, ".height")] <- tre.hpd[tre.hpd$trace == "HPD95.upper", "tree.height"]
  hpd.up[paste0(tree.prefix, ".treeLength")] <- tre.hpd[tre.hpd$trace == "HPD95.upper", "total.br.len"]

  # combine all posterior stats here
  df <- tibble(Parameter = names(hpd.lo), HPD95.lower = hpd.lo %>% unlist)
  df$HPD95.upper <- hpd.up %>% unlist

  # true value
  tru.fil <- file.path("true_logs", paste0("test-", i,".log"))
  tru.v <- read_delim(tru.fil, delim="\t", comment = "#")
  tru.v <- tru.v %>% select(!contains("siteRates") & !contains("Q_") & !Sample) %>% unlist
  # match names
  #feqs <- names(tru.v)[1:4]
  # Replace the numbers with corresponding letters
  #index <- as.numeric(gsub("frequencies_", "", feqs)) + 1
  #names(tru.v)[1:4] <- paste0("frequencies.", bases[index])

  # true tree
  tru.tre.f <- file.path("true_trees", paste0("test-", i,"_", tree.prefix, ".trees"))
  tru.tre <- read.nexus(tru.tre.f)
  # true tree stats
  tru.v[paste0(tree.prefix, ".height")] <- max(nodeHeights(tru.tre))
  tru.v[paste0(tree.prefix, ".treeLength")] <- sum(tru.tre$edge.length)

  # combine all truth here
  tru.df <- tibble(Parameter = names(tru.v), True.Value = tru.v)

  # change names to the same
  tru.df$Parameter <- sub("_", ".", tru.df$Parameter)
  # have to sysnc lphy array index to beast
  tru.df$Parameter <- sub("\\.1", "\\.2", tru.df$Parameter)
  tru.df$Parameter <- sub("\\.0", "\\.1", tru.df$Parameter)


  # merge posterior dataframe with true
  df <- merge(df, tru.df)
  stopifnot(nrow(df) == length(params))

  df_sorted <- df %>%
    mutate(Parameter = factor(Parameter, levels = params)) %>%
    arrange(Parameter) %>%
    mutate(i = i)

  # add ith processed result into the final dataframe
  res <- bind_rows(res, df_sorted)

}
# Convert HPD95.lower and HPD95.upper from character to numeric
res <- res %>%  mutate(
    i = i + 1,
    HPD95.lower = as.numeric(HPD95.lower),
    HPD95.upper = as.numeric(HPD95.upper),
    is.in = True.Value >= HPD95.lower & True.Value <= HPD95.upper
)
res

# coverage
cov <- res %>%
  group_by(Parameter) %>%        # Group by Parameter
  summarize(Coverage = sum(is.in))
cov

facet_labels <- cov %>%
  mutate(Label = paste0(Parameter, " (", Coverage, "%)")) %>%
  select(Parameter, Label) %>%
  deframe()

library(ggplot2)

p <- ggplot(res, aes(x = True.Value, colour = is.in) ) +
  geom_errorbar(aes(ymin=HPD95.lower, ymax=HPD95.upper), width = 0, alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +  # Diagonal dashed line
  #scale_x_log10() +
  #scale_y_log10() +
  facet_wrap(~Parameter, labeller = labeller(Parameter = facet_labels), scales = "free") +
  labs(x = "True value", y = "Posterior") +
  theme_minimal() + theme(legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
p

ggsave("../DihedralAngle-well-calbr.pdf", p, width = 8, height = 6)


