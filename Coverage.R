# https://github.com/walterxie/TraceR

library("TraceR")
library(ape)
require(phytools)
require(tidyverse)

setwd("~/WorkSpace/linguaPhylo.github.io/tutorials/ORC")

params <- c("Theta", "uclnMean", "sigma", "psi.height", "psi.treeLength")

res <- tibble()
for (i in 0:99) {
  cat("Start analysis ", i, " ...\n")

  fil <- file.path(paste0("ORC-", i,".log"))
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

  hpd.lo <- stats %>% select(!!params) %>% filter(row.names(stats) == "HPD95.lower")
  hpd.up <- stats %>% select(!!params) %>% filter(row.names(stats) == "HPD95.upper")

  tres <- file.path(paste0("ORC-", i,"_with_rates.trees"))
  tres.b <- readTrees(tres, burn.in = 0.1)
  tre.sta <- analyseTreeStats(tres.b)
  tre.hpd <- tre.sta %>% filter(row.names(stats) == "HPD95.lower" |
                                  row.names(stats) == "HPD95.upper")

  hpd.lo["phi.height"] <- tre.hpd[tre.hpd$trace == "HPD95.lower", "tree.height"]
  hpd.lo["phi.treeLength"] <- tre.hpd[tre.hpd$trace == "HPD95.lower", "total.br.len"]
  hpd.up["phi.height"] <- tre.hpd[tre.hpd$trace == "HPD95.upper", "tree.height"]
  hpd.up["phi.treeLength"] <- tre.hpd[tre.hpd$trace == "HPD95.upper", "total.br.len"]

  df <- tibble(Parameter = names(hpd.lo), HPD95.lower = hpd.lo %>% unlist)
  df$HPD95.upper <- hpd.up %>% unlist

  # true value
  tru.fil <- file.path("true", paste0("ORC-", i,".log"))
  tru.v <- read_delim(tru.fil, delim="\t", comment = "#")
  tru.v <- tru.v %>% select(!contains("siteRates") & !contains("Q_") & !Sample) %>% unlist
  # match names
  #feqs <- names(tru.v)[1:4]
  # Replace the numbers with corresponding letters
  #index <- as.numeric(gsub("frequencies_", "", feqs)) + 1
  #names(tru.v)[1:4] <- paste0("frequencies.", bases[index])

  tru.tre.f <- file.path("true", paste0("ORC-", i,"_psi.trees"))
  tru.tre <- read.nexus(tru.tre.f)

  tru.v["psi.height"] <- max(nodeHeights(tru.tre))
  tru.v["psi.treeLength"] <- sum(tru.tre$edge.length)

  tru.df <- tibble(Parameter = names(tru.v), True.Value = tru.v)
  df <- merge(df, tru.df)
  stopifnot(nrow(df) == length(params))

  df_sorted <- df %>%
    mutate(Parameter = factor(Parameter, levels = params)) %>%
    arrange(Parameter) %>%
    mutate(i = i)

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

ggsave("../UCLN-well-calbr.pdf", p, width = 8, height = 6)


