# analyze nanodrop data from pilot dna extractions
# kms sept 2016

############################################################
# libraries
############################################################

# install/load packages

library("tidyverse")
library("cowplot")

############################################################
# load data
############################################################

# load data
nd_dat <- read.table("data/nanodrop/pilot_nanodrop_master.txt", h = T, row.names = NULL)

# remove repeated measurements
nd_dat <- nd_dat %>%
  filter(measurement == 1)

# labels for facets
stat_names <- c(
  `r_260.280` = "260/280",
  `r_260.230` = "260/230",
  `yield` = "Total Yield (ng)",
  `1` = "Attempt #1",
  `2` = "Attempt #2",
  `3` = "Attempt #3"
)

# plot results of DNA extractions
nd_dat %>%
  filter(r_260.230 > 0, r_260.230 < 5) %>%
  gather(key = stat, value = value, -sex, -id, -protocol) %>% 
  mutate(sex = toupper(sex)) %>%
  filter(stat %in% c("yield", "r_260.280","r_260.230")) %>%
  ggplot(aes(y = value , x = sex, fill = sex, group = interaction(stat,sex)))+
  #stat_summary(fun.data = "mean_cl_boot", size = 2)+
  geom_dotplot(binaxis = "y", stackdir = "centerwhole", dotsize = 2)+
  facet_grid(stat ~ protocol, scales = "free", switch = "y", labeller = as_labeller(stat_names)) + 
  xlab("Sex")+
  ylab(NULL)+
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        strip.text = element_text(face = "bold"))





