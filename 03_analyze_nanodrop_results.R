# analyze nanodrop data from pilot dna extractions
# kms sept 2016

############################################################
# libraries
############################################################

# install/load packages

library("tidyverse")
library("cowplot")
library("ggthemes")

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
  `1` = "Base PCL #1",
  `2` = "Base PCL #2",
  `3` = "Base PCL #3",
  `4` = "Phase Lock #3"
)

# plot results of DNA extractions
nd_dat %>%
  filter(r_260.230 > 0, r_260.230 < 5) %>%
  gather(key = stat, value = value, -sex, -id, -protocol) %>% 
  mutate(sex = toupper(sex)) %>%
  filter(stat %in% c("yield", "r_260.280","r_260.230")) %>%
  mutate(phase_lock = protocol == 4) %>%
  ggplot(aes(y = value, x = phase_lock))+
  geom_jitter(color = "grey75", width = 0.5)+
  stat_summary(fun.data = "mean_cl_boot", size = 1, aes(color = as.factor(phase_lock)))+
  facet_wrap(~stat, scales = "free_y", labeller = as_labeller(stat_names)) + 
  ylab(NULL)+
  theme_bw(base_size = 18)+
  theme(strip.background = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        strip.text = element_text(face = "bold"))+
  scale_color_fivethirtyeight(guide = guide_legend(title = "Phase lock?"))+
  xlab("Phase lock?")

# qubit results

nd_dat %>%
  ggplot(aes(x = yield, y = qubit_yield)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Nanodrop Yield (ng)")+
  ylab("Qubit yield (ng)")

with(nd_dat, lm(qubit_yield ~ yield))

# ~ a factor of 5 difference
# qubit = (nanodrop yield / 4.55) - 77.52 

# repeatability of the nano-drop

tmp <- nd_dat %>%
  .[-29,] %>%
  select(id, protocol, sex, measurement, yield) %>%
  mutate(id = paste0(id, "_", protocol)) %>%
  spread(measurement, yield)

names(tmp)[4:5] <- c("measure1", "measure2")

tmp %>%
  ggplot(aes(x = measure1, y = measure2)) +
  geom_point()+
  geom_smooth(method = "lm")

with(tmp, lm(measure1 ~ measure2)) %>% 
  summary




