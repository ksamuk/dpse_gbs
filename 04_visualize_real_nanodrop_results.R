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

# nanodrop files

nanodrop_files <- list.files("data/nanodrop", pattern = ".txt$", full.names = TRUE)

# static column names (nanodrop defaults have spaces -- invalid)
col_names <- c("sample_id", "user_id", "date", "time", "ng_ul", "a260", "a280", "r260/280", "r260/230", "constant", "cursor_pos", "cursor_abs", "a340_raw", "BLANK")

# load data
nd_dat <- lapply(nanodrop_files, read.table, h = T, row.names = NULL, skip = 1, col.names = col_names)

# bind to dataframe
nd_dat <- do.call("rbind", nd_dat)

# fix column offset
nd_dat$time <- paste0(nd_dat$time, nd_dat$ng_ul)
nd_dat[5:(ncol(nd_dat) - 1)] <- nd_dat[6:(ncol(nd_dat))]
nd_dat <- nd_dat[,-ncol(nd_dat)]

# labels for facets
stat_names <- c(
	`r_260.280` = "r260/280",
	`r_260.230` = "r260/230",
	`ng/ul` = "Total Yield (ng)"
)

# plot results of DNA extractions
nd_dat %>%
	ggplot(aes(x = date, y = ng_ul))+
	geom_jitter()+
	theme_bw(base_size = 18)+
	theme(legend.position = "none",
				strip.background = element_blank(),
				axis.title.x = element_text(size = 16, face = "bold"),
				strip.text = element_text(face = "bold"))

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




