# track dna QC results

library("tidyverse")

# master dna file

dna_all <- read.csv("data/dna_extractions_master.csv", h = T, na.strings = c("NA", ""), stringsAsFactors = FALSE)
dna_all <- dna_all %>%
  mutate(unique_id = paste(population, female_id, sub_id, generation, sep = "-"))

# locate nanodrop results
nano_drop_files <- list.files("data/nanodrop", pattern = ".txt$", full.names = TRUE)

# read in nanodrop files

# static column names (nanodrop defaults have spaces -- invalid)
col_names <- c("sample_id", "user_id", "date", "time", "ng_ul", "a260", "a280", "r260/280", "r260/230", "constant", "cursor_pos", "cursor_abs", "a340_raw", "BLANK")

# load data
nd_dat <- lapply(nano_drop_files, read.table, h = T, row.names = NULL, skip = 1, col.names = col_names, stringsAsFactors = FALSE)

# bind to dataframe
nd_dat <- do.call("rbind", nd_dat)

# fix column offset
nd_dat$time <- paste0(nd_dat$time, nd_dat$ng_ul)
nd_dat[5:(ncol(nd_dat) - 1)] <- nd_dat[6:(ncol(nd_dat))]
nd_dat <- nd_dat[,-ncol(nd_dat)]

# clip dates 
nd_dat$extraction_date <- strsplit(nd_dat$sample_id, split = "-") %>% lapply(., function(x)x[5]) %>% unlist
nd_dat$sample_id <- gsub("-[a-z,A-z]*[0-9]{2}$", "", nd_dat$sample_id)

# intersect nanodrop data with full dataset

dna_all$tbd <- !(dna_all$unique_id %in% nd_dat$sample_id)

nd_dat %>%
  rename(unique_id = sample_id) %>%
  select(unique_id, ng_ul, r260.280, r260.230) %>%
  left_join(dna_all, .) %>%
  write.table("data/dna_qc_master.csv", row.names = FALSE, quote = FALSE)

# generate a list of TBD samples
dna_all %>%
  filter(tbd) %>%
  filter(generation == "F1") %>%
  select(unique_id) %>%
  write.table("data/qc_tbd.csv", row.names = FALSE, quote = FALSE)

dna_all %>%
  filter(!tbd) %>%
  filter(generation == "F1") %>%
  select(unique_id) %>%
  write.table("data/qc_done.csv", row.names = FALSE, quote = FALSE)

# generate a list of TBD samples
nd_dat %>%
  select(-user_id, -constant, -cursor_pos, -cursor_abs, -time) %>%
  write.table("data/nanodrop_results.csv", row.names = FALSE, quote = FALSE)



