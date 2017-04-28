library("tidyverse")

# read in F1 data
qc_dat <- read.csv("data/dna_qc_master_F1.csv", header = T)

head(qc_dat)

# tally the number of good extractions from all populations
qc_dat %>%
  select(-r260.280, -r260.230) %>%
  filter(!is.na(ng_ul)) %>%
  group_by(population, female_id, unique_id) %>%
  summarise(ng_ul = mean(ng_ul)) %>%
  group_by(population, female_id) %>%
  mutate(passing = as.numeric(ng_ul > 50)) %>%
  tally(passing) %>%
  View

# build list of "gold star" samples
qc_dat %>%
  select(-r260.280, -r260.230) %>%
  filter(!is.na(ng_ul)) %>%
  group_by(population, female_id, unique_id) %>%
  summarise(ng_ul = mean(ng_ul)) %>%
  ungroup %>%
  filter(ng_ul > 50) %>%
  arrange(population, female_id, desc(ng_ul)) %>% 
  group_by(population, female_id) %>%
  mutate(conc_rank = 1:length(ng_ul)) %>%
  filter(conc_rank <= 2) %>%
  write.table("data/gold_star_samples.csv", row.names = FALSE, quote = FALSE)