library("dplyr")

dna_dat <- read.csv("data/extractions_tbd.csv")

source_female <- rep(dna_dat$source_female, each = 4)
population <- rep(dna_dat$population, each = 4)

new_dat <- data.frame(population, source_female)

left_join(new_dat, dna_dat)

all_dat <- read.csv("data/extractions_all.csv")

done_dat <- read.csv("data/extractions_done.csv")

tbd_dat <- left_join(all_dat, done_dat)

head(tbd_dat)

tbd_dat <- tbd_dat %>%
  filter(is.na(extracted)) %>%
  select(-pop, -date, -extracted)

head(tbd_dat)

write.csv(tbd_dat, "data/extractions_tbd.csv")
  