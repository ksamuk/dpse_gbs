# check representation of final dna dataset

library("tidyverse")

# read in data

dat <- read.csv("data/qubit_fillout_final.csv")

dat %>%
  group_by(population, type) %>%
  tally

# generate a vector of plate positions

row_names <- LETTERS[1:8]
col_names <- 1:12

# create cell names and plate numbers
cell_names <- lapply(row_names, function(x) paste0(x, col_names)) %>% unlist
plate_num <- c(rep(1, 96), rep(2, 96))

plate_positions <- data.frame(cell_name = rep(cell_names, 2), plate_num)

# subset to number of actual samples
plate_positions <- plate_positions[1:nrow(dat),]
plate_position_vec <- paste0(plate_positions$plate_num, "_", plate_positions$cell_name)

plate_rand <- sample(plate_position_vec)

# assign random plate positiosn to DNA samples
plate_rand <- data.frame(dat, cell_name = plate_rand)

plate_rand <- plate_rand %>%
  mutate(cell_name = as.character(cell_name)) %>% 
  mutate(plate_num = gsub("_.*", "", cell_name)) %>%
  mutate(row_name = gsub("[12]_|[_0-9]", "", cell_name)) %>%
  mutate(col_name = gsub("[12]_|[A-Z]", "", cell_name))

plate_rand <- plate_rand %>%
  arrange(plate_num, row_name, as.numeric(col_name))

write_csv(plate_rand, path = "data/qubit_fillout_final_positions.csv")
