# simulate a GBS library prep using the Drosophila pseudoobscura reference genome
# examine resultant coverage/fragment distribution for different restriction enzymes
# kms sept 2016

############################################################
# libraries
############################################################

# install/load packages

 #install.packages("SimRAD")
 #source("https://bioconductor.org/biocLite.R")
 #biocLite("ShortRead")

library("SimRAD")
library("dplyr")
library("zoo")
library("ggplot2")

############################################################
# raw data
############################################################

# download reference genome
genome_url <- "ftp://ftp.flybase.net/genomes/Drosophila_pseudoobscura/dpse_r3.04_FB2016_02/fasta/dpse-all-chromosome-r3.04.fasta.gz"
dir.create("data")
download.file(genome_url, "data/dpse-all-chromosome-r3.04.fasta.gz")

# restriction enzyme sequences

# enzyme cut sites from:
# https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites

# enzyme cut site variable names are coded as follows: 
# 5'--cs_5p1  cs_3p1--3'
# 3'--cs_5p2  cs_3p2-—5'

# Csp6I
#5' ---G   TAC--- 3'
#3' ---CAT   G--- 5'
Csp6I <- list(cs_5p1 = "G", cs_3p1 = "TAC", cs_5p2 = "G", cs_3p2  = "CAT", name = "Csp6I")
	
# MspI
# 5' ---C   CGG--- 3'
# 3' ---GGC   C--- 5'
MspI <- list(cs_5p1 = "C", cs_3p1 = "CGG", cs_5p2 = "C", cs_3p2 = "GGC", name = "Csp6I")

# EcoRI
# 5' ---G   AATTC--- 3'
# 3' ---CTTAA   G--- 5'
EcoRI <- list(cs_5p1 = "G", cs_3p1 = "AATTC", cs_5p2 = "G", cs_3p2 = "CTTA", name = "Csp6I")

############################################################
# Simulating distribution of reads per inidivudal
############################################################

# estimate variance in sequencing/individual using a real GBS run
gbs_test <- read.table("data/gbs_2015_lane1.txt")

hist(gbs_test[,4] )

sd_gbs <- gbs_test[,4] %>% 
	gsub("M", "", .) %>% 
	as.numeric %>%
	length
mean_gbs <- gbs_test[,4] %>% 
	gsub("M", "", .) %>% 
	as.numeric %>%
	mean(na.rm = TRUE)

coeff_var_gbs <- mean_gbs/sd_gbs # sd ~ 0.4 * mean ("sequencing variance")


# function for full library prep simulation

rr_library_prep<- function(reference_genome, enzyme_cut, size_range = c(300, 500), num_ind = 60, 
													expected_reads = 75000000, expected_sequencing_variance = 0.4, plot = TRUE,
													read_cutoff = 10){
	
	# perform in silico digest 
	digest <- insilico.digest(reference_genome, enzyme_cut[1], enzyme_cut[2], verbose = FALSE)
	
	# size selection
	size_sel <- size.select(digest, size_range[1], size_range[2], verbose = FALSE, graph = FALSE)
	
	# simulate variable sequencing
	prop_reads <- rnorm(num_ind, mean = 100, sd = 10 * expected_sequencing_variance)
	
	#
	prop_reads <- rnorm(num_ind, mean = expected_reads, sd = expected_reads*expected_sequencing_variance)

	# transform to average reads per fragment (assuming uniform fragment representation)
	average_reads_per_frag <- prop_reads/length(size_sel)
	mean_reads_per_frag <- mean(prop_reads/length(size_sel))
	
	ifelse(reads_per_frag < 0, 0, reads_per_frag)
	
	# test if average reads per frag falls below threshold
	prop_average_reads_below_threshold <- sum(average_reads_per_frag < read_cutoff)/length(average_reads_per_frag)
	
#	data.frame(enzyme_name = enzyme_cut$name, number_of_fragments = length(size_sel),average_reads_per_frag = average_reads_per_frag, mean_reads_per_fragment = mean_reads_per_frag, prop_average_reads_below_threshold = prop_average_reads_below_threshold)
	
	data.frame(enzyme_name = enzyme_cut$name, num_ind = num_ind, number_of_fragments = length(size_sel), mean_reads_per_fragment = mean_reads_per_frag, prop_below_threshold = prop_average_reads_below_threshold)
	
}

ind_list <- list(60, 80, 100, 120, 140)
csp_preps <- lapply(ind_list, function(x) rr_library_prep(num_ind = x, reference_genome = ref_genome, enzyme_cut = Csp6I))

csp_df <- rbind_all(csp_preps)

csp_df %>%
	ggplot(aes(x = num_ind, y = prop_below_threshold)) +
	geom_line() +
	theme_bw()

CspI_lib_prep <- rr_library_prep(ref_genome, enzyme_cut = Csp6I, size_range = c(300, 500), num_ind = 60, 
																 expected_reads = 75000000, expected_sequencing_variance = 0.4, plot = TRUE,
																 read_cutoff = 10)

MspI_lib_prep <- rr_library_prep(ref_genome, enzyme_cut = MspI, size_range = c(300, 500), num_ind = 60, 
																 expected_reads = 75000000, expected_sequencing_variance = 0.4, plot = TRUE,
																 read_cutoff = 10)

EcoRI_lib_prep <- rr_library_prep(ref_genome, enzyme_cut = EcoRI, size_range = c(300, 500), num_ind = 60, 
																 expected_reads = 75000000, expected_sequencing_variance = 0.4, plot = TRUE,
																 read_cutoff = 10)

rbind_all 


