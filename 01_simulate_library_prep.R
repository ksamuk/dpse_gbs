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
library("tidyr")
library("ggplot2")

############################################################
# raw data
############################################################

# download reference genome
genome_url <- "ftp://ftp.flybase.net/genomes/Drosophila_pseudoobscura/dpse_r3.04_FB2016_02/fasta/dpse-all-chromosome-r3.04.fasta.gz"
dir.create("data")

if(!file.exists("data/dpse-all-chromosome-r3.04.fasta.gz")){
	download.file(genome_url, "data/dpse-all-chromosome-r3.04.fasta.gz")
}

ref_genome <- ref.DNAseq("data/dpse-all-chromosome-r3.04.fasta.gz", subselect.contigs = TRUE, prop.contigs = 1.0)


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

# function for full library prep simulation

# dummies for debug
reference_genome = ref_genome
enzyme_cut = MspI
size_range = c(300, 500)
num_ind = 60 
expected_reads = 75000000 
expected_sequencing_variance = 0.4 
read_cutoff = 10

rr_library_prep <- function(reference_genome, enzyme_cut, size_range = c(300, 500), num_ind = 60, 
														expected_reads = 75000000, expected_sequencing_variance = 0.4, 
														read_cutoff = 10, ind_cutoff = 0.8){
	
	# perform in silico digest 
	digest <- insilico.digest(reference_genome, enzyme_cut[1], enzyme_cut[2], verbose = FALSE)
	
	# size selection
	size_sel <- size.select(digest, size_range[1], size_range[2], verbose = FALSE, graph = FALSE)
	
	# simulate variable sequencing of each fragment
	all_reads <- rexp(length(size_sel), 1/(expected_reads/length(size_sel)))
	
	# expected reads for each individual
	expected_reads_per_ind <- expected_reads/num_ind
	
	# simulate variable sequencing per individual
	prop_reads <- rnorm(num_ind, mean = expected_reads_per_ind, sd = expected_reads_per_ind * expected_sequencing_variance) 
	prop_reads <- (prop_reads / (expected_reads_per_ind)) / num_ind
	prop_reads <- ifelse(prop_reads < 0, 0, prop_reads)
	
	# assign reads to individuals
	assigned_reads <- lapply(prop_reads, function(x) round(all_reads * x)) %>% data.frame
	
	# create data frame of assigned reads
	names(assigned_reads) <- paste0("ind_", 1:num_ind)
	assigned_reads <- data.frame(fragment = 1:length(all_reads), assigned_reads) %>%
		gather(-fragment, key = "ind", value = "frag_count")
	
	prop_fragments_useable <- assigned_reads %>%
		mutate(frag_count_acceptable = frag_count > read_cutoff) %>%
		group_by(fragment) %>%
		summarise(prop_ind_acceptable = mean(frag_count_acceptable )) %>%
		mutate(ind_count_acceptable = prop_ind_acceptable > ind_cutoff) %>%
		summarise(prop_frags_usable = mean(ind_count_acceptable)) %>%
		as.numeric


	# test if average reads per frag falls below threshold
	prop_average_reads_below_threshold <- sum(average_reads_per_frag < read_cutoff)/length(average_reads_per_frag)
	
	data.frame(enzyme_name = enzyme_cut$name, num_ind = num_ind, number_of_fragments = length(size_sel), mean_reads_per_fragment = mean_reads_per_frag, prop_below_threshold = prop_average_reads_below_threshold)
	
}

ind_list <- list(60, 80, 100, 120, 140)
csp_preps <- lapply(ind_list, function(x) rr_library_prep(num_ind = x, reference_genome = ref_genome, enzyme_cut = Csp6I))

csp_df <- rbind_all(csp_preps)

csp_df %>%
	ggplot(aes(x = num_ind, y = prop_below_threshold)) +
	geom_line() +
	theme_bw()




