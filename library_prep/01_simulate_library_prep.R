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

# create list of fasta sequences from reference genome
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
MspI <- list(cs_5p1 = "C", cs_3p1 = "CGG", cs_5p2 = "C", cs_3p2 = "GGC", name = "MspI")

# EcoRI
# 5' ---G   AATTC--- 3'
# 3' ---CTTAA   G--- 5'
EcoRI <- list(cs_5p1 = "G", cs_3p1 = "AATTC", cs_5p2 = "G", cs_3p2 = "CTTA", name = "EcoRI")

############################################################
# Simulating distribution of reads per inidivudal
############################################################

# function for full library prep simulation

reference_genome = ref_genome
enzyme_cut1 = Csp6I
size_range = c(300, 500) 
num_ind = 60 
expected_reads = 100000000
expected_sequencing_variance = 0.4
read_cutoff = 10 
ind_cutoff = 0.8
n_top_cuts = 1000
enzyme_cut2 = NULL

rr_library_prep <- function(reference_genome, enzyme_cut1, enzyme_cut2 = NULL, size_range = c(300, 500), num_ind = 60, 
														expected_reads = 100000000, expected_sequencing_variance = 0.4, 
														read_cutoff = 10, ind_cutoff = 0.8, n_top_cuts = 1000){
	
	# perform in silico restriction digest 
	digest <- insilico.digest(reference_genome, enzyme_cut1[1], enzyme_cut1[2], verbose = FALSE)
	
	# perform second restriction digest (if specified)
	if (!is.null(enzyme_cut2)){
		digest <- insilico.digest(digest, enzyme_cut2[1], enzyme_cut2[2], verbose = FALSE)
		enzyme_name <- paste0(enzyme_cut1$name, "_", enzyme_cut2$name)
	} else{
		enzyme_name <- enzyme_cut1$name
	}
	
	# size selection of fragments
	size_sel <- size.select(digest, size_range[1], size_range[2], verbose = FALSE, graph = FALSE)
	
	# simulate variable sequencing of each fragment
	# assumes an approximate exponential distribution of sequencing depth per fragment
	# tbd : specify fragment sequecning variance function
	all_reads <- rexp(length(size_sel), 1/(expected_reads/length(size_sel)))
	
	# init output df list
	out_df <- list()
	
	# if num_ind is a vector of length > 1
	# iterate over number of individuals and same summary stats
	
	for (i in 1:length(num_ind)){
		
	# expected reads for each individual
	expected_reads_per_ind <- expected_reads/num_ind[i]
	
	# simulate variable sequencing per individual
	# assumes an approximate gaussian distribution of sequencing per individual
	# tbd: specify individual sequencing variance function
	prop_reads <- rnorm(num_ind[i], mean = expected_reads_per_ind, sd = expected_reads_per_ind * expected_sequencing_variance) 
	prop_reads <- (prop_reads / (expected_reads_per_ind)) / num_ind[i]
	prop_reads <- ifelse(prop_reads < 0, 0, prop_reads)
	
	# assign reads to individuals
	assigned_reads <- lapply(prop_reads, function(x) round(all_reads * x)) %>% data.frame
	
	# create data frame of assigned reads
	names(assigned_reads) <- paste0("ind_", 1:num_ind[i])
	assigned_reads <- data.frame(fragment = 1:length(all_reads), assigned_reads) %>%
		gather(-fragment, key = "ind", value = "frag_count")
	
	# calculate the number of fragments that pass:
	# 1: the depth threshold (reads per fragment)
	# 2: the individual representation threshold (prop individuals that pass #1)
	prop_fragments_useable <- assigned_reads %>%
		mutate(frag_count_acceptable = frag_count > read_cutoff) %>%
		group_by(fragment) %>%
		summarise(prop_ind_acceptable = mean(frag_count_acceptable )) %>%
		mutate(ind_count_acceptable = prop_ind_acceptable > ind_cutoff) %>%
		summarise(prop_frags_usable = mean(ind_count_acceptable)) %>%
		as.numeric
	
	# depth of the top 1000 fragments
	top_n <- assigned_reads %>%
		group_by(fragment) %>%
		summarise(mean_depth = mean(frag_count)) %>%
		arrange(desc(mean_depth)) %>%
		select(mean_depth) %>%
		.[1:n_top_cuts,] %>%
		unlist
	
	# generate a single-row dataframe as output
	out_df[[i]] <- data.frame(enzyme_name = enzyme_name, num_ind = num_ind[i], number_of_fragments = length(size_sel), 
						 mean_reads_per_fragment_per_ind = (expected_reads/length(size_sel))/num_ind[i], prop_fragments_useable = prop_fragments_useable, 
						 number_fragments_usable = prop_fragments_useable*length(size_sel), mean_top_n_fragments = round(mean(top_n)),
						 min_top_n_fragments = round(min(top_n)), max_top_n_fragments = round(max(top_n)))
	
	}
	
	out_df
	
}

# perform in silico preps for csp6I and ecoRI + mspI
# varying number of individuals

ind_list <- seq(20, 600, by = 20)

csp_preps <- rr_library_prep(num_ind = ind_list, reference_genome = ref_genome, enzyme_cut1 = Csp6I)
eco_msp_preps <- rr_library_prep(num_ind = ind_list, reference_genome = ref_genome, enzyme_cut1 = MspI, enzyme_cut2 = EcoRI)
eco_preps <- rr_library_prep(num_ind = ind_list, reference_genome = ref_genome, enzyme_cut1 = EcoRI)


prep_df <- bind_rows(csp_preps, eco_msp_preps, eco_preps)

# bar plot of usable fragments, split by enzyme
prep_df %>%
	filter(num_ind < 600) %>%
	ggplot(aes(x = num_ind, y = number_fragments_usable)) +
	geom_bar(stat = "identity") +
	facet_wrap(~enzyme_name, scales = "free_y")+
	theme_bw()+
	xlab("Number of individuals")+
	ylab("Number of usable restriction fragments")

# line plot comparing usable fragments in both digests 
prep_df %>%
	filter(num_ind < 600) %>%
	ggplot(aes(x = num_ind, y = number_fragments_usable, color = enzyme_name)) +
	geom_line(size = 2)+
	theme_bw()+
	theme(legend.justification = c(1, 1), 
				legend.position = c(1, 1),
				legend.background = element_blank(),
				legend.key = element_blank())+
	xlab("Number of individuals")+
	ylab("Number of usable restriction fragments")

prep_df %>%
	filter(num_ind < 400) %>%
	ggplot(aes(x = num_ind, y = mean_top_n_fragments, color = enzyme_name, ymin = min_top_n_fragments, ymax = max_top_n_fragments)) +
	geom_point(size = 2)+
	geom_errorbar()+
	geom_hline(yintercept = 20, lty = 2)+
	theme_bw()+
	theme(legend.position = "none")+
	facet_wrap(~enzyme_name, scales = "free_y")+
	xlab("Number of individuals")+
	ylab("Depth of top 1000 fragments")




