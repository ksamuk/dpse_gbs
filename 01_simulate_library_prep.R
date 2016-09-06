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

# PstI 
# 5'--CTGCA  G--3'
# 3'--G  ACGTC—5'
PstI <- list(cs_5p1 = "CTGCA", cs_3p1 = "G", cs_5p2 = "ACGTC", cs_3p2  = "G")

# MspI
# 5' ---C   CGG--- 3'
# 3' ---GGC   C--- 5'
MspI <- list(cs_5p1 = "C", cs_3p1 = "CGG", cs_5p2 = "C", cs_3p2 = "GGC")

# EcoRI
# 5' ---G   AATTC--- 3'
# 3' ---CTTAA   G--- 5'
EcoRI <- list(cs_5p1 = "G", cs_3p1 = "AATTC", cs_5p2 = "G", cs_3p2 = "CTTA")

############################################################
# Basic SimRAD analysis
############################################################

# load reference genome
ref_genome <- ref.DNAseq("data/dpse-all-chromosome-r3.04.fasta.gz", prop.contigs = 1.0)

# in silico digests

# PstI
# 65482 fragments
digest_pst <- insilico.digest(ref_genome, PstI[1], PstI[2], verbose = TRUE)

# MspI
# 316549 fragments
digest_msp <- insilico.digest(ref_genome, MspI[1], MspI[2], verbose = TRUE)

# EcoRI 
# 35225 fragments
digest_eco <- insilico.digest(ref_genome, EcoRI[1], EcoRI[2], verbose = TRUE)

# EcoRI + MspI
# 351490 sites
digest_eco_msp <- insilico.digest(digest_eco, MspI[1], MspI[2], verbose = TRUE)

# size selection

pst_size <- size.select(digest_pst, 300, 500) # 5573 fragments between 300 and 500 bp 
msp_size <- size.select(digest_msp, 300, 500) # 39902 fragments between 300 and 500 bp 
eco_size <- size.select(digest_eco, 300, 500) # 1573 fragments between 300 and 500 bp 
eco_msp_size <- size.select(digest_eco_msp, 300, 500) # 47426 fragments between 300 and 500 bp 

############################################################
# Simulating distribution of reads per inidivudal
############################################################

# estimate variance in sequencing/individual using a real GBS run
gbs_test <- read.table("data/gbs_2015_lane1.txt")

sd_gbs <- gbs_test[,4] %>% 
	gsub("M", "", .) %>% 
	as.numeric %>%
	sd(na.rm = TRUE)

mean_gbs <- gbs_test[,4] %>% 
	gsub("M", "", .) %>% 
	as.numeric %>%
	mean(na.rm = TRUE)

coeff_var_gbs <- mean_gbs/sd_gbs

# illumina claims 125M for illumina 2500, but probably super variable per facility
total_sequenced_reads <- 75000000

# number of individuals
num_ind <- 60

mean_reads <- total_sequenced_reads/num_ind

# prop reads per individual
prop_reads <- rnorm(60, mean = reads_per_ind, sd = reads_per_ind/coeff_var_gbs)

hist(prop_reads)

prop_reads %>%
	ggplot(aes(x = prop_reads)) +
	geom_histogram(bins = 15)


# reads per individual






