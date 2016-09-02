# simulate a GBS library prep using the Drosophila pseudoobscura reference genome
# examine resultant coverage/fragment distribution for different restriction enzymes
# kms sept 2016

############################################################
# libraries
############################################################

# install/load packages
# install.packages("SimRAD")
library("SimRAD")

############################################################
# raw data
############################################################

# download reference genome
genome_url <- "ftp://ftp.flybase.net/genomes/Drosophila_pseudoobscura/dpse_r3.04_FB2016_02/fasta/dpse-all-chromosome-r3.04.fasta.gz"

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
# SimRAD analysis
############################################################

# load reference genome
ref_genome <- ref.DNAseq("data/dpse-all-chromosome-r3.04.fasta.gz", prop.contigs = 1.0)

# in silico digests

# PstI
# 65482 sites
digest_pst <- insilico.digest(ref_genome, PstI[1], PstI[2], verbose=TRUE)

# MspI
# 316549 sites
digest_msp <- insilico.digest(ref_genome, MspI[1], MspI[2], verbose=TRUE)

# EcoRI 
# 35225 sites
digest_eco <- insilico.digest(ref_genome, EcoRI[1], EcoRI[2], verbose=TRUE)

# EcoRI + MspI
# sites
digest_eco_msp <- insilico.digest(digest_eco, MspI[1], MspI[2], verbose=TRUE)

# size selection

tmp <- size.select(digest_msp, 300, 600)




cap6 
pstI
mspI
ecorI