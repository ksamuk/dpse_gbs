# the genomic distribution of restriction sites

############################################################
# libraries
############################################################

# install/load packages

#install.packages("SimRAD")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")

library("seqinr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("viridis")

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
ref_genome <- read.fasta("data/dpse-all-chromosome-r3.04.fasta.gz", as.string = TRUE)
ref_genome <- ref_genome[1:15]

# restriction enzyme recognition sites

Csp6I <- "gtac"
MspI <- "ccgg"
EcoRI <- "gaattc"

############################################################
# find restriction sites in the genome
############################################################


cut_sites <- list()

for (i in 1:length(ref_genome)){
  
  csp_sites <- gregexpr(Csp6I, ref_genome[i]) %>% unlist %>% data.frame(chr = names(ref_genome[i]), cut = "Csp6I", pos = .)
  msp_sites <- gregexpr(MspI, ref_genome[i]) %>% unlist %>% data.frame(chr = names(ref_genome[i]), cut = "MspI_EcoRI", pos = .)
  eco_sites <- gregexpr(EcoRI, ref_genome[i]) %>% unlist %>% data.frame(chr = names(ref_genome[i]), cut = "MspI_EcoRI", pos = .)
  
  cut_sites[[i]] <- do.call("rbind", list(csp_sites, msp_sites, eco_sites)) %>% data.frame
  
}

# bind all sites into a single df

cut_sites <- do.call("rbind", cut_sites) %>% data.frame

# visual sites in the genome
cut_sites %>%
  ggplot(aes(x = pos, xend = pos, y =  as.factor(cut), yend = (as.numeric(as.factor(cut))*2)-1)) +
  geom_segment(alpha = 0.01)+
  #geom_point(shape = "|")+
  facet_wrap(~chr, scales = "free_x")

cut_sites %>%
  #filter(chr == 2) %>%
  ggplot(aes(x = pos, y = as.factor(cut))) +
  geom_bin2d(bins = 500)+
  #geom_point(shape = "|")+
  facet_wrap(~chr) +
  scale_fill_viridis() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

