library(tidyverse)
library(taxize)


### combine DADA2 table with MOTUs

## file from 1_raw2dada2.R
DADA2_tab <- read.table("./seqtab-nochim_primer_filt_12s_7_7.txt")

# Table constructed through Local blast of ASVs
BLASTN_tab <- read.csv("./12s_pooled_aligned_with_unks_primerFilt_test2.csv")

# codes for each site name
SITE_CODES <- read.table("./metadata_file_validated.tsv", sep = "\t", header =T)
SITE_CODES$SampleID <- paste0("X", SITE_CODES$SampleID)

# Add MOTU info to DADA2 output
DADA2_tab <- DADA2_tab %>% mutate(MOTU=paste0("MOTU_", seq(1:nrow(DADA2_tab))))

# FIlter to one entry per MOTU for BLASTN results
BLASTN_taxa <- BLASTN_tab %>% group_by(query) %>%summarise (taxa = first(taxa))

##combine spp name to dada 2 tab
DADA2_tab <- left_join(DADA2_tab, BLASTN_taxa, by = c("MOTU"="query"))



read_table <- DADA2_tab %>% select(-"MOTU") %>% group_by(taxa) %>% summarise_each(funs(sum))
read_table <- t(read_table) 
colnames(read_table) <- read_table[1,]
read_table <- read_table[-1,]

read_table <- data.frame(read_table) %>% 
  rownames_to_column(var = "SampleID") %>% 
  left_join(SITE_CODES[,1:3], by = "SampleID") %>%
  select(-"V1")



sci_names=data.frame(colnames(read_table)[3:75])

list_of_common_names=NULL
for(i in 1:nrow(sci_names))
{
  list_of_common_names=rbind(list_of_common_names, sci2comm(scinames=as.character(sci_names[i,1])))
}

sci_and_common_names=cbind(as.vector(sci_names[,1]),as.vector(list_of_common_names))

sci2comm(colnames(read_table))


## output taxa matrix
write.csv(read_table, "./eDNA_preliminary_analysis/12S_real_data/primer_filt_yer_sard_becky_fullRef/taxa_matrix_7-17-20.csv", row.names =F)


