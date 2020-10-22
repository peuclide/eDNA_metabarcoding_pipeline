### apply control thresholds

library(ggridges)
library(knitr)
library(tidyverse)

## taxa names

taxa <- read.csv("taxonomy.csv")
taxa_matrix <- read.csv("taxa_matrix_7-17-20.csv")
taxa2 <- taxa %>% filter(sci_name %in% colnames(taxa_matrix)) %>% group_by(sci_name) %>% 
  summarise(com_name=first(com_name), family=first(family), genus=first(genus), fish = first(fish))


# identify all blanks/control samples in dataset - there should be field blanks corresponding to blanks taken along samples in the field, mock communities which were generated from known quantities of DNA, and NTC which are PCR negative controls. 

toMatch_blanks <- c("BLK", "BLANK", "blank","blk", "Blank", "Blk", "FB", "fb")
toMatch_mock <- c("MOCK", "mock", "Mock")
toMatch_NTC <- c("NTC", "ntc")
toMatch_ExtNeg <- c("ExtNeg", "extneg")


field_blank <- grep(paste(toMatch_blanks, collapse = "|" ), taxa_matrix$SampleName, value = T)
mock_com <- grep(paste(toMatch_mock, collapse = "|" ), taxa_matrix$SampleName, value = T)
NTC <- grep(paste(toMatch_NTC, collapse = "|" ), taxa_matrix$SampleName, value = T)
ExtNeg <- grep(paste(toMatch_ExtNeg, collapse = "|" ), taxa_matrix$SampleName, value = T)

# create long file of just sample data
sample_dat_long <- taxa_matrix %>% filter(!(SampleName %in% c(field_blank, mock_com, NTC, ExtNeg))) %>% 
  gather(.,"species", "reads", 2:(which(colnames(taxa_matrix)=="SampleName")-1)) %>% 
  left_join(taxa2, by = c("species"="sci_name"))


# create long file of just control data
control_dat_long <- taxa_matrix %>% filter(SampleName %in% c(field_blank, mock_com, NTC, ExtNeg)) %>% 
  gather(.,"species", "reads", 2:(which(colnames(taxa_matrix)=="SampleName")-1)) %>% 
  left_join(taxa2, by =c("species"="sci_name"))

# thresh_1 filter - all lab blanks

Thresh_1 <- control_dat_long %>% filter(SampleName %in% c(ExtNeg, NTC), fish == 1) %>% 
  group_by(SampleName) %>%
  #filter(reads>0) %>% 
  summarise(sumReads = round(sum(reads))) %>% 
  summarise(Thresh_1 = round(mean(sumReads)))
Thresh_1 <- as.numeric(Thresh_1[1])
# Thresh_2 - lake blanks
Thresh_2 <- control_dat_long %>% filter(SampleName %in% field_blank, fish == 1) %>% 
  group_by(Lake, SampleName) %>% 
  summarise(sumReads = round(sum(reads))) %>% 
  #filter(reads>0) %>% 
  group_by(Lake) %>% 
  summarise(Thresh_2 = round(mean(sumReads)))

# Thresh 1 and 2 combination
sample_dat_long <- sample_dat_long %>% mutate(Thresh_1 = Thresh_1)
sample_dat_long <- sample_dat_long %>% left_join(., Thresh_2, by = c("Lake"))

# convert NA to 0
sample_dat_long$Thresh_1[is.na(sample_dat_long$Thresh_1)] <- 0
sample_dat_long$Thresh_2[is.na(sample_dat_long$Thresh_2)] <- 0



# normalize reads
norm_reads = c()
thresh=c()

for(i in 1:nrow(sample_dat_long)){
  tmp <- sample_dat_long %>% filter(Lake == sample_dat_long$Lake[i], species == sample_dat_long$species[i])
  if(sum(tmp$reads>0)>1) {
    norm_reads[i] = sample_dat_long$reads[i]
    thresh[i]=1
  }else if(sample_dat_long$reads[i]<sample_dat_long$Thresh_1[i]){
    norm_reads = append(norm_reads, 0)
    thresh[i]=2
  }else if(sample_dat_long$reads[i]<=sample_dat_long$Thresh_2[i]){
    norm_reads[i] = 0
    thresh[i]=3
  }else {
    norm_reads[i] = sample_dat_long$reads[i]
    thresh[i]=4}
}

sample_dat_long$norm_reads <- norm_reads
sample_dat_long$thresh <- thresh

sample_dat_long <- sample_dat_long %>% mutate(detection = sample_dat_long$norm_reads>0)


## Optional write long form matrix
#write.csv(sample_dat_long, "./eDNA_preliminary_analysis/12S_real_data/primer_filt_yer_sard_becky_fullRef/taxa_matrix_7-23-20.norm.long.csv")

# Make detection matrix

sample_dat_long_det <- sample_dat_long %>% filter(fish ==T) %>% select("Lake", "species", "family", "genus","detection") %>% group_by(species) %>% 
  filter(sum(detection)>0) %>% group_by(Lake, species) %>% summarise(det =sum(detection)) %>% pivot_wider(names_from="Lake", values_from="det")



sample_dat_long_det2 <- left_join(sample_dat_long_det,taxa2, by= c("species"="sci_name"))

# write taxa matrix
write.csv(sample_dat_long_det2[,c(1,11:13,2:10)], "./EDNA_TAXA_NORM.csv", row.names = F)


#### detection summary

detections=sample_dat_long %>% filter(fish ==1,detection ==T) %>% group_by(Lake, SampleName)

# number of smaples with positive detections
length(unique(detections$SampleName))
length(unique(detections$species))
table(detections$thresh)


## fialed detections
failed_det <- sample_dat_long %>% filter(detection==F, fish ==1, reads >0)


table(failed_det$thresh)

