#updated 3/11
rm(list=ls())
setwd("Y:/Peter/Projects/eDNA")
library(dplyr)
#parameters for trimming
start_motu=1
end_motu=267
#make sure end motu is full count from seq table, NOT .out file
min_perc_id=98
min_align_len=95
hits_per_motu_to_return=1
hits_per_motu_to_return_sim_logic=10
align_perc_dif=2
max_mistmatches=2
#different approach trying to output multiple alignments if alignment is ambiguous, start with if align diff is > 2% only print out first, if <2%, output every alignment within 2% 
use_similarity_logic="TRUE"


## blast output from blastn set to -outfmt 6

blast_output=read.table(file="./results_primerFilt_Yer_Sard_Becky_peter.out", sep="\t")
colnames(blast_output)=c("query","match","perc_ID","align_len","mismatch","gap","query_start","query_stop",
                         "match_start","match_stop","eval","bit")
filtered_subset_dat=NULL
for(i in start_motu:end_motu)
{
  current_motu=paste("MOTU",i,sep="_")
  current_dat=filter(blast_output, query==current_motu)
  no_hit_line=cbind(current_motu,"no_hits","no_hits","no_hits","no_hits","no_hits","no_hits","no_hits","no_hits","no_hits","no_hits","no_hits")
  colnames(no_hit_line)=c("query","match","perc_ID","align_len","mismatch","gap","query_start","query_stop",
                          "match_start","match_stop","eval","bit")
  if(nrow(current_dat)==0){filtered_subset_dat=rbind(filtered_subset_dat, no_hit_line); next}
  filtered_current_dat=filter(current_dat, perc_ID>=min_perc_id, align_len>=min_align_len)
  if(nrow(filtered_current_dat)==0){filtered_subset_dat=rbind(filtered_subset_dat, no_hit_line);next}
  ordered_dat=arrange(filtered_current_dat, desc(perc_ID))
  if(use_similarity_logic != "TRUE")
    {
      if(nrow(ordered_dat)<hits_per_motu_to_return){filtered_subset_dat=rbind(filtered_subset_dat, ordered_dat)}
      else{filtered_subset_dat=rbind(filtered_subset_dat, ordered_dat[1:hits_per_motu_to_return,])}
    }
  if(use_similarity_logic == "TRUE")
  {
    #if only 1 matches within the parameters then output that one match
    if(nrow(filtered_current_dat)==1){filtered_subset_dat=rbind(filtered_subset_dat, ordered_dat);next}
    top_hit=ordered_dat[1,]
    second_hit=ordered_dat[2,]
    if(second_hit[,5]-top_hit[,5] >= max_mistmatches){filtered_subset_dat=rbind(filtered_subset_dat, top_hit)}
    #print all hits within align_perc_dif of the top hit
    else if(second_hit[,5]-top_hit[,5] < max_mistmatches)
    {
      min_hit_perc=top_hit[,3]-align_perc_dif
      hits_close_to_top=filter(ordered_dat, perc_ID>=min_hit_perc)
      if(nrow(hits_close_to_top)<hits_per_motu_to_return_sim_logic){filtered_subset_dat=rbind(filtered_subset_dat, hits_close_to_top)}
      else{filtered_subset_dat=rbind(filtered_subset_dat, hits_close_to_top[1:hits_per_motu_to_return_sim_logic,])}
    }

  }
}

write.csv(filtered_subset_dat, file="./12s_pooled_aligned_with_unks_primerFilt_test2.csv", row.names=F)









