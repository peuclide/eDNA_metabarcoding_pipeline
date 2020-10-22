


library(dada2)

# path to fastq files
path <- "Y:/Peter/Projects/eDNA/eDNA_preliminary_analysis/12S_real_data/DST-raw-data-files/adapt_filt"

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.filt.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.filt.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

## filter and trim


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt2.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(18), maxLen = 125, minLen = 100,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

fnFsFilt <- sort(list.files(paste0(path, "/filtered"), pattern="_F_filt2.fastq", full.names = TRUE))
plotQualityProfile(fnFsFilt[1:4])


## with 0 files - remove these files

"183_S179_L001_R1_001.fastq.gz"	
"183_S179_L001_R1_001.filt.fastq.gz"



### Errors

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errR, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


## merge reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

## remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# write filter results
write.csv(track, "./eDNA_preliminary_analysis/12S_real_data/tracked_dada2_filters.csv")


### sequence tracking plots



## assign taxa
# ref <- "Y:/Peter/Projects/eDNA/eDNA_preliminary_analysis/12S_real_data/yer_reference/all_12s_GreatLakes.fasta"
# taxa <- assignTaxonomy(seqtab.nochim, ref, multithread=TRUE)


## export for local blast

write.table(t(seqtab.nochim), "Y:/Peter/Projects/eDNA/eDNA_preliminary_analysis/12S_real_data/seqtab-nochim_primer_filt_12s_7_7.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


View(t(seqtab.nochim))


seq <- rownames(t(seqtab.nochim))
seq_ID <- paste0("MOTU_", 1:length(seq))


### make fasta of ASV

make_fasta <- function(seq,s_ID){
  m <- matrix(nrow=length(seq)*2, ncol =1)
  name_pos <- seq(1,(nrow(m)-1),2)
  seq_pos <- seq(2,nrow(m),2)
  m[name_pos,] <- paste0(">",seq_ID)
  m[seq_pos,] <- seq
  m
  
}

fasta <- make_fasta(seq, seq_ID)
View(fasta)

write.table(fasta, "Y:/Peter/Projects/eDNA/eDNA_preliminary_analysis/12S_real_data/seqtab-nochim_primer_filt_12s_7_7.fasta", row.names=F, col.names=F, quote=FALSE) 
            