library(dada2); packageVersion("dada2")
library(ggplot2)
################RUN 1 ###################################################################################################
# File parsing 
filtpathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/REV/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
#Plot & save errors 
Run1_errF<-plotErrors(errF, nominalQ=TRUE)
Run1_errR<-plotErrors(errR, nominalQ=TRUE)
ggsave(filename=paste(filtpathF, "run1_errorFwd.pdf", sep="/"), plot=Run1_errF)
ggsave(filename=paste(filtpathR, "run1_errorRev.pdf", sep="/"), plot=Run1_errR)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
dadaFs<-vector("list", length(sample.names))
dadaRs<-vector("list", length(sample.names))
names(mergers) <- sample.names
names(dadaFs)<-sample.names
names(dadaRs)<-sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  dadaFs[[sam]]<-ddF
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  dadaRs[[sam]]<-ddR
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
# Construct sequence table and remove chimeras
seqtab_run1 <- makeSequenceTable(mergers)
saveRDS(seqtab_run1, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/seqtab_run1.rds")# CHANGE ME to where you want sequence table saved
#Sanity check -- loss through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)
track_run1<-track
rm(derepF); rm(derepR)
print("Run1 Done!")
##################################END Run 1#################################################
################RUN 2 ###################################################################################################
# File parsing 
filtpathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/REV/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
#Plot & save errors 
Run2_errF<-plotErrors(errF, nominalQ=TRUE)
Run2_errR<-plotErrors(errR, nominalQ=TRUE)
ggsave(filename=paste(filtpathF, "run2_errorFwd.pdf", sep="/"), plot=Run2_errF)
ggsave(filename=paste(filtpathR, "run2_errorRev.pdf", sep="/"), plot=Run2_errR)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
dadaFs<-vector("list", length(sample.names))
dadaRs<-vector("list", length(sample.names))
names(mergers) <- sample.names
names(dadaFs)<-sample.names
names(dadaRs)<-sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  dadaFs[[sam]]<-ddF
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  dadaRs[[sam]]<-ddR
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
# Construct sequence table and remove chimeras
seqtab_run2 <- makeSequenceTable(mergers)
saveRDS(seqtab_run2, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/seqtab_run2.rds")# CHANGE ME to where you want sequence table saved
#sanity check - track loss 
#Sanity check -- loss through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)
track_run2<-track
rm(derepF); rm(derepR)
print("Run2 Done!")
##################################END Run 2#################################################

################RUN 3 ###################################################################################################
# File parsing 
filtpathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/FWD/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/REV/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
#Plot & save errors 
Run3_errF<-plotErrors(errF, nominalQ=TRUE)
Run3_errR<-plotErrors(errR, nominalQ=TRUE)
ggsave(filename=paste(filtpathF, "run3_errorFwd.pdf", sep="/"), plot=Run3_errF)
ggsave(filename=paste(filtpathR, "run3_errorRev.pdf", sep="/"), plot=Run3_errR)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
dadaFs<-vector("list", length(sample.names))
dadaRs<-vector("list", length(sample.names))
names(mergers) <- sample.names
names(dadaFs)<-sample.names
names(dadaRs)<-sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  dadaFs[[sam]]<-ddF
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  dadaRs[[sam]]<-ddR
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
# Construct sequence table and remove chimeras
seqtab_run3 <- makeSequenceTable(mergers)
saveRDS(seqtab_run3, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/seqtab_run3.rds")# CHANGE ME to where you want sequence table saved
#Sanity check -- loss through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)
track_run3<-track
rm(derepF); rm(derepR)
print("Run3 Done!")
##################################END Run 3#################################################
load(file="/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Results/dada2results/cav_out_filter_trim_dada2.Rdta")
CombinedStats<-rbind(track_run1, track_run2, track_run3)
TrackReads<-cbind(FilterAndTrimRunStats, CombinedStats)
saveRDS(TrackReads, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Results/dada2results/track_reads.rds")

