library(dada2); packageVersion("dada2")

################RUN 1 ###################################################################################################
# File parsing 
pathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
Run1Out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft = c(8, 8), truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
Run1Out
print("Run1 Complete!")
################RUN 1 end ########################################################################
###########################
################RUN 2 ###################################################################################################
# File parsing 
pathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
Run2Out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft = c(8, 8), truncLen=c(240,220), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
Run2Out
print("Run2 Complete!")
################RUN 2 end ########################################################################
###########################
################RUN 3 ###################################################################################################
# File parsing 
pathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/REV" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
#Replace trimLeft with Primer Lengths, truncLen based on quality plots 
Run3Out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft = c(8, 8), truncLen=c(240,210), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
Run3Out
print("Run3 Complete!")
################RUN 3 end ########################################################################
###########################
FilterAndTrimRunStats<-rbind(Run1Out, Run2Out, Run3Out)
save(FilterAndTrimRunStats, file="/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Results/dada2results/cav_out_filter_trim_dada2.Rdta")

