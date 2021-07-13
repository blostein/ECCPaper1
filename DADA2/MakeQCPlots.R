library(dada2); packageVersion("dada2")
library(ggplot2)
####
resultpath<- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Results/dada2results/qcplots"
################RUN 1 ###################################################################################################
pathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/FWD/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/REV/" # CHANGE ME ...
fastqFs <- paste(pathF, sort(list.files(pathF, pattern="fastq.gz")), sep="")
fastqRs <- paste(pathR, sort(list.files(pathR, pattern="fastq.gz")), sep="")
#plot forwards for run1 
max <- 9
x <- seq_along(fastqFs)
FwdList<-split(fastqFs, ceiling(x/max))
FWDplots <- lapply(FwdList, plotQualityProfile)
names(FWDplots)=seq(1:length(FwdList))
lapply(names(FWDplots), 
       function(x) ggsave(filename=paste(resultpath, "Run1_F", x,".pdf",sep=""), plot=FWDplots[[x]]))
#plot reverse for run1 (samples 1-48 in PTB)
x <- seq_along(fastqRs)
RevList<-split(fastqRs, ceiling(x/max))
Revplots <- lapply(RevList, plotQualityProfile)
names(Revplots)=seq(1:length(RevList))
lapply(names(Revplots), 
       function(x) ggsave(filename=paste(resultpath, "Run1_R", x,sep=""), device="pdf", plot=Revplots[[x]]))
################RUN 1 end ###################################################################################################


################RUN 2 ###################################################################################################
pathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/FWD/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/REV/" # CHANGE ME ...
fastqFs <- paste(pathF, sort(list.files(pathF, pattern="fastq.gz")), sep="")
fastqRs <- paste(pathR, sort(list.files(pathR, pattern="fastq.gz")), sep="")
#plot forwards for run2 (samples 49-96 in PTB)
x <- seq_along(fastqFs)
FwdList<-split(fastqFs, ceiling(x/max))
FWDplots <- lapply(FwdList, plotQualityProfile)
names(FWDplots)=seq(1:length(FwdList))
lapply(names(FWDplots), 
       function(x) ggsave(filename=paste(resultpath, "Run2_F", x,sep=""), device="pdf", plot=FWDplots[[x]]))
#plot reverse for run1 (samples 49-96 in PTB)
x <- seq_along(fastqRs)
RevList<-split(fastqRs, ceiling(x/max))
Revplots <- lapply(RevList, plotQualityProfile)
names(Revplots)=seq(1:length(RevList))
lapply(names(Revplots), 
       function(x) ggsave(filename=paste(resultpath, "Run2_R", x,sep=""), device="pdf", plot=Revplots[[x]]))
################RUN 2 end ###################################################################################################

################RUN 3 ###################################################################################################
pathF <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/FWD/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/REV/" # CHANGE ME ...
fastqFs <- paste(pathF, sort(list.files(pathF, pattern="fastq.gz")), sep="")
fastqRs <- paste(pathR, sort(list.files(pathR, pattern="fastq.gz")), sep="")
#plot forwards for run3
x <- seq_along(fastqFs)
FwdList<-split(fastqFs, ceiling(x/max))
FWDplots <- lapply(FwdList, plotQualityProfile)
names(FWDplots)=seq(1:length(FwdList))
lapply(names(FWDplots), 
       function(x) ggsave(filename=paste(resultpath, "Run3_F", x,sep=""), device="pdf", plot=FWDplots[[x]]))
#plot reverse for run3
x <- seq_along(fastqRs)
RevList<-split(fastqRs, ceiling(x/max))
Revplots <- lapply(RevList, plotQualityProfile)
names(Revplots)=seq(1:length(RevList))
lapply(names(Revplots), 
       function(x) ggsave(filename=paste(resultpath, "Run3_R", x,sep=""), device="pdf", plot=Revplots[[x]]))
################RUN 3 end ###################################################################################################


