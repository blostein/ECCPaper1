################################################
################################################SETUP
library(here)
library(phyloseq)
library(stringr)
library(plyr)
library(dplyr)
library(decontam)
library(stringr)
library(ggplot2)
#set paths
dada.data<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/DADA2Pipe/PRJEB35824/Results/"
#meta data directory 
#meta.data<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData"
#prj directory
prjdir<-gsub('/Code', '', here::here())
data.out<-file.path(prjdir, 'Data', 'PRJ')
code.dir<-file.path(prjdir, 'Code')
output.out<-file.path(prjdir, 'Output', 'PRJ')
#set dataset names
taxname<-"PRJ_tax_homd.rds"
seqname<-"PRJ_seqtab_homd.rds"
ampdata=F
#ampdata<-file.path(meta.data, "OriginalMetaData/SeqToCOHRAKey.Rdata")
#set variables
decontam=F #can't do decontam on these samples 
decontam.method<-'either' #see ?isContaminant for method options
decontam.thresh<-c(0.1, 0.3) #see ?isContaminant for threshold
readlimit<-1000
abund.filter<-0.05 #what abundance filter? [0 = 0% abundance to 1 = 100% abundance]
#abund.filter is used for filtering based on per sample abundance, not total sample abundance
prev<-0.05# what prevalence filter? [0=ASV present in 0% of samples to 1=ASV present in 100% of samples]
clean_meta='no_meta'#options: 'from_file', 'no_meta'
################################################################################################
###########################Get Amplification status
#not run as file output saved 
#source(file.path(code.dir, 'Preprocess', 'AmplificationStatusMeta.R'))
################################################################################################
###########################Create Phyloseq Object 
source(file.path(code.dir, 'Preprocess', 'CreatePhyloSeq.R'))
###########################Run Decontam
if(decontam==T){
  source(file.path(code.dir, 'Preprocess', 'Decontam.R'))
}
##########################Filter Samples
lessthan1000reads<-subset_samples(phy.asv, sample_sums(phy.asv)<readlimit)
prj_samplefilter<-subset_samples(phy.asv, sample_sums(phy.asv)>=readlimit)
save(prj_samplefilter, file=file.path(data.out, "phyloseq", "phyasv_postsamplefilter.Rdata"))
ps_truesamples=prj_samplefilter
#########################Filter Taxa
source(file.path(code.dir, 'Preprocess', 'FilterTaxa.R'))
prj.filter=ps.filter
save(prj.filter, file=file.path(data.out, 'phyloseq', "phyasv_postTax.Rdata"))
########################Clean Metadata
#no metadata 