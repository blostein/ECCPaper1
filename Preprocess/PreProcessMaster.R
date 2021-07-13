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
################################################
################################################Directory based coding 
#prj directory
prjdir<-gsub('/Code', '', here::here())
#original file folders - COHRA
og.data<-file.path(prjdir, 'OGData', 'COHRA')
meta.data<-file.path(og.data, 'Metadata')
dada.data<-file.path(og.data, 'DADA2')
#output and code folders
code.dir<-file.path(prjdir, 'Code')
data.out<-file.path(prjdir, 'Data')
output.out<-file.path(prjdir, 'Output')
################################################
################################################Global options - both COHRA and PRJ
readlimit<-1000
abund.filter<-0.05 #what abundance filter? [0 = 0% abundance to 1 = 100% abundance]
#abund.filter is used for filtering based on per sample abundance, not total sample abundance
prev<-0.05# what prevalence filter? [0=ASV present in 0% of samples to 1=ASV present in 100% of samples]
################################################
################################################COHRA options - decontam, filenames
#set dataset names
taxname<-"tax_final_homd_mocks.rds"
seqname<-"seqtab_final_homd.rds"
#have amp data yes/no - if yes load data
amp_data=T
if(amp_data==T){load(file.path(meta.data, "SeqToCOHRAKey.Rdata"))}

#set variables
decontam.method<-'either' #see ?isContaminant for method options
decontam.thresh<-c(0.1, 0.3) #see ?isContaminant for threshold

################################################################################################
#Run COHRA
################################################################################################
###########################Create Phyloseq Object 
source(file.path(code.dir, 'Preprocess', 'CreatePhyloSeq.R'))
###########################Run Decontam
source(file.path(code.dir, 'Preprocess', 'Decontam.R'))
##########################Filter Samples
source(file.path(code.dir, 'Preprocess', 'FilterSamples.R'))
#########################Filter Taxa
source(file.path(code.dir, 'Preprocess', 'FilterTaxa.R'))
########################Clean Metdata
################This step is highly dependent on the dataset of interest. Please use this place to load your metadata or run
#a custom metadata cleaning script. 
source(file.path(code.dir, 'Preprocess', 'CreateMetadata.R'))
########################Merge Metadata to Phloseq
meta<-as.data.frame(as.matrix(ps.filter@sam_data))
MetaVisit_noDups<-MetaVisit%>%filter(count!= 'Double; drop')
meta2<-left_join(meta, MetaVisit_noDups)
row.names(meta2)<-meta2$FoxCavID
ps_meta_nodups<-ps.filter
ps_meta_nodups@sam_data<-sample_data(meta2)
save(ps_meta_nodups, file=file.path(data.out, 'phyloseq', 'phyasv_meta.Rdata'))

#save data on read counts, sample loss and parameters 
list_ps=list(phy.asv, phy.asv_decontam, lessthan1000reads, morethan1000reads, ps_truesamples, ps.filter)
taxa_count=lapply(list_ps, function(x) ntaxa(x)); names(taxa_count)=c('all', 'decontam', 'lessthan1000', 'morethan1000', 'truesamples', 'taxafilter')
sample_count=lapply(list_ps, function(x) table(x@sam_data$SampleType))
names(sample_count)=c('all', 'decontam', 'lessthan1000', 'morethan1000', 'truesamples', 'taxafilter')
save(list=c('readlimit', 'abund.filter', 'prev', 'decontam.thresh', 'taxa_count', 'sample_count'), 
     file=file.path(output.out, 'Preprocess', paste0('preprocessParamsStats_', readlimit, '.Rdata')), compress=T)


################################################################################################
#Run PRJ
################################################################################################
output.out<-file.path(prjdir, 'Output', 'PRJ')
data.out<-file.path(prjdir, 'Data', 'PRJ')
dada.data<-file.path(prjdir, 'OGData', 'PRJ', 'DADA2')
#set dataset names
taxname<-"PRJ_tax_homd.rds"
seqname<-"PRJ_seqtab_homd.rds"
#no amplification data for PRJ
amp_data=F
###########################Create Phyloseq Object 
source(file.path(code.dir, 'Preprocess', 'CreatePhyloSeq.R'))
###########################No decontam step
##########################Filter Samples
lessthan1000reads<-subset_samples(phy.asv, sample_sums(phy.asv)<readlimit)
prj_samplefilter<-subset_samples(phy.asv, sample_sums(phy.asv)>=readlimit)
save(prj_samplefilter, file=file.path(data.out, "phyloseq", "phyasv_postsamplefilter.Rdata"))
ps_truesamples=prj_samplefilter
#########################Filter Taxa
source(file.path(code.dir, 'Preprocess', 'FilterTaxa.R'))
prj.filter=ps.filter
save(prj.filter, file=file.path(data.out, 'phyloseq', "phyasv_postTax.Rdata"))
