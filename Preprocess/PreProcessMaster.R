################################################
################################################SETUP
#set path
dada.data<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/DADA2Pipe/Data/"
prjdir<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/Paper1"
data.out<-file.path(prjdir, 'Data')
code.dir<-file.path(prjdir, 'Code')
#set dataset names
taxname<-"tax_final_homd_mocks.rds"
seqname<-"seqtab_final_homd.rds"
ampdata<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData/OriginalMetaData/SeqToCOHRAKey.Rdata"
#set variables
readlimit<-1000
abund.filter<-0.05 #what abundance filter? [0 = 0% abundance to 1 = 100% abundance]
#abund.filter is used for filtering based on per sample abundance, not total sample abundance
prev<-0.05# what prevalence filter? [0=ASV present in 0% of samples to 1=ASV present in 100% of samples]
################################################################################################
###########################Get Amplification status
#not run as file output saved 
#source(file.path(code.dir, 'Preprocess', 'AmplificationStatusMeta.R'))
################################################################################################
###########################Create Phyloseq Object 
source(file.path(code.dir, 'Preprocess', 'CreatePhyloSeq.R'))
###########################Run Decontam



##########################Filter Samples



#########################Filter Taxa


########################Clean Metdata
################This step is highly dependent on the dataset of interest. Please use this place to load your metadata or run
#a custom metadata cleaning script 

########################Merge Metadata to Phloseq