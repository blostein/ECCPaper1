##############################Create a SRA table for VIPS data
library(dplyr)
library(gdata)
library(purrr)
#load in dbGAP ids
dbGAPS<-read.xls("/Volumes/EPID/FACULTY/FOXMAN/COHRA II/dbGAPIDs/dbGaP_IDs_COHRA2.xlsx")

#load in VIPS meta data
my_dir="/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/ECC/Data/"
dir.create(my_dir, 'SRA')
load(file.path(my_dir, "phyloseq", "phyasv_truesamples.Rdata"))
load(file.path(my_dir, 'Metadata', 'MetaVisit.Rdata'))
meta<-as.data.frame(as.matrix(ps_truesamples@sam_data)) %>% select(COHRAID, FoxCavID, Plate, Batch, SampleType)%>%filter(SampleType!='Empty')
MV=MetaVisit%>%
  mutate(host_age=ifelse(AgeAtExamMonths<36, '<3 years', '3-5 years'),
         host_sex=BabySex)%>%
  dplyr::select(COHRAID, FoxCavID, BabySubjectID, Site, VisitDate, host_sex, host_age)%>%unique()

meta=left_join(meta, MV)
#include the following variables in the SRA table
#-SRAID
#-dbGAPID
#Sex
#collection date

dbGAPS$BabySubjectID<-dbGAPS$COHRA2_ID
dbGAPS<-dbGAPS %>% filter(BabySubjectID %in% meta$BabySubjectID) %>% unique()
SRA<-left_join(meta, dbGAPS)
dim(SRA %>% unique())
SRA$sample_name=SRA$FoxCavID
SRA$Site=paste0('USA: ', SRA$Site)

#Creat the biosample table table: 
SRA2=SRA %>% ungroup() %>% select(-BabySubjectID, -COHRA2_ID, -COHRAID) %>% 
  mutate(sample_name = sample_name, #THIS IS THE SAMPLEID FROM YOUR METADATA
         sample_title = "",
         bioproject_accession = "",
         organism = "human metagenome", 
         strain = "not applicable", 
         isolate = "not applicable", 
         cultivar = "not applicable", 
         ecotype = "not applicable", 
         collection_date = VisitDate, 
         env_broad_scale = "not applicable", 
         env_local_scale = "not applicable" , 
         env_medium = "saliva", 
         geo_loc_name = Site, 
         host = "Homo sapiens", 
         isol_growth_condt = "not applicable", 
         lat_lon = "not applicable", 
         ethnicity='Caucasian',
         host_age=host_age,
         host_sex=host_sex)

#after this everything is metadata, so can add whatever columns you want
#I already have relationship, visit, and site in here so I won't add those
#add in all other relavent metadata 



SRA3<-SRA2 %>% select(-FoxCavID, -Site, -VisitDate) %>% select(-host_age, -host_sex, -dbGaP_ID, -Plate, -Batch, -SampleType, everything())
head(SRA3)
SRA3<-SRA3%>%arrange(sample_name)
write.csv(SRA3, file=file.path(my_dir, 'SRA', "ECC1_biosample.csv"))
write.table(SRA3, file=file.path(my_dir, 'SRA', "ECC1_biosample.tsv"), row.names=FALSE, sep="\t")

## make SRA metadata
path_16s="/Volumes/EPID/FACULTY/FOXMAN/COHRA II/Cavities Project 2018/16S Data"
path1 <- file.path(path_16s, "MiSeq M04695 2018 Run 150 - Foxman FoxCav 1-4/")
path2 <- file.path(path_16s, "MiSeq M02127 2018 R337 - Foxman FoxCav 5-8/")
path3 <- file.path(path_16s, 'MiSeq for FoxCav 9 and 10/')
paths=list(path1, path2, path3)
Fs=map_dfr(paths, ~(data.frame('read1'=list.files(.x, pattern="R1_001.fastq.gz", full.names = F))))
Rs=map_dfr(paths, ~(data.frame('read2'=list.files(.x, pattern="R2_001.fastq.gz", full.names = F))))

files<-cbind(Fs, Rs)
files$FoxCavID=gsub('_.*', '', files$read1)

files2<-left_join(files, meta)%>%mutate(sample_name=FoxCavID) %>%filter(SampleType!='Empty')
nrow(files2)==nrow(SRA3)
files2$library_ID<-files2$FoxCavID
files2<-files2 %>% select(sample_name, library_ID, read1, read2)

files3=files2 %>% 
  mutate(sample_name = sample_name, #THIS IS THE SAMPLEID FROM YOUR METADATA
         library_ID = library_ID,
         title= 'V4 of 16S rRNA amplicon',
         library_strategy='AMPLICON', 
         library_source='METAGENOMIC',
         library_selection='PCR',
         library_layout='paired',
         platform='ILLUMINA',
         instrument_model='Illumina MiSeq', 
         design_description= 'Amplified V4 region of 16S rRNA amplicon',
         filetype='fastq',
         filename=read1,
         filename2=read2)

write.csv(files3, file=file.path(my_dir, 'SRA', "ECC1_SRA.csv"))
write.table(files3, file=file.path(my_dir, 'SRA', "ECC1_SRA.tsv"), row.names=FALSE, sep="\t")


#filenames for submission          
filenames<-c(files2$read1, files2$read2)
filenames<-paste("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/OGFiles/*/", filenames, sep='')
write.table(filenames, file.path(my_dir, 'SRA',"ECC16SFilenames.txt"), quote = F, row.names = F)


### also include controls (have to do separately as did not realize)
load(file.path(my_dir, "phyloseq", "phyasv.Rdata"))
meta<-as.data.frame(as.matrix(phy.asv@sam_data)) %>% select(COHRAID, FoxCavID, Plate, Batch, SampleType)%>%filter(!(SampleType%in%c('Empty', 'True Sample')))

SRA2=meta %>% ungroup() %>% 
  mutate(sample_name = paste0('FoxCav_', COHRAID), #THIS IS THE SAMPLEID FROM YOUR METADATA
         sample_title = "",
         bioproject_accession = "",
         organism = "human metagenome", 
         strain = "not applicable", 
         isolate = "not applicable", 
         cultivar = "not applicable", 
         ecotype = "not applicable", 
         collection_date = "not applicable", 
         env_broad_scale = "not applicable", 
         env_local_scale = "not applicable" , 
         env_medium = "saliva", 
         geo_loc_name = 'not applicable', 
         host = "Homo sapiens", 
         isol_growth_condt = "not applicable", 
         lat_lon = "not applicable", 
         ethnicity='not applicable')

#after this everything is metadata, so can add whatever columns you want
#I already have relationship, visit, and site in here so I won't add those
#add in all other relavent metadata 



SRA3<-SRA2 %>% select(-FoxCavID, -COHRAID) %>% select(-Plate, -Batch, -SampleType, everything())
head(SRA3)
SRA3<-SRA3%>%arrange(sample_name)
write.csv(SRA3, file=file.path(my_dir, 'SRA', "ECC1_controls_biosample.csv"))
write.table(SRA3, file=file.path(my_dir, 'SRA', "ECC1_controls_biosample.tsv"), row.names=FALSE, sep="\t")

#sra file 
#need to get the right merge, which is difficult because I renamed the foxcav number to run3CONTROL for run 3 since exact same names
#were used for run 2 and run 3
#extract which sample names were the run 3 ones
Flist=map(paths, ~(data.frame('read1'=list.files(.x, pattern="R1_001.fastq.gz", full.names = F))))
my_problems=c(Flist[[3]]$read1[str_detect(Flist[[3]]$read1, 'FoxCav', negate=T)])
files2=files %>% mutate(FoxCavID=ifelse(read1 %in% my_problems, paste0('run3', FoxCavID), FoxCavID))

files2<-left_join(meta, files2, by="FoxCavID")%>%mutate(sample_name=COHRAID)
nrow(files2)==nrow(SRA3)
files2$library_ID<-files2$FoxCavID
files2<-files2 %>% select(sample_name, library_ID, read1, read2)

files3=files2 %>% 
  mutate(sample_name =paste0('FoxCav_', sample_name), #THIS IS THE SAMPLEID FROM YOUR METADATA
         library_ID = paste0('FoxCav_controls_', seq(1:nrow(files2))),
         title= 'V4 of 16S rRNA amplicon',
         library_strategy='AMPLICON', 
         library_source='METAGENOMIC',
         library_selection='PCR',
         library_layout='paired',
         platform='ILLUMINA',
         instrument_model='Illumina MiSeq', 
         design_description= 'Amplified V4 region of 16S rRNA amplicon',
         filetype='fastq',
         filename=read1,
         filename2=read2)

write.csv(files3, file=file.path(my_dir, 'SRA', "ECC1_controls_SRA.csv"))
write.table(files3, file=file.path(my_dir, 'SRA', "ECC1_controls_SRA.tsv"), row.names=FALSE, sep="\t")


#filenames for submission          
filenames<-c(files2$read1, files2$read2)
filenames<-paste("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/OGFiles/*/", filenames, sep='')
write.table(filenames, file.path(my_dir, 'SRA',"ECC16SFilenames_controls.txt"), quote = F, row.names = F)
