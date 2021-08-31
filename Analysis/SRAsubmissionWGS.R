##############################Create a SRA table for VIPS data
library(dplyr)
library(gdata)
library(purrr)
library(tidyr)
#load in dbGAP ids
dbGAPS<-read.xls("/Volumes/EPID/FACULTY/FOXMAN/COHRA II/dbGAPIDs/dbGaP_IDs_COHRA2.xlsx")

#load in VIPS meta data
my_dir="/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/ECC/Data/"
if(dir.exist(file.path(my_dir, 'SRA'))==FALSE){dir.create(my_dir, 'SRA')}
my_dir_OG="/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/ECC/OGData/COHRA/Metadata/"

load(file.path(my_dir_OG, 'metag_metashort.Rdata'))
load(file.path(my_dir, 'Metadata', 'MetaVisit.Rdata'))

MV=MetaVisit%>%
  mutate(host_age=ifelse(AgeAtExamMonths<36, '<3 years', '3-5 years'),
         host_sex=BabySex)%>%
  dplyr::select(COHRAID, FoxCavID, BabySubjectID, Site, VisitDate, host_sex, host_age)%>%unique()

meta=left_join(metag_metashort%>%select(-Site)%>%
               mutate(BabySubjectID=as.integer(BabySubjectID)),
               MV)
#include the following variables in the SRA table
#-SRAID
#-dbGAPID
#Sex
#collection date

dbGAPS$BabySubjectID<-dbGAPS$COHRA2_ID
dbGAPS<-dbGAPS %>% filter(BabySubjectID %in% meta$BabySubjectID) %>% unique()
SRA<-left_join(meta, dbGAPS)
dim(SRA %>% unique())
SRA$sample_name=SRA$SampleName
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
         env_medium = SampleType, 
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



SRA3<-SRA2 %>% select(-FoxCavID, -Site, -VisitDate, -SampleType) %>% select(-host_age, -host_sex, -dbGaP_ID, everything())
head(SRA3)
SRA3<-SRA3%>%arrange(sample_name)
write.csv(SRA3, file=file.path(my_dir, 'SRA', "ECC_WGS_biosample.csv"))
write.table(SRA3, file=file.path(my_dir, 'SRA', "ECC_WGS_biosample.tsv"), row.names=FALSE, sep="\t")

## make SRA metadata
files_wgs=read.table(file.path(my_dir, 'SRA', 'wgs_samples.txt'))
files2=files_wgs %>% pivot_wider(names_from=V3, values_from=V2)
colnames(files2)=c('SampleName', 'read1', 'read2')
files2<-left_join(files2, meta)%>%mutate(sample_name=SampleName)
nrow(files2)==nrow(SRA3)
files2$library_ID<-files2$OriginalID
files2<-files2 %>% select(sample_name, library_ID, read1, read2)%>%
  mutate(read_unpaired=gsub('R2_paired_decon.fq.gz', 'unpaired_decon.fq.gz', read2))

files3=files2 %>% 
  mutate(sample_name = sample_name, #THIS IS THE SAMPLEID FROM YOUR METADATA
         library_ID = library_ID,
         title= 'Whole genome shotgun sequencing',
         library_strategy='WGS', 
         library_source='METAGENOMIC',
         library_selection='RANDOM PCR', 
         library_layout='paired',
         platform='ILLUMINA',
         instrument_model='Illumina HiSeq 4000', #see email from COSMOS ID 
         design_description= 'whole genome shotgun sequencing',
         filetype='fastq',
         filename=read1,
         filename2=read2, 
         filename3=read_unpaired)

write.csv(files3, file=file.path(my_dir, 'SRA', "ECC1_WGS_SRA.csv"))
write.table(files3, file=file.path(my_dir, 'SRA', "ECC1_WGS_SRA.tsv"), row.names=FALSE, sep="\t")


#filenames for submission          
filenames<-c(files2$read1, files2$read2)
filenames<-paste("/scratch/bfoxman_root/bfoxman/blostein/meta_CAVITIES/Files/CleanPipe/PostHR/", filenames, sep='')
write.table(filenames, file.path(my_dir, 'SRA',"ECCWGS.txt"), quote = F, row.names = F)
