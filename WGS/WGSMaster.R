#libraries
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(splitstackshape)
library(dplyr)
library(tidyr)
library(tibble)
#set cutoff for deseq read limit (filter taxa/KEGG orthlogs with count<read_limit before running deseq)
read_limit=500
#prj directory
prjdir<-gsub('/Code', '', here::here())
code.dir<-file.path(prjdir, 'Code')
data.out<-file.path(prjdir, 'Data')
output.out<-file.path(prjdir, 'Output')
#WGS files in directory
wgs.path=file.path(prjdir, 'OGData', 'COHRA', 'WGS')
humann3.path=file.path(wgs.path, 'Humann3Files')
#WGS file path out directory
humann3.out= file.path(data.out, 'H3')
dir.create(file.path(data.out, 'MS'))
p=dir.create(file.path(data.out, 'MS', 'plaque'))
s=dir.create(file.path(data.out, 'MS', 'saliva'))

#################################################################clean humann3 data
#load metagenomics metadata
load(file.path(prjdir, 'OGData', 'COHRA', 'Metadata', 'metag_metashort.Rdata'))
source(file.path(code.dir, 'WGS', 'cleanHumannData.R'))

################################################################make slim sqm objects
source(file.path(code.dir, 'WGS', 'slimSQMObjects.R'))

###############################################################Run Deseq2 on sqm 
source(file.path(code.dir, 'WGS', 'runWGSDeSeq.R'))
