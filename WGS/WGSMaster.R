#libraries
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(splitstackshape)
library(dplyr)
library(tidyr)
library(tibble)
library(SQMtools)
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

###############################################################Pull out and save
#taxonomic annotations from the top significant and oxidative phosphorylation related 
#KEGG orthologs (to save time when running Results.RMD and Supplement.RMD)
#pull out taxa from highly significant KEGGS
sig_lessthan10=r3%>%filter(-log10(padj)>10)%>%pull(KEGG)
sig_lessthan10list=list()
for(i in 1:length(sig_lessthan10)){
  print(paste('starting', i))
  df=subsetFun(plaque, fun=paste(sig_lessthan10[i]))
  df=as.data.frame(df$taxa$species$abund)%>%rownames_to_column('Species')%>%pivot_longer(starts_with('Sample'), names_to='SampleName', values_to='Abundance')
  sig_lessthan10list[[i]]=df; rm(df)
}
names(sig_lessthan10list)=sig_lessthan10
sig_taxa=lapply(sig_lessthan10list, function(x) x %>%select(Species)%>%unique())%>%bind_rows(.id='p')%>%group_by(p)%>%dplyr::summarise(Taxa=paste(Species, collapse='; '))%>%mutate(KEGG=p)

op_lists=r3%>%filter(padj<0.05 & PathAnnotation_01_3=='Oxidative phosphorylation')%>%mutate(p=word(kname, 1, 3))%>%group_by(p)%>%dplyr::summarise(k=list(KEGG))
op_taxa_list=list()
for(i in 1:nrow(op_lists)){
  print(paste('starting', i))
  df=subsetFun(plaque, fun=paste(op_lists$k[[i]], collapse='|'))
  df=as.data.frame(df$taxa$genus$abund)%>%rownames_to_column('Species')%>%pivot_longer(starts_with('Sample'), names_to='SampleName', values_to='Abundance')
  op_taxa_list[[i]]=df; rm(df)
}
names(op_taxa_list)=word(op_lists$p, 1)

save(list=c('sig_taxa', 'op_taxa_list'), file=file.path(data.out, 'MS', 'taxpathannotation.Rdata'))

