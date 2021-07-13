library(tidyr)
library(dplyr)
library('stringr')
library('stringi')
library(Biostrings)
library(phyloseq)
library('msa')

meta.data<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData"
#prj directory
prjdir<-gsub('/Code', '', here::here())
data.out<-file.path(prjdir, 'Data')
load(file.path(data.out, 'phyloseq', 'phyasv_meta.Rdata'))
ps_sm=phyloseq::subset_taxa(ps_meta_nodups, Species=='mutans')

ref=readDNAStringSet('~/Downloads/smut25175.fasta.txt')
  
contigs.table=read.table('/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/GreatLakesWGSPipe/metasqueeze/plaque_assembly/intermediate/19.plaque_assembly.contigsinbins')
tx=readLines('/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/GreatLakesWGSPipe/metasqueeze/plaque_assembly/results/02.plaque_assembly.rnas')
rRNAs=readDNAStringSet('/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/GreatLakesWGSPipe/metasqueeze/plaque_assembly/results/02.plaque_assembly.rnas')

seq_name=names(rRNAs)
seq_name=as.data.frame(do.call(rbind, strsplit(seq_name, '\t|;|\\|')), fill=T)
names(rRNAs)=sub("^(megahit[^_]*_[^\t]*).*", "\\1", names(rRNAs))

smut=contigs.table %>% filter(V3 %in% c('maxbin.001.fasta.contigs', 'maxbin.003.fasta.contigs'))
smut_contigs=smut$V1
smut_16S=seq_name %>% filter(V3 %in% smut_contigs & str_detect(V2, '16S'))

smut16s=rRNAs[names(rRNAs)%in% smut_16S$V1]

t=c(smut16s, ref)
t2=c(ps_sm@refseq, ref)
t3=c(smut16s[1], ps_sm@refseq, ref)
mult=msa(t, method="ClustalW", type="dna", order="input")
mult2=msa(t2, method="ClustalW", type="dna", order="input")
mult3=msa(t3, method="Muscle", type="dna", order="input")
