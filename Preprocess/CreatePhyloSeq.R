# packages
library(phyloseq)
library(stringr)
library(plyr)
library(dplyr)
library(stringr)
#Code
tax<-readRDS(file.path(dada.data, taxname))
seqtab<-readRDS(file.path(dada.data, seqname))
amp_meta<-load(file.path(ampdata))

phy.asv<- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                   tax_table(tax), sample_data(metadata))
dna <- Biostrings::DNAStringSet(taxa_names(phy.asv))
names(dna) <- taxa_names(phy.asv)
phy.asv <- merge_phyloseq(phy.asv, dna)
taxa_names(phy.asv) <- paste0("ASV", seq(ntaxa(phy.asv)))
dir.create(file.path(data.out, 'phyloseq'))
save(phy.asv, file=file.path(data.out, 'phyloseq', "phyasv.Rdata"))