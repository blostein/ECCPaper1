
#Code
tax<-readRDS(file.path(dada.data, taxname))
seqtab<-readRDS(file.path(dada.data, seqname))

if(amp_data==T){
phy.asv<- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                   tax_table(tax), sample_data(metadata))
}
if(amp_data==F){
  phy.asv<- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                     tax_table(tax), sample_data(data.frame('SampleID'=row.names(seqtab), row.names = row.names(seqtab))))
}
dna <- Biostrings::DNAStringSet(taxa_names(phy.asv))
names(dna) <- taxa_names(phy.asv)
phy.asv <- merge_phyloseq(phy.asv, dna)
taxa_names(phy.asv) <- paste0("ASV", seq(ntaxa(phy.asv)))
#add plate data
phy.asv@sam_data$Run=ifelse(phy.asv@sam_data$Plate %in% paste('Plate', c(1, 2, 3, 4)), 'Run 1',
                       ifelse(phy.asv@sam_data$Plate %in% paste('Plate', c(5, 6, 7, 8)), 'Run 2', 
                        ifelse(phy.asv@sam_data$Plate%in% paste('Plate', c(9, 10, 11, 12)), 'Run 3',
                                     NA)))

dir.create(file.path(data.out, 'phyloseq'))
save(phy.asv, file=file.path(data.out, 'phyloseq', "phyasv.Rdata"))