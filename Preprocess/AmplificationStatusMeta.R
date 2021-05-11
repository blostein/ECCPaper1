library(stringr)
metadata<-read.csv("/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData/OriginalMetaData/SeqToCOHRAKey.csv")
#AmplificationStatus is based on the color code in the platemaps from the code. change blanks to pass
#metadata$AmplificationStatus<-ifelse(metadata$AmplificationStatus!="Fail", "Pass", as.character(metadata$AmplificationStatus))
#read in found datastatus 
#amplificationstatus<-read.csv("/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData/OriginalMetaData/AmplificationKey.csv")
#join together
#metadata<-left_join(metadata, amplificationstatus)
#check status (an old variable describing amplification + found status, see /Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/Cavities/Code/Markdown_phylseq.RMD) vs 
with(metadata, xtabs(~Status+AmplificationStatus)) #matches previous description 
#extract IDS in dada2 data, not in metadata 
#rownames(seqtab)[!(rownames(seqtab) %in% metadata$FoxCavID)] #none, pass check
#create metadata variable for true sample or control 
metadata$SampleType<-case_when(str_detect(metadata$COHRAID, "Mock")~"Mock", 
                               str_detect(metadata$COHRAID, "Extraction")~"Extraction",
                               str_detect(metadata$COHRAID, "Negative")~"Water",
                               metadata$FoundStatus=="NotFound" ~ "Empty", 
                               metadata$FoundStatus=="Found" ~ "True Sample")
#create metadata variable for batch 
metadata$Batch = case_when(metadata$Plate %in% c(1, 2, 3, 4)==TRUE ~ 1, 
                           metadata$Plate %in% c(5, 6, 7, 8)==TRUE ~ 2, 
                           metadata$Plate %in% c(9, 10, 11, 12)==TRUE ~3)
with(metadata, xtabs(~Batch+Plate))
#rename plate variable
metadata$Plate<-factor(paste("Plate ", metadata$Plate, sep=''), levels=c("Plate 1", "Plate 2", "Plate 3", "Plate 4", "Plate 5", "Plate 6", "Plate 7", "Plate 8", "Plate 9", "Plate 10", "Plate 11", "Plate 12"))

#Multiply quant by 1000 and bottom at 0 
metadata$DNAQuant1000<-ifelse(metadata$DNAQuant<=0.01, 0.01, metadata$DNAQuant)*100
rownames(metadata)<-metadata$FoxCavID
save(metadata, file="/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData/OriginalMetaData/SeqToCOHRAKey.Rdata")
