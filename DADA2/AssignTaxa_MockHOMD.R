library(dada2); packageVersion("dada2")
# load seqtab
seqtab <- readRDS("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/seqtab_final_homd.rds")
# Assign taxonomy to genus, then to species
taxa <- assignTaxonomy(seqtab, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/ref_tax/HOMD_assigntaxa_togenus_plusmocks.fasta.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/ref_tax/HOMD_AssignSpecies_plusmocks.fasta.gz")
#assign taxonomy all the way to species
taxa_toSp<-assignTaxonomy(seqtab, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/ref_tax/HOMD_AssignTaxaToSpecies_plusmocks.fasta.gz", multithread=TRUE)
# Write to disk
#saveRDS(seqtab, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/seqtab_final_homd.rds") # CHANGE ME to where you want sequence table saved

saveRDS(taxa, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/tax_final_homd_mocks.rds") # CHANGE ME ...
saveRDS(taxa_toSp, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/tax_toSp_final_homd_mocks.rds") # CHANGE ME ...

