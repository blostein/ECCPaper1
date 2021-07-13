library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run1/seqtab_run1.rds")
st2 <- readRDS("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run2/seqtab_run2.rds")
st3 <- readRDS("/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/Run3/seqtab_run3.rds")
rownames(st3)[rownames(st3) %in% rownames(st2)] <- paste('run3', rownames(st3)[rownames(st3) %in% rownames(st2)], sep="")
st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
dim(st.all)
dim(seqtab)
sum(st.all)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- rowSums(seqtab)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track_postchimera<-track
saveRDS(track_postchimera, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/track_postchimera_HOMD.rds")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/ref_tax/HOMD_assigntaxa_togenus.fasta.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/ref_tax/HOMD_AssignSpecies.fasta.gz")
# Write to disk
saveRDS(seqtab, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/seqtab_final_homd.rds") # CHANGE ME to where you want sequence table saved

saveRDS(taxa, "/scratch/bfoxman_root/bfoxman/blostein/CAVITIES/Files/tax_final_homd.rds") # CHANGE ME ...
