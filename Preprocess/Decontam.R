#########################################decontam processing#################################################################
#####################################################CAVITIES####################################################################################
#load chosen phyloseq object  & rename
#load(file.path(data.out, 'phyloseq', "phyasv.Rdata"))
ps<-phy.asv
#remove samples with 0 reads from consideration
ps<-subset_samples(ps, sample_sums(ps)>0)

#######library sizes
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
libsize_g<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType)) + geom_point()

#Frequency based method 
contamdf.freq <- isContaminant(ps, batch="Batch", method="frequency", conc="DNAQuant1000")
set.seed(100)
contam.freq_g<-plot_frequency(ps, taxa_names(ps)[which(contamdf.freq$contaminant)], conc="DNAQuant1000") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

#Prevalence based method
sample_data(ps)$is.neg <- sample_data(ps)$SampleType %in% c("Extraction", "Water", "Empty")
#table(sample_data(ps)$is.neg, sample_data(ps)$SampleType)
contamdf.prev <- isContaminant(ps, batch="Batch", method="prevalence", neg="is.neg")

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleType %in% c("Extraction", "Water", "Empty"), ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleType %in% c("Mock", "True Sample"), ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
prev_g<-ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#Run selected method
contamdf<-isContaminant(ps, batch="Batch", method=decontam.method, neg="is.neg", conc="DNAQuant1000", threshold=decontam.thresh)

#subsetting out taxa 
contaminants<-subset_taxa(ps, rownames(ps@tax_table) %in% rownames(contamdf[contamdf$contaminant==TRUE, ]))
notcontaminants<-subset_taxa(ps, !(rownames(ps@tax_table) %in% rownames(contamdf[contamdf$contaminant==TRUE, ])))

#save new phyloseq object
phy.asv_decontam<-notcontaminants
save(phy.asv_decontam, file=file.path(data.out, 'phyloseq', "phyasv_decontam.Rdata"))

#Save information from this run 
dir.create(file.path(output.out, 'Preprocess'))
save(list=c('contaminants', 'notcontaminants', 'libsize_g', 'prev_g'), file=file.path(output.out, 'Preprocess', paste0('decontam_', decontam.method, '_', paste0(decontam.thresh, collapse = '&'), '.Rdata')), compress=T)

#Print information 
print(paste0('Decontam run with the method ', decontam.method, ' and the thresholds ', paste0(decontam.thresh, collapse=' & '), ' removed ', ntaxa(contaminants), ' taxa (contaminants), of ', ntaxa(ps), ' original taxa'))

#Examples of other methods 
#combining both methods to identify contaminants
#contamdf.comb<-isContaminant(ps, batch="Batch", method="combined", neg="is.neg", conc="DNAQuant1000", threshold=0.1)
#contamdf.both<-isContaminant(ps, batch="Batch", method="both", neg="is.neg", conc="DNAQuant1000", threshold=c(0.1, 0.5))
#contamdf.either<-isContaminant(ps, batch="Batch", method="either", neg="is.neg", conc="DNAQuant1000", threshold=c(0.1, 0.3))
