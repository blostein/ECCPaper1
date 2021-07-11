#########################################filter samples#################################################################
#####################################################CAVITIES####################################################################################
#################################################################################################################################################
#load chosen phyloseq object  & rename
#load(file.path(data.out, 'phyloseq', "phyasv_decontam.Rdata"))
ps<-phy.asv_decontam
#chose variables 
#Samples with missing metadata (entirely missing) ####CAVITIES dependent
missing_metadata<-c("11000016-7", "21000212-8", "21000438-3")

#Checks
sampVfound<-with(ps@sam_data, xtabs(~SampleType+FoundStatus))
sampVamp<-with(ps@sam_data, xtabs(~SampleType+AmplificationStatus))

#summary of read count in fail vs pass samples 
readcountsByAmp<-ggplot(data.frame(readcounts=sample_sums(ps), AmplificationStatus=ps@sam_data$AmplificationStatus), aes(fill=AmplificationStatus, x=readcounts))+geom_histogram()+facet_wrap(~AmplificationStatus, ncol=1)+geom_vline(xintercept=readlimit, color='red', linetype=2)+theme_classic()
readcountsByType<-ggplot(data.frame(readcounts=sample_sums(ps), AmplificationStatus=ps@sam_data$SampleType), aes(fill=AmplificationStatus, x=readcounts))+geom_histogram()+facet_wrap(~AmplificationStatus, ncol=1)+geom_vline(xintercept=readlimit, color='red', linetype=2)+theme_classic()
  
weirdFail<-subset_samples(ps, AmplificationStatus=="Fail" & sample_sums(ps)>=readlimit)

#less than 1000 reads 
lessthan1000reads<-subset_samples(ps, sample_sums(ps)<readlimit)
morethan1000reads<-subset_samples(ps, sample_sums(ps)>=readlimit)
#with(lessthan1000reads@sam_data, xtabs(~SampleType+AmplificationStatus))
#with(morethan1000reads@sam_data, xtabs(~SampleType+AmplificationStatus))
#passed but shouldn't have 
EmptyPass<-subset_taxa(subset_samples(ps, sample_sums(ps)>=readlimit & SampleType=="Empty"), taxa_sums(subset_samples(ps, sample_sums(ps)>=1000 & SampleType=="Empty"))>500)
Extractionkit<-subset_taxa(subset_samples(ps, SampleType=="Extraction"), taxa_sums(subset_samples(ps, SampleType=="Extraction"))>0)
#plot_bar(EmptyPass, fill="Genus");plot_bar(Extractionkit, fill='Genus')

#final filter
#exclude samples>1000 reads that are empty, mock or water 
ps_truesamples<-subset_samples(morethan1000reads, SampleType=="True Sample")
ps_truesamples<-subset_taxa(ps_truesamples, taxa_sums(ps_truesamples)>0)

#exclude samples with entirely missing metadata
ps_truesamples<-subset_samples(ps_truesamples, !(COHRAID%in%missing_metadata))
#save true samples
save(ps_truesamples, file=file.path(data.out, 'phyloseq', "phyasv_truesamples.Rdata"))

#Count readlimit fail by sample type 
ps@sam_data$ReadLimit=ifelse(sample_sums(ps)>=readlimit, 'pass', 'fail')
sampleTypeVsReadLimit=ps@sam_data %>% group_by(SampleType, ReadLimit) %>% dplyr::summarise(n=n())

#save
save(list=c('lessthan1000reads', 'sampleTypeVsReadLimit', 'readcountsByAmp', 'readcountsByType'), file=file.path(output.out, 'Preprocess', paste0('filtersamples_', readlimit, '.Rdata')), compress=T)

#print
print(paste0('Filtering out all samples with total reads of <', readlimit, ' as well as any samples with that many reads but that were not true samples, and samples that were missing all metadata, of an original ', nsamples(ps), ' samples, ', nsamples(ps_truesamples), ' are left', 
      ' See SampleTypeVsReadLimit tibble for more info'))
