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
readcountsByType<-ggplot(data.frame(readcounts=sample_sums(ps), SampleType=ps@sam_data$SampleType)%>%filter(SampleType!='Empty'), aes(fill=SampleType, x=readcounts))+geom_histogram()+facet_wrap(~SampleType, ncol=1)+geom_vline(xintercept=readlimit, color='red', linetype=2)+theme_classic()
readcountsbyType_preDecontam<-ggplot(data.frame(readcounts=sample_sums(phy.asv), SampleType=phy.asv@sam_data$SampleType)%>%filter(SampleType!='Empty'), aes(fill=SampleType, x=readcounts))+geom_histogram()+facet_wrap(~SampleType, ncol=1)+geom_vline(xintercept=readlimit, color='red', linetype=2)+theme_classic()

  
weirdFail<-subset_samples(ps, AmplificationStatus=="Fail" & sample_sums(ps)>=readlimit)

#less than 1000 reads 
lessthan1000reads<-subset_samples(ps, sample_sums(ps)<readlimit)
morethan1000reads<-subset_samples(ps, sample_sums(ps)>=readlimit)
#with(lessthan1000reads@sam_data, xtabs(~SampleType+AmplificationStatus))
#with(morethan1000reads@sam_data, xtabs(~SampleType+AmplificationStatus))
#passed but shouldn't have 
EmptyPass<-subset_taxa(subset_samples(ps, sample_sums(ps)>=readlimit & SampleType=="Empty"), taxa_sums(subset_samples(ps, sample_sums(ps)>=1000 & SampleType=="Empty"))>500)
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

#plot positive controls
mock_true<-read.csv(file.path(meta.data, "ZymoMockD6305.csv"))
mock_true$SciName=paste(mock_true$Genus, mock_true$Species)
mock_true$Name="ZymoBIOMIC D6306"
ps_mock =transform_sample_counts(subset_taxa(subset_samples(phy.asv, SampleType=="Mock"), taxa_sums(subset_samples(phy.asv, SampleType=="Mock"))>50), function(OTU) OTU/sum(OTU))
ps_mock_bar=plot_grid(
  plot_bar(ps_mock, 'Plate', 'Abundance', 'Genus')+facet_wrap(~Run, scales='free_x')+theme_bw()+theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 45, vjust = .4)), 
  ggplot(data=mock_true, aes(x=Name, y=X16S.Only/100, fill=SciName))+theme_bw()+geom_bar(position="stack", stat="identity", width=.4, color="black")+scale_x_discrete(name="", labels=c("True mock"))+theme(legend.text = element_text(face="italic"), legend.title=element_blank(), legend.position = "right")+scale_y_continuous(name="")+facet_wrap(~Name)+theme(axis.text.x = element_text(angle = 45, vjust = .4)),
  align='h', axis = 'tb', rel_widths = c(1, 1))
weirdmock=subset_taxa(subset_samples(phy.asv, AmplificationStatus=='Fail' & SampleType=='Mock'), taxa_sums(subset_samples(phy.asv, AmplificationStatus=='Fail' & SampleType=='Mock'))>0)
weirdmock_bar=plot_bar(weirdmock, 'Abundance', 'Genus', 'Genus')+theme_bw()+theme(legend.position = 'none', axis.text.y = element_text(face='italic'))

#plot extraction controls
Extractionkit<-subset_taxa(subset_samples(ps, SampleType=="Extraction"), taxa_sums(subset_samples(ps, SampleType=="Extraction"))>0)
#plot weird water sample
weirdwater=subset_taxa(subset_samples(morethan1000reads, SampleType=='Water'), taxa_sums(subset_samples(morethan1000reads, SampleType=='Water'))>0)
weirdwater_bar=plot_bar(weirdwater, 'Abundance', 'Genus', 'Genus')+theme_bw()+theme(legend.position = 'none', axis.text.y = element_text(face='italic'))
#save
save(list=c('lessthan1000reads', 'sampleTypeVsReadLimit', 'readcountsByAmp', 'readcountsByType', 'readcountsbyType_preDecontam', 'ps_mock_bar', 'weirdwater_bar'), 
     file=file.path(output.out, 'Preprocess', paste0('filtersamples_', readlimit, '.Rdata')), compress=T)

#print
print(paste0('Filtering out all samples with total reads of <', readlimit, ' as well as any samples with that many reads but that were not true samples, and samples that were missing all metadata, of an original ', nsamples(ps), ' samples, ', nsamples(ps_truesamples), ' are left', 
      ' See SampleTypeVsReadLimit tibble for more info'))
