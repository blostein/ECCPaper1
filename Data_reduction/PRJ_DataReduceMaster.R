################################################
################################################SETUP
library(DirichletMultinomial)
#no meta data for prj
#meta.data<-"/Volumes/EPID/FACULTY/FOXMAN/Freida Blostein/CAVITIES_meta/MetaData"
#prj directory
prjdir<-gsub('/Code', '', here::here())
data.out<-file.path(prjdir, 'Data', 'PRJ')
code.dir<-file.path(prjdir, 'Code')
output.out<-file.path(prjdir, 'Output', 'PRJ')
dir.create(file.path(output.out, 'Network')); dir.create(file.path(data.out, 'Network'))
dir.create(file.path(output.out, 'CST')); dir.create(file.path(data.out, 'CST'))
#dir.create(file.path(output.out, 'RF')); dir.create(file.path(data.out, 'RF'))
################################################################################  
my_data='PRJ'
#set variables
transform=c('hellinger', 'clr') #one of 'hellinger', 'clr', 'compositional' or other options in ?microbiome::transform

#WCGNA specific
minModuleSize = 5;
#deep split ds
splitSize=4

#RF specific RF not run for Proj, no metadata
#duplicates=c(T, F)
#incident_visit=c(F)
#rps<-5
#num<-10

#CST specific 
k= c(4, 5, 6)
thresh.8<-0.8
thresh.9<-0.9
thresh2<-0.1

################################################################################  
##############load datasets
load(file.path(data.out, 'phyloseq', 'phyasv_postTax.Rdata'))
ps_compositional=prj.filter %>% microbiome::transform(., 'compositional') # transform to relative abundance
  
for(i in 1:length(transform)){
  #apply transform
  transform_i=transform[i]; print(transform_i)
  ps<-microbiome::transform(prj.filter, transform_i) # change to reflect right ps name
  #loop over WCGNA options
  for(j in 1:length(minModuleSize)){
    for(h in 1:length(splitSize)){
      ##############Run WCGNA Network
      ################################################################################  
      #create outputdir name
      transform_method_params=paste0(transform_i, "_ms_", minModuleSize, "_ds_", splitSize)
      net.out<-file.path(output.out, 'Network', transform_method_params)
      dir.create(net.out)
      print(paste0("Starting WCGNA network on", transform[i], 'data'))
      ms=minModuleSize[j]
      ds=splitSize[h]
      source(file.path(code.dir, 'Data_reduction', 'CreateNetwork.R'))
    }      
  }
  #No RF for PRJ
  # for(j in 1:length(duplicates)){
  #   for(h in 1:length(incident_visit)){
  #     ##############Run RF
  #     ################################################################################  
  #     print(paste0("Starting random forest on", transform[i], 'data'))
  #     dups=duplicates[j]
  #     iVisit=incident_visit[h]
  #     source(file.path(code.dir, 'Data_reduction', 'CreateRF.R'))
  #   }
  # }
  print(paste0("Done with transform:", transform_i))
}  

################################################################################  
##############Run CST 
####Please note: for this to work the DMM_flux.R script must have been run on a remote server, and the resulting
#file DMM_fit_asvtrim.rda saved in the Created data repository 
load(file.path(data.out, 'CST', 'PRJ_DMM_fit_asvtrim.rda'))
### The return value can be queried for measures of fit (Laplace, AIC, BIC); these are plotted for different k
fit<-fit_asvtrimmed
lplc <- sapply(fit, laplace)
AIC<-sapply(fit, AIC)
BIC<-sapply(fit, BIC)
pdf(file.path(output.out,  'CST', paste0(my_data, '_laplacefit.pdf')))
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
dev.off()
data=as.matrix(as.data.frame(prj.filter@otu_table))
#loop through k CST #s
for(i in 1:length(k)){
  k_i=k[i]
  source(file.path(code.dir, 'Data_reduction', 'CreateCST.R'))
}

#diversity statistics & smutans as an indicator taxa
#want to use ALL taxa for diversity statistics (pre tax filter)
#don't need smutans for PRJ
load(file.path(data.out, 'phyloseq', 'phyasv_postsamplefilter.Rdata'))
alpha<-estimate_richness(prj_samplefilter, measures=c('Shannon', 'Chao1', 'Observed'))
alpha$SampleID=rownames(alpha)
# smut<- microbiome::transform(ps_truesamples, 'compositional') %>% 
#   subset_taxa(Species=="mutans")%>%psmelt()%>%
#   mutate(smutAbundance=Abundance, smutPrevalence=ifelse(Abundance>0, "Yes", "No"))%>%
#   select(Sample, COHRAID, FoxCavID,smutAbundance, smutPrevalence)
# diversity_smutans=left_join(alpha, smut)
saveRDS(alpha, file=file.path(data.out, 'phyloseq', 'diversity.Rds'))
