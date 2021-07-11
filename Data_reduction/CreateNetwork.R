#Network creation 
####################
####################Prepare data
#extract portions of ps object
my.tax<-as.data.frame(tax_table(ps))
my.tax$ASV<-row.names(my.tax)
my.tax$ScientificName<-paste(my.tax$Genus, ifelse(!is.na(my.tax$Species), my.tax$Species, my.tax$ASV), sep=" ")
Twhole.s<-as.data.frame(otu_table(ps))
###########
################

####################
####################Detect sample outliers
gsg = goodSamplesGenes(Twhole.s, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(Twhole.s), method = "average");
pdf(file.path(net.out, 'SampleTree.pdf'))
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
###########
################




####################
####################Plot soft thresholding powers
#lnames = load(file=paste0(output.dir, "TwholeAndDatatraits.Rdata"))
powers = c(c(1:10), seq(from = 11, to=30, by=1))
sft = pickSoftThreshold(Twhole.s, powerVector = powers, verbose = 5, networkType = "signed")
pdf(file.path(net.out, 'SoftPower.pdf'))
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
###########
################

###########################################################
#PAUSE, change soft power threshold depending on graph 
softPower = sft$powerEstimate;
###########################################################


####################
#####################create adjacency matrix & topological overlap matrix
adjacency = adjacency(Twhole.s, power = softPower, type = "signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM
TaxaTree = hclust(as.dist(dissTOM), method = "average");
###########
################

####################
#####################Plot Taxa Tree
sizeGrWindow(12,9)
pdf(file.path(net.out, 'TaxaTree.pdf'))
plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()
###########
################

####################
#####################Create modules
dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                            deepSplit = ds, pamRespectsDendro = FALSE, minClusterSize = ms);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
####################
#####################Plot Taxa tree with modules 
sizeGrWindow(8,6)
pdf(file.path(net.out, 'TaxaTreewModules.pdf'))
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")
dev.off()
####################
#####################Identify module eigenASVs and cluster dissimilarity
MEList = moduleEigengenes(Twhole.s, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
####################
#####################Plot clustering of eigenASVs
sizeGrWindow(7, 6)
pdf(file.path(net.out, 'eigenASVsClustering.pdf'))
plot(METree, main = "Clustering of module eigenASVs",
     xlab = "", sub = "")
dev.off()
#################################################################
#PAUSE change thresholding
MEDissThres = 0.20
#################################################################
merge = mergeCloseModules(Twhole.s, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
####################
#####################Plot Taxa Tree with merged modules
sizeGrWindow(12, 9)
plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors
####################
#####################estimate correlation between module eigene ASVs and traits of interest
ntaxa=ncol(Twhole.s)
nSamples=nrow(Twhole.s)
MEs0 = moduleEigengenes(Twhole.s, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#moduleTraitCor = cor(MEs, dataTraits, use = "p");
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
####################
#####################estimate correlation between module eigene ASVs and taxa
modNames = substring(names(MEs), 3)
TaxaModuleMembership = as.data.frame(cor(Twhole.s, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
names(TaxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
###########################################################
#############Creating data for export
#ASVs belong to which modules? 
probes = names(Twhole.s)
TaxaInfo0 = data.frame(OTU = probes,moduleColor = moduleColors)
my.tax$OTU<-my.tax$oligonumber
#get long form of PS object
pslong<- ps_compositional%>%
  psmelt() %>%                                         # melt to long format
  arrange(OTU)%>%
  mutate(ScientificName=paste(Genus, ifelse(!is.na(Species), Species, OTU), sep=" "))
#Add module data to species data
long.data.taxainf0<-left_join(pslong, TaxaInfo0)
network.data<-long.data.taxainf0 %>% select(-c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)) %>%unique()
#summarize per module abundance
if(my_data=='cohra'){
  module.data <- network.data %>% dplyr::group_by(Sample, moduleColor) %>% dplyr::mutate(modAbund=sum(Abundance)) %>% select(-Abundance) %>%unique()
  module.data<-module.data %>% select(COHRAID, FoxCavID, moduleColor, modAbund)%>%unique() #This will break with non COHRA data
}
if(my_data=='PRJ'){
  module.data <- network.data %>% dplyr::group_by(SampleID, moduleColor) %>% dplyr::mutate(modAbund=sum(Abundance)) %>% select(-Abundance) %>%unique()
  module.data<-module.data %>% select(SampleID, moduleColor, modAbund)%>%unique() #This will break with non PRJ data
}
#hubs 
hubs=chooseTopHubInEachModule(Twhole.s, dynamicColors, omitColors = FALSE)
#hubs2 = chooseTopHubInEachModule(Twhole.s, dynamicColors, omitColors = FALSE, type='signed')
#my.tax[my.tax$ASV %in% hubs2,]

#intramodular connectivity
intraConnect_a<-intramodularConnectivity(adjacency, dynamicColors)
intraConnect_a$ASV<-row.names(intraConnect_a)
intraConnect_a<-left_join(intraConnect_a, my.tax)
intraConnect_a$module=moduleColors

#name networks based on hub and top abundant taxa 
network_names=left_join(
  #hubs
  intraConnect_a%>%group_by(module)%>%slice_max(order_by=kWithin, n=2)%>%dplyr::mutate(id=paste0("CentralTaxa", row_number()))%>%
    select(module, ScientificName, id)%>% pivot_wider(id_cols=module, names_from=id, values_from=ScientificName)%>%
    dplyr::mutate(moduleColor=module, module_label=paste0('Central taxa:', CentralTaxa1, ' + ', CentralTaxa2)),
  #top taxa
  long.data.taxainf0%>%group_by(moduleColor, ScientificName)%>%dplyr::summarise(mean=mean(Abundance))%>%group_by(moduleColor)%>%
    slice_max(order_by=mean, n=2)%>%select(ScientificName, moduleColor)%>%dplyr::mutate(id=paste0('TopTaxa', row_number()))%>%
    pivot_wider(id_cols=moduleColor, values_from=ScientificName, names_from=id)) %>% 
  dplyr::mutate(module_label2= paste0(TopTaxa1, ' & ', TopTaxa2, ' network'), module_label3= paste0(module_label2, '\n (central taxa: ', CentralTaxa1, ')'))

#add network names to module data
module.data=left_join(module.data, network_names)

#record parameters: 
all_params=c(transform_i, ms, ds, softPower, MEDissThres)
parameters=paste0(transform_i, "_ms_", ms, '_ds_', ds) 
names(all_params)=c("Transform method", "min Module size", "deep split number", "soft power", "thresholding")
params_df=as.data.frame(all_params); colnames(params_df)=parameters

#group data into network object
network_data<-list("adj"=adjacency, "intraStats"=intraConnect_a, "module.data"=module.data, 
               "hubs"=hubs, 'taxaColors'=moduleColors, 'networkNames'=network_names,
               "moduleNum"=moduleColors%>%unique()%>%length(), 
               "params"=params_df)

assign(paste0(transform_i, '_network'), network_data)
###########################################################
#############ExportData

save(list=c(paste0(transform_i, '_network')), file=file.path(data.out, 'Network', paste0(parameters, '.Rdata')))

