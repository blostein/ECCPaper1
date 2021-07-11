#A script to make Humann3 data palatable for R manipulation & ggplot graphing 
metag_metashort=metag_metashort%>% mutate_at(vars(c(BabySubjectID, Visit)), as.integer)

#Species data
#load in data
species_level=read.csv(file.path(humann3.path, 'merged_abundance_table_species.txt'), sep="\t")
#change row and column names
colnames(species_level)<-sub('^(.*?_.*?)_.*', "\\1", colnames(species_level))
row.names(species_level)<-species_level$clade_name
#split into info and abundance data
species_info<-species_level %>%select(clade_name, NCBI_tax)
species_level<-species_level %>% select(-clade_name, -NCBI_tax)
#pivot and dataframe
species_level<-data.frame(t(as.matrix(species_level)))
#create species_info object
species_info=list('otu_table'=species_level, 'taxa_table'=species_info, 'sample_table'=metag_metashort)
#pivot long and add metadata
species_info$species_long=left_join(species_info$otu_table%>%rownames_to_column(var='SampleName')%>%pivot_longer(cols=-SampleName, names_to='Species', values_to='WGSAbundance'), species_info$sample_table)

#path level data
path_level=read.csv(file=file.path(humann3.path, "ECC_pathabundCPM_unstratified.tsv"), sep="\t")
path_strat=read.csv(file=file.path(humann3.path, "ECC_pathabund_stratCPM.tsv"), sep="\t")
path_coverage=read.csv(file=file.path(humann3.path, "ECC_pathcoverage.tsv"), sep="\t")

pathdata=list('Unstratified'=path_level, 'Stratified'=path_strat, 'Coverage'=path_coverage)
#check that column names are the same
colnames(path_strat)==colnames(path_level)
new_names=sub('^(.*?_.*?)_.*', "\\1", colnames(path_level))
#change column and row names
pathdata=lapply(pathdata,  setNames, new_names)
pathdata=lapply(pathdata, function(x) pivot_longer(x, cols= -X..Pathway,names_to="SampleName", values_to = "Abundance_CPM"))
#species name for stratified only
for(i in c(2, 3)){
  pathdata[[i]]$Species=gsub(".*\\|", '', pathdata[[i]]$X..Pathway)
  pathdata[[i]]$X..Pathway=gsub("\\|.*", '', pathdata[[i]]$X..Pathway)
}

# pathd description and path id 
pathdata=lapply(pathdata, function(x) x%>%
                  mutate(PathID=gsub(":.*", "", X..Pathway), 
                         PathDescription=gsub(".*: ", "", X..Pathway))%>%
                  left_join(metag_metashort))
#add metadata to object
pathdata[['sample_table']]=metag_metashort; 
#rename pathdata$Coverarage$Abundance to coverage
pathdata$Coverage=pathdata$Coverage %>%mutate(Coverage=Abundance_CPM)%>%select(-Abundance_CPM)

#save Humann3 objects 
dir.create(humann3.out)
save(list=c('species_info', 'pathdata'), file=file.path(humann3.out, 'Humann3.Rdata'))
