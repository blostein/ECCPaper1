#split up metasqueeze SQM objects
#load metasqueeze objects
load(file.path(wgs.path, 'saliva_sqm.Rdata'))
load(file.path(wgs.path, 'plaque_sqm.Rdata'))

#Kegg
saliva_kegg_abund=saliva$functions$KEGG$abund
plaque_kegg_abund=plaque$functions$KEGG$abund
saliva_kegg_paths=saliva$misc$KEGG_paths
plaque_kegg_paths=plaque$misc$KEGG_paths
saliva_kegg_names=saliva$misc$KEGG_names
plaque_kegg_names=plaque$misc$KEGG_names

#species
saliva_species_abund=saliva$taxa$species$abund
plaque_species_abund=plaque$taxa$species$abund
saliva_species_percent=saliva$taxa$species$percent
plaque_species_percent=plaque$taxa$species$percent
#save as slim files 
save(list=c('saliva_kegg_abund', 'plaque_kegg_abund', 
       'saliva_species_percent', 'plaque_species_percent', 
       'saliva_species_abund', 'plaque_species_abund', 
       'saliva_kegg_paths', 'plaque_kegg_paths', 
       'saliva_kegg_names', 'plaque_kegg_names'), 
     file=file.path(data.out, 'MS', 'kegg_and_taxa_slim.Rdata'))

rm(saliva); rm(plaque)