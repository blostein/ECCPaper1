#################Make a tidy graph object with all the information you could want
#to run need:
  #a network object from CreateNetwork (nw, minimal)
  #a list of taxa you want to label (topTax)
  #a list of genera you want to color by (topGenus)
  #a df with additional per taxa variables containing a column 'name' which is the OTU/ASV #
  #here called 'my_df'

#Labels --> sciname notation + filter label to only be taxa of interest
nw$intraStats=nw$intraStats%>%
    #first shorten to scientific notation
    mutate(ScientificNameLabel=paste(ifelse(is.na(nw$intraStats$Species), 
                                               nw$intraStats$Genus, 
                                               paste0(substring(nw$intraStats$Genus, 1, 1), '.')), 
                                        ifelse(is.na(nw$intraStats$Species), 
                                               nw$intraStats$ASV, 
                                               nw$intraStats$Species), sep=' '))%>%
  #pull which taxa to display in network graph (topTax)  
  mutate(ScientificNameDisplay=ifelse(ScientificName%in%topTax, 
                                        ScientificNameLabel, ''))%>%
  #pull which genera to color in netwrok graph (topGenus)
  mutate(GenusLabel=factor(ifelse(Genus %in% topGenus, Genus, 'Other'), 
                             levels=c(topGenus, 'Other')))

#convert to igraph formats
g_whole<-graph_from_adjacency_matrix(nw$adj, mode = "undirected", weighted = TRUE,
                                     diag = FALSE, add.colnames = NULL, 
                                     add.rownames = NA)
g2a<-igraph::as_data_frame(g_whole, what='both')

#add vertex info
g2a$vertices<-g2a$vertices%>%
  left_join(nw$intraStats, c('name'='ASV'))%>% # add Scientific labels, network stats
  left_join(nw$networkNames)%>% #add module labels
  left_join(my_df)#add other information about taxa 

#covert graph format 
g2<-graph_from_data_frame(g2a$edges, directed=F, vertices = g2a$vertices)
gg_whole<-as_tbl_graph(g2)
rm(nw) #remove temporary network object
