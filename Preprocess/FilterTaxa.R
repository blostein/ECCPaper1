##########################################filter taxa#################################################################
#####################################################CAVITIES####################################################################################
#################################################################################################################################################
##################################
#load chosen phyloseq object  & rename
#load(file.path(data.out, 'phyloseq', "phyasv_truessamples.Rdata"))
ps<-ps_truesamples

#chose variables 
prev.filter<-round(prev*nrow(ps@sam_data),0)

#table of ASVs at phyla level
#table(tax_table(ps)[, "Phylum"], exclude = NULL)
#223 are NA at the Phyla level, should be filtered

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

#plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#curious about Cyanobacteria
#prevdf %>% filter(Phylum=="Cyanobacteria")
#Arthospira -> Spironela? weird

# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !is.na(Phylum))

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
prevVsabun_g<-ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
#no obvious separation point - note Cyanobacteria all  low abundance and fairly low prevalence

#determine which ASVs have low prevalence 
lowprev<-prevdf1 %>% filter(Prevalence<prev.filter)
lowprev$OTU = rownames(lowprev)

#create RA dataset
ps_ra<-transform_sample_counts(ps1, function(OTU) OTU/sum(OTU)) %>% psmelt()
#identify the maximum RA for low prevalence ASVs & filter to low abundance (below filter) ASVs
lowabundance.lowprev<-ps_ra %>% filter(OTU %in% rownames(lowprev))%>% dplyr::group_by(OTU)%>% dplyr::mutate(maxAbund = max(Abundance)) %>% filter(maxAbund<abund.filter) %>% ungroup()%>%dplyr::select(OTU)%>%unique()
#identify the maximum RA for low prevalence ASVs & filter to high abundance (above filter) ASVs
hiabundance.lowprev<-ps_ra %>% filter(OTU %in% rownames(lowprev))%>% dplyr::group_by(OTU)%>% dplyr::mutate(maxAbund = max(Abundance)) %>% filter(maxAbund>abund.filter) %>% ungroup()%>%dplyr::select(OTU)%>%unique()
#not low prev 
hiprev<-ps_ra %>% filter(!(OTU %in% rownames(lowprev))) %>% dplyr::select(OTU) %>% unique()

#check that your numbers add up 
check1<-dim(lowabundance.lowprev)[1]+dim(hiabundance.lowprev)[1]+dim(hiprev)[1]==dim(ps1@otu_table)[2]

#examine hiabundance.lowprev otus 
ps_hiAloP<-ps_ra %>% filter(OTU %in% hiabundance.lowprev$OTU)
#ps_hiAloP %>% ggplot(., aes(x=OTU, y=Abundance))+geom_violin()+facet_wrap(~Genus, scales='free_x')

#make vector of otus to keep
keepTaxa= rownames(prevdf1)[!(row.names(prevdf1) %in% lowabundance.lowprev$OTU)]
ps.filter = prune_taxa(keepTaxa, ps1)
ps.filter #274 taxa 

#check that numbers add up 
check2<-dim(hiabundance.lowprev)[1]+dim(hiprev)[1]==dim(ps.filter@otu_table)[2]

#Warn user if either check failed:
if(check1==F | check2==F){print("WARNING: ASV losses don't add up!")}

#save ps object
save(ps.filter, file=file.path(data.out, 'phyloseq', "phyasv_postTax.Rdata"))
