
#all BabySubjectIDs & Unique Person IDs for tracking 
unique_kids_n=MetaVisit %>%select(BabySubjectID)%>%unique()%>%nrow()
unique_ids_n=MetaVisit%>%select(UniquePersonID)%>%unique()%>%nrow()
all_ids=MetaVisit%>%pull(UniquePersonID)%>%unique()
all_babyids=MetaVisit%>%pull(BabySubjectID)%>%unique()
all_ids_pa=MetaVisit%>%filter(Site=='PA')%>%pull(UniquePersonID)%>%unique()
all_babyids_pa=MetaVisit%>%filter(Site=='PA')%>%pull(BabySubjectID)%>%unique()
site_df=MetaVisit%>%select(UniquePersonID, Site)%>%unique()%>%group_by(Site)%>%dplyr::summarise(n=n())
site_df_person=MetaVisit%>%select(BabySubjectID, Site)%>%unique()%>%group_by(Site)%>%dplyr::summarise(n=n())

#missing metadata by visit
unique_id_byvisit=MetaVisit %>% filter(Visit<=IncidentVisit2)%>% group_by(Visit)%>%select(Visit, UniquePersonID)%>%split(f=.$Visit)
unique_babyid_byvisit= MetaVisit %>% filter(Visit<=IncidentVisit2)%>% group_by(Visit)%>%select(Visit, BabySubjectID)%>%split(f=.$Visit)

#incident visit ids
incident_id_byvisit=MetaVisit %>% filter(Visit==IncidentVisit2)%>% group_by(Visit)%>%select(Visit, UniquePersonID)%>%split(f=.$Visit)

#by visit by case status regardless of if has saliva sample 
ids_notin=lapply(unique_id_byvisit, function(x) data.frame('UniquePersonID'=all_ids[!(all_ids %in% x$UniquePersonID)]))
ids_notin=lapply(ids_notin, function(x) x%>%mutate(IncidentVisit=gsub('.*_', '', x$UniquePersonID)))
ids_notin_df=bind_rows(ids_notin, .id='Visit')%>%filter(as.integer(Visit)<as.integer(IncidentVisit))
ids_notin_n=ids_notin_df%>%group_by(Visit) %>% filter(!(Visit %in% c(2, 4)))%>%dplyr::summarise(n=n())%>%mutate(Case='Missing metadata', Type=Case, Visit=as.integer(Visit))

#by visit by case status regardless of if has saliva sample no duplicates 
ids_notin_b=lapply(unique_babyid_byvisit, function(x) data.frame('BabySubjectID'=all_babyids[!(all_babyids %in% x$BabySubjectID)]))
MVnodups=MetaVisit %>% filter(count!="Double; drop")%>%select(BabySubjectID, IncidentVisit2)%>%unique()
ids_notin_b=lapply(ids_notin_b, function(x) x%>%left_join(MVnodups))
ids_notin_b_df=bind_rows(ids_notin_b, .id='Visit')%>%filter(as.integer(Visit)<as.integer(IncidentVisit2))
ids_notin_b_n=ids_notin_df%>%group_by(Visit) %>% filter(!(Visit %in% c(2, 4)))%>%dplyr::summarise(n=n())%>%mutate(Case='Missing metadata', Type=Case, Visit=as.integer(Visit))


#by visit by case status in visit 2 and 4 (PA only visit)
ids_notin_PAonly=lapply(unique_id_byvisit, function(x) data.frame('UniquePersonID'=all_ids_pa[!(all_ids_pa %in% x$UniquePersonID)]))
ids_notin_PAdf=bind_rows(ids_notin_PAonly, .id='Visit')%>%filter(as.integer(Visit)%in% c(2, 4))%>%
                group_by(Visit)%>%dplyr::summarise(n=n())%>%mutate(Case='Missing metadata', Type=Case, Visit=as.integer(Visit))

#and again no dups
ids_notin_b_PAonly=lapply(unique_babyid_byvisit, function(x) data.frame('UniquePersonID'=all_babyids_pa[!(all_babyids_pa %in% x$BabySubjectID)]))
ids_notin_b_PAdf=bind_rows(ids_notin_PAonly, .id='Visit')%>%filter(as.integer(Visit)%in% c(2, 4))%>%
  group_by(Visit)%>%dplyr::summarise(n=n())%>%mutate(Case='Missing metadata', Type=Case, Visit=as.integer(Visit))


#count of cases and controls with meta data at each visit
hasMeta<-MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2)%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case, fill = list(n = 0))

hasMeta_nodup<-MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2 & count!="Double; drop")%>%
  select(BabySubjectID, Visit, Case)%>%
  unique()%>%group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case, fill = list(n = 0))

#count of incident cases and controls 
incidentCount<-MetaVisit %>% filter(Visit==IncidentVisit2)%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case, fill = list(n = 0))

incidentCount_nodups<-MetaVisit %>% filter(Visit==IncidentVisit2 & count!='Double; drop')%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case, fill = list(n = 0))

incidentCount_w16S<-MetaVisit %>% filter(Visit==IncidentVisit2 & HasCleanSample==1)%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case, fill = list(n = 0))

incidentCount16Snodups<-MetaVisit %>% filter(Visit==IncidentVisit2 & HasCleanSample==1 & count!='Double; drop')%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case, fill = list(n = 0))

all_data_byvisit=rbind(ids_notin_n, ids_notin_PAdf, hasMeta%>%mutate(Type='Attended visit'))
all_data_byvisit_nodups=rbind(ids_notin_b_n, ids_notin_b_PAdf, hasMeta_nodup %>% mutate(Type='Attended visit'))
#check that this all adds up to total unique visit ids
check_data=left_join(all_data_byvisit%>%group_by(Visit)%>%
                       dplyr::summarise(n_tot=sum(n)), incidentCount%>%
                       group_by(Visit)%>%dplyr::summarise(incident_n=sum(n)))%>%
            mutate(new_n=lag(n_tot-incident_n))
check_data%>%filter(n_tot!=new_n)

check_data_nodups=left_join(all_data_byvisit_nodups%>%group_by(Visit)%>%
                              dplyr::summarise(n_tot_nodups=sum(n)), incidentCount_nodups%>%
                              group_by(Visit)%>%dplyr::summarise(incident_n_nodups=sum(n)))%>%
  mutate(new_n_nodups=lag(n_tot_nodups-incident_n_nodups))

#has metadata
sample_row_meta<-all_data_byvisit%>%filter(Case!='Missed visit/no visit data')%>%group_by(Case)%>%mutate(sum=sum(n))%>%
  pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%relocate(Total, .after = last_col())

#missing metadata
sample_row_nometa<-all_data_byvisit%>%filter(Case=='Missed visit/no visit data')%>%group_by(Case)%>%
  pivot_wider(names_from=Visit, values_from=n)

#has 16S clean sample
has16S=MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2 &HasCleanSample==1)%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
  group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case,
           fill = list(n = 0))%>%group_by(Case)%>%mutate(sum=sum(n))
sample_row_16S<-has16S%>%
  pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%relocate(Total, .after = last_col())

has16S_nodup=MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2 &HasCleanSample==1 & count!='Double; drop')%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
  group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case,
           fill = list(n = 0))%>%group_by(Case)%>%mutate(sum=sum(n))

sample_row_16S<-has16S%>%
  pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%relocate(Total, .after = last_col())


#missing 16S clean sample
no16S=MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2 & is.na(HasCleanSample))%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
  group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case,
           fill = list(n = 0))%>%group_by(Case)%>%mutate(sum=sum(n))
sample_row_no16S=no16S%>%
  pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%relocate(Total, .after = last_col())

no16S_nodups=MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2 & is.na(HasCleanSample &count!='Double; drop'))%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
  group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case,
           fill = list(n = 0))%>%group_by(Case)%>%mutate(sum=sum(n))
sample_row_no16S=no16S%>%
  pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%relocate(Total, .after = last_col())


#has 16S clean sample no dups 
sample_row_16S_nodups<-MetaVisit %>% filter(Visit<=10 & Visit<=IncidentVisit2 &HasCleanSample==1)%>%filter(count!="Double; drop")%>%
  select(UniquePersonID, Visit, Case)%>%
  unique()%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
  group_by(Visit, Case)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, Case,
           fill = list(n = 0))%>%group_by(Case)%>%mutate(sum=sum(n))%>%
  pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%relocate(Total, .after = last_col())

sample_row_clean_both=paste0(sample_row_16S, ifelse(sample_row_16S!=sample_row_16S_nodups, paste0('\n(', sample_row_16S_nodups, ')'), ''))

#incident row
incident_row<-incidentCount%>%group_by(Case)%>%
  mutate(sum=sum(n))%>%pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%
  select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%
  relocate(Total, .after = last_col())
incident_rownodups=incidentCount_nodups%>%group_by(Case)%>%
  mutate(sum=sum(n))%>%pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%
  select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%
  relocate(Total, .after = last_col())
#incident row with 16S 
incident_row16S<-incidentCount_w16S%>%group_by(Case)%>%
  mutate(sum=sum(n))%>%pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%
  select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%
  relocate(Total, .after = last_col())
#incident row 16S no dups
incident_row16Snodups= incidentCount16Snodups%>%group_by(Case)%>%
  mutate(sum=sum(n))%>%pivot_wider(names_from=Case, values_from=c(n, sum))%>%
  mutate(label=paste0(n_Case, '/', n_Control), Total=paste0(sum_Case, '/', sum_Control))%>%
  select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%
  relocate(Total, .after = last_col())

#incident row both 16S and no 16S
incident_row_both16S<-c('', '', '', paste0(incident_row16S, ifelse(incident_row16S!=incident_row16Snodups, paste0(' (', incident_row16Snodups, ')'), '')))
incident_row_both<-c('', '', '', paste0(incident_row, ifelse(incident_row!=incident_rownodups, paste0(' (', incident_rownodups, ')'), '')))

#wgs data
sample_row_wgs<-species_info$sample_table%>%select(COHRAID, Visit, CaseStatus)%>% unique() %>% mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
  group_by(Visit, CaseStatus)%>%
  dplyr::summarise(n=n())%>%ungroup() %>%
  complete(Visit, CaseStatus,
           fill = list(n = 0))%>%group_by(CaseStatus)%>%mutate(sum=sum(n))%>%
  pivot_wider(names_from=CaseStatus, values_from=c(n, sum))%>%
  mutate(label=paste0(n_case, '/', n_control), Total=paste0(sum_case, '/', sum_control))%>%
  select(Visit, label, Total)%>%pivot_wider(names_from=Visit, values_from=label)%>%
  relocate(Total, .after = last_col())

#table
metadataTab=all_data_byvisit%>%select(-Case)%>%group_by(Visit, Type)%>%dplyr::summarise(n=sum(n))%>%pivot_wider(id_cols=Visit, names_from=Type, values_from=n)%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))
metadataTab_nodups=all_data_byvisit_nodups%>%select(-Case)%>%group_by(Visit, Type)%>%dplyr::summarise(n=sum(n))%>%pivot_wider(id_cols=Visit, names_from=Type, values_from=n)%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))

sixteenSTab=left_join(left_join(has16S%>%group_by(Visit)%>%dplyr::summarise(`Has clean 16S saliva sample`=sum(n)), 
                      no16S%>%group_by(Visit)%>%dplyr::summarise(`Missing saliva sample/failed 16S QC`=sum(n))), 
              incidentCount_w16S%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
                group_by(Visit)%>%dplyr::summarise(`incident N with 16S`=sum(n)))

sixteenSTab_nodups=left_join(left_join(has16S_nodup%>%group_by(Visit)%>%dplyr::summarise(`Has clean 16S saliva sample (no duplicates)`=sum(n)), 
                                no16S_nodups%>%group_by(Visit)%>%dplyr::summarise(`Missing saliva sample/failed 16S QC (no duplicates)`=sum(n))), 
                      incidentCount16Snodups%>%mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
                        group_by(Visit)%>%dplyr::summarise(`incident N with 16S (no duplicates)`=sum(n)))

lostSamplesTab=left_join(check_data%>%
                           mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), 
                                               labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
                           select(-new_n), 
                         left_join(metadataTab, sixteenSTab))
lostSampleTab_nodup=left_join(check_data_nodups%>%
                                mutate(Visit=factor(Visit, levels=c(2, 3, 4, 5, 7, 8, 9, 10), 
                                                    labels=c('Birth', '3', 'First tooth', '12', '24', '36', '48', '60')))%>%
                                select(-new_n_nodups), 
                              left_join(metadataTab_nodups, sixteenSTab_nodups))%>%mutate(`Attended visit (no duplicates)`=`Attended visit`, `Missing metadata (no duplicates)`=`Missing metadata`)%>%select(-`Missing metadata`, -`Attended visit`)

lostSamples=left_join(lostSamplesTab, lostSampleTab_nodup)%>%select(Visit, n_tot, n_tot_nodups, incident_n, incident_n_nodups, `Missing metadata`, 
                                                                    `Missing metadata (no duplicates)`, `Attended visit`, `Attended visit (no duplicates)`, 
                                                                    `Missing saliva sample/failed 16S QC`, `Missing saliva sample/failed 16S QC (no duplicates)`, 
                                                                    `Has clean 16S saliva sample`, `Has clean 16S saliva sample (no duplicates)`, 
                                                                    `incident N with 16S`, `incident N with 16S (no duplicates)`)
colnames(lostSamples)=c('Visit', 'Total N (unique records)', 'Total N (unique children, no duplicates)', 'Incident N', 'Incident N (excluding duplicates)', colnames(lostSamples)[6:15])
#make histogram figure


sample_table2<-rbind(incident_row_both, sample_row_clean_both, sample_row_wgs)
colnames(sample_table2)<-c('Birth', '  2 mos ', 'First tooth', '   12 mos', '   24 mos ', '  36 mos*', '  48 mos ', '  60 mos ', '  Total*')
rownames(sample_table2)<-c('Incident cases/controls',  'Case/control samples \n with 16S data', 'Case/control samples \n with WGS data')



sample_hist<-MetaVisit %>% 
  mutate(IncidentVisit=ifelse(IncidentVisit==5, 6, IncidentVisit2))%>%
  select(UniquePersonID, IncidentVisit, Case)%>%
  unique()%>%
  ggplot(aes(x=IncidentVisit, fill=Case)) +
  #geom_density(alpha=.2, adjust=1.5)+scale_x_continuous("Age (months)", breaks = c( 3, 4, 5, 6, 7, 8, 9, 10), labels=c('Birth', "3", 'First tooth', "12", "24", "37", "49", "60"), limits = c(3, 11))+
  geom_histogram(aes(y = ..density../5), position='dodge')+
  scale_y_continuous("Counts", breaks = c(0, 0.1, 0.2, 0.3, 0.4), labels=round(c(0, 0.1*94, 0.2*94, 0.3*94, 0.4*94)))+
  scale_fill_manual(values=c("darkorange3", "royalblue3"), labels=c("ECC case", "Control"))+theme_classic()+
  theme(text=element_text(size = 24))+theme(legend.text=element_text(size=12))+theme(legend.position = 'left')+
  theme(rect = element_rect(fill = "transparent") # all rectangles
  )+theme(axis.text.x=element_blank())+theme(axis.title.x = element_blank())+theme(plot.margin = unit(c(2, 0, -1.5, 0), "cm"))


samplehist2<-MetaVisit %>% 
  mutate(IncidentVisit=ifelse(IncidentVisit==5, 6, IncidentVisit2))%>%
  #filter(Visit==IncidentVisit2 & HasCleanSample==1)%>%
  filter(Visit==IncidentVisit2)%>%
  select(UniquePersonID, IncidentVisit, Case)%>%
  unique()%>%
  ggplot(aes(x=IncidentVisit, fill=Case)) +
  scale_x_continuous("Age (months)", breaks = c( 3, 4, 5, 6, 7, 8, 9, 10), labels=c('Birth', "3", 'First tooth', "12", "24", "37", "49", "60"), limits = c(3, 11))+
  scale_y_continuous("Incident visit counts")+
  geom_histogram(position='dodge')+
  scale_fill_manual(values=c("darkorange3", "royalblue3"), labels=c("ECC case", "Control"))+theme_classic()+
  theme(text=element_text(size = 24))+theme(legend.text=element_text(size=12))+theme(legend.position = 'left')+
  theme(rect = element_rect(fill = "transparent") # all rectangles
  )+theme(axis.text.x=element_blank())+theme(axis.title.x = element_blank())+theme(plot.margin = unit(c(2, 0, -1.5, 0), "cm"))

# d_tooth_file<-file.path("~/Downloads/noun_Decayed Tooth_970447.png")
# h_tooth_file<-file.path("~/Downloads/noun_Tooth_970480.png")
# 
# 
# myt <- ttheme_minimal(core = list(fg_params=list(hjust = 1, x=1, cex=1.3),                 bg_params=list(fill=c("lightgrey", "grey"))), 
#                       colhead=list(fg_params=list(cex=1.3)),rowhead=list(fg_params=list(cex=1.3)))
# 
# 
# t1<-tableGrob(sample_table2, theme=myt)
# t1<-gtable_add_grob(t1, grobs = rectGrob(gp = gpar(fill = NA, lwd = 4)),t = 1, b = 1, l = 2, r = ncol(t1))
# 
# x1=0.735
# y1= 1.03; y2=y1-0.055; y3=y1-0.07; y4=y1-0.12
# h1<-plot_grid(sample_hist, t1, ncol=1)
# h1a<-plot_grid(samplehist2, t1, ncol=1)
# h2<-h1+draw_image(d_tooth_file, x = x1, y = y1, hjust = 1, vjust = 1, width = 0.07, height = 0.15)+draw_image(h_tooth_file, x = x1, y = y2, hjust = 1, vjust = 1, width = 0.07, height = 0.15)+annotate("segment", x = x1-0.05, xend = 0.3, y = y3, yend = y3, colour = "darkorange3", size=1, alpha=0.6, arrow=arrow())+annotate("segment", x = x1-0.05, xend = 0.3, y = y4, yend = y4, colour = "royalblue3", size=1, alpha=0.6, arrow=arrow())
# h3<-h2+annotate('text', x=x1, y=y1-0.11, label='1) For ECC cases at each visit, \n (incident visit) select subsample of \n children who remained caries-free \n at that age (controls)', hjust=0)+annotate('text', x=0.3, y=0.84, label='2) Sequence (16S-V4) all available \n prior saliva samples up to and \n including the incident visit', hjust=0)+annotate('text', x=x1+0.03, y=y1-0.3, label='3) For subsample of cases & \n controls, sequence (WGS) \n       incident visit saliva \n           & plaque samples', hjust=0)


