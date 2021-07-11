##########################################Random Forest#################################################################
#####################################################CAVITIES####################################################################################
#################################################################################################################################################
#set seed
set.seed(49125)
#make transform directory
parameters=paste0(transform_i, '_iVisit_', iVisit, '_duplicates_', dups)
 
#chose your variables repeated cv variables
numFolds <- trainControl(method = "repeatedcv", number = num, repeats=rps, summaryFunction = twoClassSummary,
                         classProbs = TRUE, savePredictions = T) #how many folds in X validation?

#subset out taxtable
taxa<-as.data.frame(ps@tax_table)
taxa$ASV<-row.names(taxa)
taxa$SciName<-paste(taxa$Genus, ifelse(is.na(taxa$Species), taxa$ASV, taxa$Species), sep=' ')
n_taxa=length(taxa_names(ps))
########Code in this section only applicable to COHRA ECC RF

########################################################################
#metadata processing and choices
load(file.path(data.out, 'Metadata', 'MetaVisit.Rdata'))
MetaVisit=left_join(MetaVisit, ps@sam_data)
#Due to case-control incident density sampling design, some individuals
#were sampled as both cases and controls
if(dups==T){MV<-MetaVisit}
if(dups==F){MV=MetaVisit%>%filter(count!="Double; drop")}

#Filter to only samples with microbiome data
MV=MV %>%filter(!is.na(FoxCavID))

#################################################################################
#Run RF function
runRF=function(ps_object, v, iVisit){
  #ps_visit
  ps_df=ps_object%>%otu_table()%>%as.data.frame()
  if(v=='pre-incident'){ids=MV %>% filter(Visit==IncidentVisit-1)%>%pull(FoxCavID)
  ps_visit=ps_df[row.names(ps_df)%in% ids, ]
  ps_visit=left_join(MV %>% filter(FoxCavID %in% ids & Visit==IncidentVisit-1)%>%select(FoxCavID, COHRAID, CaseEver), 
                     ps_visit%>%mutate(FoxCavID=row.names(ps_visit)))}
  if(is.numeric(v)==T){ids=MV%>%filter(Visit==v)%>%pull(FoxCavID)
  ps_visit=ps_df[row.names(ps_df)%in% ids, ]
  ps_visit=left_join(MV %>% filter(FoxCavID %in% ids)%>%select(FoxCavID, COHRAID, CaseEver), 
                                           ps_visit%>%mutate(FoxCavID=row.names(ps_visit)))}
  
  #incident visit & iVisit filtering
  #incident visit
  if(v!='pre-incident'){
    incident=MV %>%filter(IncidentVisit==v & Visit==v)%>%pull(COHRAID)%>%unique()
    if(iVisit==F){
      ps_visit=ps_visit %>% filter(!(COHRAID %in% incident))
    }
  }
  ps_visit<-data.frame(case=ps_visit[, 'CaseEver'], as.matrix(ps_visit%>%select(-FoxCavID, -COHRAID, -CaseEver)))
  if(ncol(ps_visit)!=n_taxa+1){stop('incorrect taxa number')}
  #checks
  if(v!='pre-incident'){
    if(iVisit==T & dim(ps_visit)[1]!=MV %>% filter(Visit==v)%>%nrow()){stop("incorrect n!")}
    if(iVisit==F & dim(ps_visit)[1]!=MV %>% filter(Visit==v & !(COHRAID %in%incident))%>%nrow()){stop("incorrect n!")}
  }
  if(v=='pre-incident'){
    if(dim(ps_visit)[1]!=MV %>% filter(Visit==IncidentVisit-1)%>%nrow()){stop("incorrect n!")}
  }
  #run RF and evaluate 
  rfFit<-train(case ~ ., data =ps_visit, method='rf', proximity=TRUE, trControl=numFolds) #no preprocess options: don't scale or center RF data
  rfEval<-evalm(rfFit, positive='Case')
  
  #Smutans only random forest
  smut_asv=taxa%>%filter(Genus=='Streptococcus' & Species=='mutans') %>%row.names()%>%as.name()
  rfFit_SM<-eval(substitute(train(case ~ x, data =ps_visit, method='rf', proximity=TRUE, trControl=numFolds), list(x=as.name(smut_asv))))
  rfEval_SM<-evalm(rfFit_SM, positive='Case')
  
  #return
  f=list('rfFits'=list('All taxa'=rfFit, 'S. mutans only'=rfFit_SM), 'rfEvals'=list('All taxa'=rfEval, 'S. mutans only'=rfEval_SM))
  names(f$rfFits)=paste0('V', v, names(f$rfFits)); names(f$rfEvals)=paste0('V', v, names(f$rfEvals))
  return(f)
}
#################################################################################
#run RF at prefered visits
rf_12=runRF(ps, 5, iVisit)
rf_24=runRF(ps, 7, iVisit)
rf_1yr=runRF(ps, 'pre-incident', iVisit)
#all evaluations together
allEval<-evalm(list(rf_12$rfFits$`V5All taxa`, rf_24$rfFits$`V7All taxa`, 
                    rf_12$rfFits$`V5S. mutans only`, rf_24$rfFits$`V7S. mutans only`), 
               positive='Case', gnames = c("All taxa 12-month", "All taxa 24-month all taxa", 
                                           "S. mutans only 12-month", "S. mutans only 24-month"))

#regrooupd into fits and evals
rfFits<-c(rf_12$rfFits, rf_24$rfFits, rf_1yr$rfFits)
rfEvals<-c(rf_12$rfEvals, rf_24$rfEvals, rf_1yr$rfEvals)

#rename objects to reflect transform method
assign(paste0("rfFits_", transform_i), rfFits)
assign(paste0("rfEvals_", transform_i), rfEvals)
assign(paste0("rfAllEvals_", transform_i), allEval)

#save 
save(list=c(paste0("rfFits_", transform_i), paste0("rfEvals_", transform_i), paste0("rfAllEvals_", transform_i)), 
     file=paste0(data.out, '/RF/', "rfResults_", parameters, ".Rdata", sep=''))

