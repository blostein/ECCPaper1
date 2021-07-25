library(psych)
library(plyr)
library(dplyr)
library(stringr)
library(chron)
library(lubridate)
library(naniar)
library(compareGroups)
#function for na recoding 
na_codes <- function(x, ...) {
  x[x %in% c(...)] <- NA
  x
}

#############Read in datasets###############################################################
USCUR<-read.csv(file.path(meta.data, "USCURPostnatal.csv"))
Postnatal<-read.csv(file.path(meta.data, "DFTandPostnatal.csv"))
Prenatal<-read.csv(file.path(meta.data, "Prenatal.csv"))
#Incident Visit info
ControlPrime<-read.csv(file.path(meta.data,"ControlPrime.csv"))
CasePrime<-read.csv(file.path(meta.data, "CasePrime.csv"))
colnames(CasePrime)<-c("BabySubjectID", "IncidentVisit"); CasePrime$Case<-"Case"
colnames(ControlPrime)<-c("BabySubjectID", "IncidentVisit"); ControlPrime$Case<-"Control"
IncidentVisits<-rbind(CasePrime, ControlPrime)

#############Prenatal###############################################################
#Mother Birth Date, Baby birth date to Date format
Prenatal <-Prenatal %>% mutate_at(vars(MotherBirthDate, BabyBirthdate), funs(as.Date(., format="%m/%d/%y")))
#categorical variables
Prenatal$Race<-ifelse(str_detect(as.character(Prenatal$Race), '1'), "White", "Not White")
Prenatal$BabySex<-recode(Prenatal$BabySex, `1`="Male", `2`="Female", .default=NA_character_)
Prenatal$Delivery<-recode(Prenatal$Delivery, `1`="Vaginal", `2`="C-section", .default=NA_character_)
Prenatal$Education_HS<-ifelse(Prenatal$Education<0, "Missing", ifelse(Prenatal$Education<5, "High school degree or less", "Associates degree or higher"))
Prenatal$Education_Cat<-with(Prenatal, ifelse(Education<0, "Missing", ifelse(Education<5, "High school or less", ifelse(Education<6, "Associates degree", ifelse(Education==6, "Undergraduate degree", ifelse(Education>6, "Graduate degree", NA))))))
Prenatal$Education<-recode(Prenatal$Education, `1`="8th grade or less", `2`="9th-12th grade, no diploma", `3`="High school graduate/GED", `4`="Some college credit, no degree", `5`="Associate degree", `6`="Bachelor's degree", `7`="Master's degree", `8`="Doctorate/professional degree")
Prenatal$MomsAgeatBirth<-with(Prenatal, as.period(interval(start = MotherBirthDate, end = BabyBirthdate)))$year
Prenatal<-Prenatal %>% rename_at(.vars=c("Prim_Tot_Teeth_Present", "Prim_swt", "Prim_d1ft", "Prim_d2ft", "PERM_TOT_TEETH_PRESENT", "PERM_SWT", "PERM_D1MFT", "PERM_D2MFT"), function(x) paste0("MomV1_", x))
Prenatal_trimmed = Prenatal %>% select(MotherSubjectID, BabySubjectID, Site, Race, MomV1_PERM_TOT_TEETH_PRESENT, MomV1_PERM_SWT, MomV1_PERM_D1MFT, MomV1_PERM_D2MFT, Education, BabyBirthdate, BabySex, Delivery)

######################Postnatal: Double counting###########################################################
#add in information about cases and controls - this creates two records for any observations that are:
#a) a case and a control (n=1)
#b) selected as a control twice (n=1)
#because there are now duplicate cohraids (subjectid-visit number) I create a unique id (subjectid_visit_{ca|co}_incidentvisit)
#where ca is case and co is control 
#I also create flags for easy removal of these observations as necessary
Postnatal<-left_join(Postnatal, IncidentVisits)
Postnatal$UniqueID<-paste(Postnatal$BabySubjectID, Postnatal$Visit, substr(Postnatal$Case, 1, 2), Postnatal$IncidentVisit, sep="_")
Postnatal$UniquePersonID<-paste(Postnatal$BabySubjectID, substr(Postnatal$Case, 1, 2), Postnatal$IncidentVisit, sep="_")
Postnatal$COHRAID<-paste(Postnatal$BabySubjectID, Postnatal$Visit, sep="-")

#Identify double records
Doubles<-Postnatal %>% group_by(COHRAID)%>% select(COHRAID, UniqueID) %>% dplyr::summarise(n=n()) %>%filter(n>1)%>%dplyr::mutate(BabySubjectID=gsub('-.*', '', COHRAID))%>%pull(BabySubjectID)%>%unique()
Doubles<-Postnatal %>% filter(BabySubjectID %in% Doubles)%>%select(BabySubjectID, IncidentVisit, Case, UniquePersonID)%>%unique()%>%arrange(BabySubjectID, Case, -IncidentVisit)%>%group_by(BabySubjectID)%>%dplyr::mutate(count=row_number())%>%ungroup()
#Since doubles is arranged by BabySubjectID, then Case, then (reverse) Incident Visit, we can be assured that by selecting
#the first row of any babysubjectID we will select case over controls or the last incident selection of control otherwise, just by selecting the first sample
Postnatal <- left_join(Postnatal, Doubles) %>% dplyr::mutate(count=case_when(count==1~'Double; keep', count==2~"Double; drop", is.na(count)~"Single"))
#to delete the double counted rows, can run this: Postnatal %>% filter(count !='Doubles; drop') 

######################Postnatal: Make new variables###########################################################
#Create a Decay Yes/no variable
Postnatal$Decay<-ifelse(Postnatal$Prim_d1ft>0, "Decay present (d1ft>0)", ifelse(Postnatal$Prim_d1ft<0, "No teeth", ifelse(Postnatal$Prim_d1ft==0, "No decay (d1ft=0)", NA)))
#Create a Case/Not Case Variable
EverDecayIDS<-Postnatal%>%filter(Decay=="Decay present (d1ft>0)")%>% pull(BabySubjectID)
Postnatal<-Postnatal %>% dplyr::mutate(CaseStatus=case_when(Case=="Case" & Visit==IncidentVisit~"Case (incident visit)", 
                                                     Case=="Control" & Visit==IncidentVisit~"Control (incident visit)", 
                                                     Case=="Case" & Visit<IncidentVisit~"Case (pre-incident visit)", 
                                                     Case=="Control" & Visit<IncidentVisit~"Control (pre-incident visit)", 
                                                     Case=="Control" & Visit>IncidentVisit~"Control (post-incident visit)", 
                                                     Case=="Case" & Visit>IncidentVisit~"Case (post-incident visit)"))
Postnatal$CaseEver<-with(Postnatal, ifelse(Case=="Control", "Control", "Case"))
#add in site variable
Postnatal<-left_join(Postnatal, Prenatal %>% select(BabySubjectID, MotherSubjectID, Site))

#Date variables
Postnatal$VisitDate<-as.Date(Postnatal$VisitDate, format="%m/%d/%y")
Postnatal$LastAteDate<-as.Date(Postnatal$LastAteDate, format="%m/%d/%y")
#Time since last ate
Postnatal <-Postnatal%>% mutate_at(vars(VisitDate, LastAteDate), funs(as.Date(., format="%m/%d/%y"))) %>% mutate_at(vars(VisitTime, LastAteTimeTime), funs(paste0(sub("^(.*?:.*?):.*", "\\1", .), ":00")))%>%
  dplyr::mutate(VisitDateTime=ymd_hms(paste(VisitDate, VisitTime)), AteDateTime=ymd_hms(paste(LastAteDate, LastAteTimeTime)), TimeSinceAte=as.duration(AteDateTime%--% VisitDateTime))%>%
     dplyr::mutate(AteFlag=case_when(TimeSinceAte<0 | TimeSinceAte >60*60*48 ~ 'Problematic date times', TimeSinceAte<60*60*2~'<2 hrs ago', TimeSinceAte>=60*60*2~'>=2 hrs ago'), 
            AteFlag1hr=case_when(TimeSinceAte<0 | TimeSinceAte >60*60*48 ~ 'Problematic date times', TimeSinceAte<60*60*1~'<1 hrs ago', TimeSinceAte>=60*60*1~'>=1 hrs ago'))
#Possible mistakes in recording
#'last ate' at a time after visit time
#Postnatal %>% filter(TimeSinceAte<0)%>%select(BabySubjectID, Visit, VisitDateTime, AteDateTime, TimeSinceAte)#19
#'last ate >2 days ago
#Postnatal %>% filter(TimeSinceAte>60*60*48)%>%select(BabySubjectID, Visit, VisitDateTime, AteDateTime, TimeSinceAte)#11

#numeric variables: ageatexamyears & months, salivapH prim_swt prim_d12ft, prim_d2ft, prim_Tot_Teeth_Present
#PERM_SWT PERM_D1MFT, PERM_D2MFT, PERM_TOT_TEETH_PRESENT SMSdfs PrimTotalSurf 
#ProportionSMSdfs SMSDMFS PERMTotalSurf PRroportionSMSDMFS BabyWeigth Lbs
#BabyWeightOz, BabyWeightKg, BabyHeightFt, BabyHeightIn, BabyHeightCm
Postnatal<-Postnatal %>% replace_with_na_all(condition = ~.x %in% c(-9999, -8888, -4444, -5555, -6666, -7777, -9, -9.999, -99.9, -0.99))
#categorical variables Education GeneralHealthBaby Dentobuff
#Education
Postnatal$Education_post<-Postnatal$Education
Postnatal$Education_HS_post<-ifelse(Postnatal$Education<4, "High school education or less", "Greater than high school education")
Postnatal$Education2_post<-recode_factor(Postnatal$Education, `1`="High school or less", `2` = "High school or less", `3`="High school or less", `4`="Associates degree or some college", `5`="Associates degree or some college", `6`="Bachelors degree or more", `7`="Bachelors degree or more", `8`="Bachelors degree or more")
Postnatal$Education_post<-recode_factor(Postnatal$Education, `1`="Eight grade or less", `2`='9th-12th grade', `3`="High school graduate or GED", `4` = 'Some college', `5`='Associates degree', `6`="Bachelors degree", `7`="Masters degree", `8`="Doctoral or professional degree")
#General health baby
Postnatal$GeneralHealthBaby<-recode_factor(Postnatal$GeneralHealthBaby, `1`="Excellent", `2`="Good", `3`="Fair", `4`="Poor")
#Dentobuff
Postnatal$Dentobuff<-recode_factor(Postnatal$Dentobuff, `111`="Failed", `1`="Blue", `2`="Green", `3`="Yellow")
Postnatal<-Postnatal %>% select(-Education)
#How many visits does each person have (limit to up to visit 10, as that is the limit of the saliva samples)? 
Postnatal<-Postnatal %>% group_by(BabySubjectID) %>% filter(Visit<=10) %>% add_tally() 


##############Postnatal and prenatal merge
#make metadata dataset with prenatal and postnatal information 
MetaVisit<-left_join(Prenatal, Postnatal)
MetaVisit$COHRAID<-paste(MetaVisit$BabySubjectID, '-', MetaVisit$Visit, sep='')
#how many controls eventually had decay at or before the V10 but were NOT included in this study as cases (hadn't had those visits at the time of ?
MetaVisit %>% filter(CaseStatus=="Control (post-incident visit)" & Decay=="Decay present (d1ft>0)") %>% select(BabySubjectID, IncidentVisit, Visit, Decay, Prim_d1ft) %>% unique()
EverDecay<-MetaVisit %>% filter(Visit<=10 & Prim_d1ft>0)%>%pull(BabySubjectID)%>%unique()
MetaVisit$DecayBy10=ifelse(MetaVisit$BabySubjectID %in% EverDecay, 'Yes', 'No')
#MetaVisit saliva pH
MetaVisit$SalivapH=ifelse(MetaVisit$SalivapH<2.5, NA, MetaVisit$SalivapH)
dir.create(file.path(data.out, 'Metadata'))
#add foxcavid 
load(file.path(data.out, 'phyloseq', 'phyasv_truesamples.Rdata'))
m=as.data.frame(as.matrix(ps_truesamples@sam_data))%>%select(FoxCavID, COHRAID)
MetaVisit_new=left_join(MetaVisit, m)
#check sample size
if(nrow(MetaVisit_new)==nrow(MetaVisit)){MetaVisit=MetaVisit_new}; if(nrow(MetaVisit_new)!=nrow(MetaVisit)){print("STOP: check n")}
#remove one subject who has only one metavisit record & no FoxCavID (i.e. no metagenomics samples) 
StrangeID=MetaVisit%>%group_by(BabySubjectID)%>%dplyr::summarise(n=n())%>%filter(n==1)%>%pull(BabySubjectID) 
MetaVisit=MetaVisit %>%filter(BabySubjectID!=StrangeID)
#chek sample size
if(nrow(MetaVisit_new)-1 != nrow(MetaVisit)){print("STOP:check n")}
#Make a column that denotes if there is a clean 16S saliva sample 
HasCleanSample<-data.frame('COHRAID'=c(ps_truesamples@sam_data$COHRAID))%>%mutate(HasCleanSample=1)
MetaVisit<-left_join(MetaVisit, HasCleanSample)
MetaVisit$UniquePersonID %>% unique() %>% length() #191 
MetaVisit$BabySubjectID %>% unique() %>% length() # 189 (one person a control twice, one person a case & a control)
#Also note that 1 individual was recorded as being a control at Visit 8 but actually has saliva samples out to Visit 10 and never had caries
test<-MetaVisit %>% filter(Visit<=10 & HasCleanSample==T & Visit<=IncidentVisit)%>%filter(!(UniquePersonID %in% c("11000244_Co_8", "22000022_Co_8")))%>%pull(COHRAID)
#this doesn't affect any of the statistical analysis in this paper: for RFs, @ 12 and 24 months, individual is a control no matter what, similarly for chisq/mean comparisons at 12 & 24 months, same  
weirdIDs<-ps_truesamples@sam_data$COHRAID[!(ps_truesamples@sam_data$COHRAID %in% test)]
weirdPerson<-gsub('-.*', '', weirdIDs, )%>%unique()
MetaVisit$IncidentVisit2<-ifelse(MetaVisit$BabySubjectID == weirdPerson, 10, MetaVisit$IncidentVisit)
#save MetaVisit
save(MetaVisit, file=file.path(data.out, 'Metadata', 'MetaVisit.Rdata'))
