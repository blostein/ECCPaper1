library(psych)
library(plyr)
library(dplyr)
library(stringr)
library(chron)
library(lubridate)
library(naniar)
library(compareGroups)
library(readxl)
#function for na recoding 
na_codes <- function(x, ...) {
  x[x %in% c(...)] <- NA
  x
}

#############Read in datasets###############################################################
USCUR<-read.csv(file.path(meta.data, "USCURPostnatal.csv"))
Postnatal<-read.csv(file.path(meta.data, "DFTandPostnatal.csv"))
Prenatal<-read.csv(file.path(meta.data, "Prenatal.csv"))
additional_1=readxl::read_xlsx(file.path(meta.data, 'MaternalFactors_UpdateFromIDLIst_20211005.xlsx'), sheet = 1)
additional_2=readxl::read_xlsx(file.path(meta.data, 'MaternalFactorsCariesDAR_addedTotalTeeth_20210811.xlsx'), sheet = 1)
additionalMeta=rbind(additional_1, additional_2)
additionalMeta_cc<-additionalMeta %>%filter(BabysubjectID%in% USCUR$BabySubjectID)

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

#additional meta data - child race 
babyRaceData=additionalMeta%>%
  mutate(BabySubjectID=BabysubjectID)%>%
  select(BabySubjectID, BabyRace, BabyEthnicity, FatherRace, FatherEthnicity)%>%
  unique()%>%
  mutate(BabyRace=case_when(BabyRace==1~ 'White', 
                            BabyRace==6~ 'Bi- or Multi-racial',
                            BabyRace==-8888~'Unknown'))

##############Postnatal and prenatal merge
#make metadata dataset with prenatal and postnatal information 
MetaVisit<-left_join(left_join(Prenatal, Postnatal), babyRaceData)
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

########################################
#add in antibiotic data
birth_antibiotic = read_xlsx(file.path(meta.data, 'BlosteinDAR_BabyAntibiotics_20220303.xlsx'), sheet='BirthAntibiotic')
visit_antibiotic=read_xlsx(file.path(meta.data, 'BlosteinDAR_BabyAntibiotics_20220303.xlsx'), sheet='PostnatalF2F')
call_antibiotic=read_xlsx(file.path(meta.data, 'BlosteinDAR_BabyAntibiotics_20220303.xlsx'), sheet='ShortPhoneInterview')
MV = MetaVisit %>% select(MotherSubjectID, BabySubjectID, Site, Visit, VisitDate)

visit_antibiotic = left_join(left_join(MV, visit_antibiotic), birth_antibiotic%>%
                               rename_with(~paste0('birth_', .x), BabyAntibiotics:BabyAntibiotic2Days))

visit_antibiotic = visit_antibiotic %>% 
  mutate(antibiotic_within_week=case_when(ABStop %in% c(1, 2, 3) | 
                                            birth_BabyAntibiotics==1 ~ 'Antibiotics reported w/in 1 week of visit',
                                          !(ABStop %in% c(1, 2, 3)) & (birth_BabyAntibiotics !=1 | is.na(birth_BabyAntibiotics)) ~ 'No antibiotics reported w/in 1 week of visit', 
                                          is.na(ABStop) & is.na(birth_BabyAntibiotics) ~ 'No antibiotics reported w/in 1 week of visit'))
all_antibiotic = visit_antibiotic %>% 
  group_by(BabySubjectID)%>%
  arrange(BabySubjectID, VisitDate)%>%
  mutate(last_VisitDate = lag(VisitDate), last_ABStop=lag(ABStop))%>%
  left_join(call_antibiotic %>% 
              rename_with(~paste0('call_', .x), Illness:ABFreq))

all_antibiotic = all_antibiotic%>%
  filter(CallDate<VisitDate & last_VisitDate<CallDate)%>%
  mutate(time_since_call=VisitDate-CallDate, 
         time_since_last_visit=CallDate-last_VisitDate)%>%
  #if there are any calls that happened w/in a week before the visit
  #AND antibiotics were reported AT that call that are either still being taken or were finished w/in 2 days
  #then we can say antibiotics were within a week of that call
  mutate(antibiotic_within_week=case_when(time_since_call<=7 & call_ABStop %in% c(1, 2)~'Antibiotics reported w/in 1 week of visit', 
                                          TRUE~antibiotic_within_week), 
         #let's also flag any such observations
         within_week_flag=case_when(time_since_call<=7 & call_ABStop %in% c(1, 2)~1, 
                                    TRUE~0))%>%
  #otherwise, lets make a call antibiotic variable
  mutate(any_antibiotic_since_last_visit='Antibiotic reported since last visit', 
         time_since_antibiotic = case_when(call_ABStop==1 ~ 0, 
                                           call_ABStop==2~ 2, 
                                           call_ABStop==3~7, 
                                           call_ABStop==4~30,
                                           call_ABStop==5~90, 
                                           call_ABStop==6~120,
                                           TRUE~9999),
         time_since_antibiotic=time_since_call+time_since_antibiotic,
         apprx_days_since_most_recent_antibiotic=case_when(time_since_antibiotic<=90~as.character(time_since_antibiotic), 
                                                           time_since_antibiotic>=9999~'Missing time since antibiotic', 
                                                           time_since_antibiotic>90 & time_since_antibiotic<9999~'>90'), 
         antibiotics_within_3mos_of_visit=ifelse(time_since_antibiotic>90, 'No', 'Yes' ))

recent_call_antibiotic = all_antibiotic %>% 
  group_by(BabySubjectID, Visit)%>%
  filter(time_since_antibiotic==min(time_since_antibiotic))%>%
  select(BabySubjectID, Visit, any_antibiotic_since_last_visit, antibiotics_within_3mos_of_visit, time_since_antibiotic,
         apprx_days_since_most_recent_antibiotic, within_week_flag)%>%
  unique()

recent_visit_antibiotic = visit_antibiotic %>%
  group_by(BabySubjectID, Visit)%>%
  filter(is.na(AB) | ABStop==min(ABStop))%>%
  select(BabySubjectID, Visit, birth_BabyAntibiotics, antibiotic_within_week, AB, ABStop, antibiotic_within_week)%>%
  unique()%>%
  mutate(AB=list(AB))%>%
  unique()

if(
  recent_call_antibiotic%>%unique()%>%group_by(BabySubjectID, Visit)%>%dplyr::summarise(n=n())%>%filter(n>1)%>%nrow() != 0 |
  recent_visit_antibiotic%>%unique()%>%group_by(BabySubjectID, Visit)%>%dplyr::summarise(n=n())%>%filter(n>1)%>%nrow() != 0
){print('>1 record per person-visit for antibiotic data. Check')}

#combine call and visit antibiotic data 
antibiotic_data=left_join(recent_visit_antibiotic, recent_call_antibiotic)
antibiotic_data = antibiotic_data %>% 
  mutate(antibiotic_within_week_wCall=case_when(within_week_flag==1~'Antibiotics reported w/in 1 week of visit', 
                                                within_week_flag==0 | is.na(within_week_flag)~antibiotic_within_week), 
         any_antibiotic_since_last_visit=case_when(!is.na(AB)~'Antibiotic reported since last visit', 
                                                   is.na(AB) & !is.na(any_antibiotic_since_last_visit)~any_antibiotic_since_last_visit, 
                                                   is.na(any_antibiotic_since_last_visit)~'No antibiotic reported since last visit'),
         antibiotics_within_3mos_of_visit= case_when(ABStop %in% c(1, 2, 3, 4, 5)~ 'Yes', 
                                                     !is.na(antibiotics_within_3mos_of_visit) ~ antibiotics_within_3mos_of_visit, 
                                                     is.na(antibiotics_within_3mos_of_visit) & !(ABStop %in% c(1, 2, 3, 4, 5))~'No'), 
         time_since_antibiotic_visit = case_when(ABStop==1 ~ 0, 
                                                 ABStop==2~ 2, 
                                                 ABStop==3~7, 
                                                 ABStop==4~30,
                                                 ABStop==5~90, 
                                                 ABStop==6~120),
         apprx_days_since_most_recent_antibiotic=case_when(is.na(time_since_antibiotic_visit)~apprx_days_since_most_recent_antibiotic, 
                                                           !is.na(time_since_antibiotic_visit) & is.na(apprx_days_since_most_recent_antibiotic)~as.character(time_since_antibiotic_visit), 
                                                           !is.na(time_since_antibiotic_visit) & !is.na(time_since_antibiotic)~
                                                             as.character(min(as.integer(time_since_antibiotic), time_since_antibiotic_visit))))

antibiotic_data_condensed=antibiotic_data%>%
  mutate(antibiotic_within_week=antibiotic_within_week_wCall)%>%
  select(BabySubjectID, Visit, antibiotic_within_week, any_antibiotic_since_last_visit, 
         antibiotics_within_3mos_of_visit, apprx_days_since_most_recent_antibiotic)

antibiotic_data_condensed%>%group_by(BabySubjectID, Visit)%>%dplyr::summarise(n=n())%>%filter(n>1)%>%nrow()==0
#save(antibiotic_data_condensed, file=file.path(datadir, 'metadata', 'antibiotics.Rdata'))

#merge in antibiotic data
MetaVisit=left_join(MetaVisit, antibiotic_data_condensed)

#add an indicator for if in WGS subset 
load(file.path(meta.data, 'metag_metashort.Rdata'))
MetaVisit = MetaVisit %>% 
  mutate(person_in_wgs=ifelse(BabySubjectID %in% metag_metashort$BabySubjectID, 'Child has shotgun sequenced samples from incident visit', 'Child has no shotgun sequenced samples from incident visit'), 
         sample_in_wgs=ifelse(COHRAID %in% metag_metashort$COHRAID, 'Sample in shotgun sequenced subset', 'Sample not in shotgun sequenced subset'))

#add breastfeeding data
load(file.path(meta.data, '2022_04_04Breastfeeding.Rdata'))
MetaVisit=MetaVisit%>%left_join(BF_all)%>%mutate(BFCurrent=ifelse(VisitDate<=BFDate, 'Currently breastfeeding', 'Not currently breastfeeding'))

#Check dimensions and case control counts
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



#Additional metadata comparing whole cohort to nested subset
additionalMeta = additionalMeta %>% arrange(SUBJECTID, PhCall)%>%
  mutate(PrenatalEducation=ifelse(Education<0, NA, Education))%>%
  dplyr::group_by(SUBJECTID)%>%
  tidyr::fill(PrenatalEducation, .direction=c('down'))%>%
  mutate(VisitDate=as.Date(VisitDate, format='%b %d %Y'))

additionalMeta1=additionalMeta %>% 
  mutate(in_ecc=ifelse(BabysubjectID %in% MetaVisit$BabySubjectID, 
                       'In nested case-control sample', 
                       'Not in nested case-control sample'), 
         any_d1mft=ifelse(Prim_d1ft>0, 'Yes', 'No'), 
         in_ecc_case=case_when(BabysubjectID %in% MetaVisit$BabySubjectID &
                              BabysubjectID %in% CasePrime$BabySubjectID~ 
                              'Selected case',
                              BabysubjectID %in% MetaVisit$BabySubjectID &
                                !(BabysubjectID %in% CasePrime$BabySubjectID)~
                                    'Selected control', 
                                  TRUE~'Not selected'))%>%
  dplyr::mutate(across(.cols = c(Prim_Tot_Teeth_Present, Prim_d1ft, Prim_d1fs), 
                ~na_if(na_if(., -9999), -6666)))%>%
  mutate(BabySex=recode(babysex, `1`="Male", `2`="Female", .default=NA_character_),
         BabyRace=recode(BabyRace, `1`='White', `2`='Not white', `6`='Not white', .default=NA_character_),
         Delivery=recode(Delivery, `1`="Vaginal", `2`="C-section", .default=NA_character_), 
         BabyAgeAtExam=ifelse(BabyAgeAtExam>0, BabyAgeAtExam*12, NA), 
         CurrentBreastfed=case_when(Breastfed==2 & BFCurrent==-6666~'Not currently breastfed',
                                    BFCurrent==2~'Not currently breastfed', 
                                    BFCurrent==1~'Currently breastfed', 
                                    BFCurrent==-9999~NA_character_),
        Education_HS=ifelse(PrenatalEducation<0, "Missing", ifelse(PrenatalEducation<5, 
                                                           "High school degree or less", 
                                                           "Associates degree or higher")), 
        Education_Cat=ifelse(PrenatalEducation<0, "Missing", 
                             ifelse(PrenatalEducation<5, 
                                    "High school or less",
                                    ifelse(PrenatalEducation<6, "Associates degree", 
                                           ifelse(PrenatalEducation==6, "Undergraduate degree", 
                                                  ifelse(PrenatalEducation>6, "Graduate degree", NA))))), 
        PrenatalEducation=recode(PrenatalEducation, `1`="8th grade or less", `2`="9th-12th grade, no diploma", 
                                 `3`="High school graduate/GED", `4`="Some college credit, no degree", 
                                 `5`="Associate degree", `6`="Bachelor's degree", `7`="Master's degree", 
                                 `8`="Doctorate/professional degree", .default=NA_character_))%>%
  left_join(BF_all)%>%
  mutate(ogBFCurrent=BFCurrent,
         BFCurrent=ifelse(VisitDate<=BFDate, 'Currently breastfeeding', 'Not currently breastfeeding'))

save(additionalMeta1, file=file.path(data.out, 'Metadata', 'CohortData.Rdata'))


