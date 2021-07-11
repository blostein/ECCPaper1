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
Postnatal %>% filter(n<3) %>% select(BabySubjectID, Visit, Prim_d1ft)
#individual 21000438 has visits 8 and 7 and no decay

##############Postnatal and prenatal merge
#make metadata dataset with prenatal and postnatal information 
MetaVisit<-left_join(Prenatal, Postnatal)
MetaVisit$COHRAID<-paste(MetaVisit$BabySubjectID, '-', MetaVisit$Visit, sep='')
#how many controls eventually had decay at or before the V10 but were NOT included in this study as cases?
MetaVisit %>% filter(CaseStatus=="Control (post-incident visit)" & Decay=="Decay present (d1ft>0)") %>% select(BabySubjectID, IncidentVisit, Visit, Decay, Prim_d1ft) %>% unique()
#   BabySubjectID IncidentVisit Visit   Decay                      Prim_d1ft
#2      11000326             8    10    Decay present (d1ft>0)         4
#3      11000320             8    10    Decay present (d1ft>0)         4
#4      11000444             8    10    Decay present (d1ft>0)         1
#5      21000386             8     9    Decay present (d1ft>0)         8
EverDecay<-MetaVisit %>% filter(Visit<=10 & Prim_d1ft>0)%>%pull(BabySubjectID)%>%unique()
MetaVisit$DecayBy10=ifelse(MetaVisit$BabySubjectID %in% EverDecay, 'Yes', 'No')
#MetaVisit saliva pH
#ggplot(MetaVisit, aes(x=1, y=SalivapH))+geom_boxplot()
MetaVisit$SalivapH=ifelse(MetaVisit$SalivapH<2.5, NA, MetaVisit$SalivapH)
dir.create(file.path(data.out, 'Metadata'))
#add foxcavid 
load(file.path(data.out, 'phyloseq', 'phyasv_truesamples.Rdata'))
m=as.data.frame(as.matrix(ps_truesamples@sam_data))%>%select(FoxCavID, COHRAID)
MetaVisit_new=left_join(MetaVisit, m)
if(nrow(MetaVisit_new)==nrow(MetaVisit)){MetaVisit=MetaVisit_new}; if(nrow(MetaVisit_new)!=nrow(MetaVisit)){print("STOP: check n")}
#remove one subject who has only one metavisit record & no FoxCavID (i.e. no metagenomics samples) 
StrangeID=MetaVisit%>%group_by(BabySubjectID)%>%dplyr::summarise(n=n())%>%filter(n==1)%>%pull(BabySubjectID) 
MetaVisit=MetaVisit %>%filter(BabySubjectID!=StrangeID)
if(nrow(MetaVisit_new)-1 != nrow(MetaVisit)){print("STOP:check n")}
#Make a column that denotes if there is a saliva sample 
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


#####################Phone call interview data
#USCUR 
USCUR<-USCUR %>% dplyr::mutate(PhoneCall=factor(recode(as.numeric(PhCall), '1' = '10 wk call', '2'= '6 mo call', '3'='12 mo call', '4'='18 mo call', 
                                                '5'='24 mo call', '6'= '30 mo call', '7'='36 mo call', '8'='42 mo call', '9'='48 mo call', 
                                                '10'='54 mo call', '11'='60 mo call', '12'='66 mo call', '13'= '72 mo call'), 
                                         levels=c("10 wk call", "6 mo call", "12 mo call", "18 mo call", "24 mo call", "30 mo call", "36 mo call",
                                                  "42 mo call", "48 mo call", "54 mo call", "60 mo call", "66 mo call", "72 mo call")))%>%
  dplyr::mutate(CallDate=as.Date(CallDate, format="%m/%d/%y"))

# 1-4 houskeeping
#5-9 bf and bm 
#10-49 diet
#50-58 tooth hygiene 

#Date CallDate

#BF variables Q35 
BF<-USCUR %>% select(MotherSubjectID, BabySubjectID, PhCall, CallDate, Q35_Breastfed, Q35a_BFCurrent, Q35b_BFStopMnth, Q35b_BFStopWks, Q37a_BMCurrent)
BF %>% select(BabySubjectID) %>% unique() %>% dim() # 189 total ids 
BF<-left_join(BF, Prenatal %>% select(MotherSubjectID, BabySubjectID, BabyBirthdate))
BF$CallDate<-as.Date(BF$CallDate, format="%m/%d/%y")
class(BF$CallDate)
BF %>% filter(is.na(CallDate)) %>% dim() # 2 missing call date 
BF %>% filter(is.na(CallDate))
BF<-BF %>% filter(!is.na(CallDate)) # remove
table(BF$Q35b_BFStopMnth)
BF$Q35b_BFStopMnth<-str_replace_all(BF$Q35b_BFStopMnth, c(" months" = "", " month"="", " months"= "", " MONTHS"=""))
table(BF$Q35b_BFStopMnth)
table(BF$Q35b_BFStopWks)
BF<-BF %>% mutate_at(vars(contains("Q35b_BFStop")), funs(na_codes(.,-9999, 88)))%>% mutate_at(vars(contains("Q35b_BFStop")), funs(as.numeric(.)))
BF_StopDate<-BF %>% filter((!is.na(Q35b_BFStopMnth)) | !is.na(Q35b_BFStopWks))
BF_StopDate %>% dim()
BF_StopDate %>% select(BabySubjectID) %>% unique() %>% dim() # no repeated ids, ever
#individual with a stop date given has only one stop date given 
a<- BF %>% filter(!(BabySubjectID %in% BF_StopDate$BabySubjectID))
No <- a %>% filter(Q35_Breastfed==2)
Yes <- a  %>% filter(Q35_Breastfed==1)
Mixed1<-unique(No[No$BabySubjectID %in% Yes$BabySubjectID,1])
Mixed2<-unique(Yes[Yes$BabySubjectID %in% No$BabySubjectID,1])
Mixed1==Mixed2 # good
Never <- unique(No[!(No$BabySubjectID %in% Yes$BabySubjectID), 1])
Always<- unique(Yes[!(Yes$BabySubjectID %in% No$BabySubjectID), 1])
length(Never)+length(Always)+length(Mixed1)==dim(a %>% select(BabySubjectID) %>% unique())[1]
length(Never)+length(Always)+length(Mixed1)+(BF_StopDate %>% select(BabySubjectID) %>% unique() %>% dim())[1] ==189
#make new variables Type and BFDate
BF_StopDate<-left_join(BF_StopDate %>% dplyr::mutate(Type="BF stop date known", BFDate = BabyBirthdate %m+% months(Q35b_BFStopMnth)+ ifelse(is.na(Q35b_BFStopWks), 0, weeks(Q35b_BFStopWks))) %>% select(BabySubjectID, Type, BFDate), BF)
BF_Never<-BF %>% filter(MotherSubjectID %in% Never) %>% dplyr::mutate(Type="BF never", BFDate=BabyBirthdate)
LastYesDate<-BF %>% filter(MotherSubjectID %in% Mixed1 & Q35_Breastfed==1) %>%
  group_by(BabySubjectID) %>%
  top_n(n=1, wt=CallDate)
FirstNoDate<-BF %>% filter(MotherSubjectID %in% Mixed1 & Q35_Breastfed==2) %>% 
  group_by(BabySubjectID) %>% 
  top_n(n=1, wt=rev(CallDate))
dim(LastYesDate)
dim(FirstNoDate)
FirstNoDate$CallDate>LastYesDate$CallDate
WrongMixed<-FirstNoDate[(FirstNoDate$CallDate>LastYesDate$CallDate)==F, ]
BF %>% filter(BabySubjectID == WrongMixed$BabySubjectID) # all 1s except the one that is a 2 
# move to always 
LastAlways<- BF %>% filter(MotherSubjectID %in% c(Always, WrongMixed$MotherSubjectID) & Q35a_BFCurrent==1)%>%
  group_by(BabySubjectID) %>% 
  top_n(n=1, wt=CallDate)%>% dplyr::mutate(Type="BF always", BFDate=CallDate) %>% select(BabySubjectID, Type, BFDate)
BF_Always<-left_join(BF %>% filter(MotherSubjectID %in% c(Always, WrongMixed$MotherSubjectID)), LastAlways)
# make mixed and remove 
MixedDate=LastYesDate$CallDate + floor((FirstNoDate$CallDate-LastYesDate$CallDate)/2)
BF_Mixed<-left_join(left_join(FirstNoDate %>% dplyr::mutate(Type="BF stop date unknown") %>% 
                                filter(BabySubjectID != WrongMixed$BabySubjectID) %>% select(BabySubjectID, CallDate, Type), 
                              LastYesDate %>% dplyr::mutate(LastYesDate=CallDate) %>% select(BabySubjectID, LastYesDate))%>%
                      dplyr::mutate(BFDate=LastYesDate+floor((CallDate-LastYesDate)/2)) %>% select(BabySubjectID, BFDate, Type), BF)

(BF_Mixed %>% select(BabySubjectID) %>% unique()%>%dim())[1]==length(Mixed1[Mixed1!=WrongMixed$MotherSubjectID])
#rejoin
BF2<-rbind(BF_Always, BF_Never, BF_StopDate, BF_Mixed)
BF2$BFDuration<-BF2$BFDate-BF2$BabyBirthdate
BF2$BFDurationCat<-cut(as.numeric(BF2$BFDate-BF2$BabyBirthdate), c(0, 30, 180, 365, 5000), include.lowest = T, labels=c("Never or for <=1 month", ">1 month & <=6 months", ">6 months & <=12 months", ">12 months"))
with(BF2, table(Type, BFDurationCat))
BF3<-BF2 %>% select(MotherSubjectID, BabySubjectID, Type, BFDate, BFDuration, BFDurationCat)%>%unique()
with(BF3, table(Type, BFDurationCat))


#oral hygiene 
OH<-USCUR %>% select(MotherSubjectID, BabySubjectID, PhCall, CallDate, Q48a_WipeTeeth, Q48a_WipeTeeth1, Q48a_WipeTeeth2, 
                     Q48a_WipeTeeth3, Q48a_WipeTeeth4, Q48b_WipeFreq, Q48f_BrushSelf, Q48f_BrushSelfFreq, Q48c_Toothpaste, Q48d_Fluoride)
#figuring out the weirdness of Q48a_WipeTeeth1-4:
# in DD, 1 is "Yes brush", 2 is "No", 3 is "DK", 4 is "Yes Wipe" 
#for all, 1 is yes and 0 is no 
with(OH, table(Q48a_WipeTeeth1, Q48a_WipeTeeth2))
# no one who says 'yes' to 1 also says yes to 2
with(OH, table(Q48a_WipeTeeth1, Q48a_WipeTeeth3))
#2 people who say yes to 1 say yes to 3
with(OH, table(Q48a_WipeTeeth1, Q48a_WipeTeeth4))
#6 people who say yes to 1 also say yes to 4
#conclusion: just use 1... need the b variable (frequency) 
table(OH$Q48a_WipeTeeth, OH$Q48a_WipeTeeth1)# for 48a 1 is yes,2  is no, and 9 or -9999 is missing, for 48a1, 0 is no, 1 is yes
#use either yes to Q48a or yes to 48a1 
OH$BrushTeeth<-with(OH, ifelse(Q48a_WipeTeeth==1 | Q48a_WipeTeeth1==1, "Yes", ifelse(Q48a_WipeTeeth==2 | Q48a_WipeTeeth1==0, "No", "Missing")))
OH<-OH %>% dplyr::mutate(BrushTeethFreq=factor(case_when(BrushTeeth=="Yes" & Q48b_WipeFreq<5 & Q48b_WipeFreq>0 ~ "Less than once a day", 
                                                  BrushTeeth=="Yes" & Q48b_WipeFreq==5~ "Once a day", 
                                                  BrushTeeth=="Yes" & Q48b_WipeFreq %in% c(6, 7)~ "More than once a day", 
                                                  BrushTeeth=="Yes" & Q48b_WipeFreq==0~ "Yes, frequency not given", 
                                                  BrushTeeth=="Yes" & Q48b_WipeFreq==9~ "Yes, frequency not given",
                                                  BrushTeeth=="Yes" & Q48b_WipeFreq==-9999~ "Yes, frequency not given",
                                                  BrushTeeth=="No" & Q48b_WipeFreq==-9999 ~"No", 
                                                  BrushTeeth=="Missing"~"Missing"), levels=c("Missing", "No", "Less than once a day", "Once a day", "More than once a day", "Yes, frequency not given")))
with(OH, table(BrushTeeth, BrushTeethFreq))
with(OH, table(Q48f_BrushSelf, Q48f_BrushSelfFreq))
OH<-OH %>% dplyr::mutate(BrushSelf = factor(case_when(Q48f_BrushSelf==-9999 | Q48f_BrushSelf==9 ~ "Missing", 
                                               Q48f_BrushSelf==2 ~ "No", 
                                               Q48f_BrushSelf==1 & Q48f_BrushSelfFreq<5 ~ "Less than once a day", 
                                               Q48f_BrushSelf==1 & Q48f_BrushSelfFreq==5 ~ "Once a day", 
                                               Q48f_BrushSelf==1 & Q48f_BrushSelfFreq>5 ~ "More than once a day"), 
                                     levels=c("Missing", "No", "Less than once a day", "Once a day", "More than once a day")))%>%
  dplyr::mutate(Toothpaste=factor(case_when(Q48c_Toothpaste==-9999| Q48c_Toothpaste==9 ~ "Missing", 
                                     Q48c_Toothpaste==2 ~ "No", 
                                     Q48c_Toothpaste==1 & Q48d_Fluoride==1~ "Yes, w/ fluoride", 
                                     Q48c_Toothpaste==1 & Q48d_Fluoride==2~ "Yes, w/o fluoride", 
                                     Q48c_Toothpaste==1 & Q48d_Fluoride==9~ "Yes, unknown fluoride status"), levels = c("Missing", "No", "Yes, w/ fluoride", "Yes, w/o fluoride", "Yes, unknown fluoride status")))
with(OH, table(BrushSelf, Q48c_Toothpaste)) # 1 is yes, 2 is no 9 and -9999 is missing
with(OH, table(BrushTeeth, Q48c_Toothpaste))
with(OH, table(BrushTeeth, BrushSelf, Q48c_Toothpaste))
with(OH, table(Q48c_Toothpaste, Q48d_Fluoride))
with(OH, table(BrushTeeth, BrushSelf, Toothpaste))

with(OH, table(BrushTeethFreq, BrushSelf))
OH<-OH %>% dplyr::mutate(TeethBrushed=factor(case_when(BrushSelf=="More than once a day" | BrushTeethFreq=="More than once a day"~"More than once a day",
                                                BrushSelf!="More than once a day" & BrushTeethFreq!="More than once a day" & (BrushSelf=="Once a day"|BrushTeethFreq=="Once a day")~"Once a day", 
                                                BrushSelf!="More than once a day" & BrushTeethFreq!="More than once a day" & BrushSelf!="Once a day" & BrushTeethFreq!="Once a day"& (BrushSelf=="Less than once a day"| BrushTeethFreq=="Less than once a day")~"Less than once a day", 
                                                BrushTeethFreq=="No" & BrushSelf%in%c("Missing", "No")~"No", 
                                                BrushTeethFreq=="Missing" & BrushSelf=="Missing" ~"Missing"), 
                                      levels=c("Missing", "No", "Less than once a day", "Once a day", "More than once a day")))

#Diet 
Diet<-USCUR[,c(1:5, 11:50, 61)]
Diet<-Diet %>% mutate_at(vars(matches("Q38|Q43|Q42")), funs(as.factor(.)))
table(Diet$PhoneCall)
Diet %>% group_by(BabySubjectID) %>% dplyr::summarise(n=n()) %>% select(n) %>% table()
#Diet %>% group_by(BabySubjectID) %>% dplyr::summarise(n=n()) %>% select(n)%>%hist()

BabyFood<-Diet %>% select(BabySubjectID, Q42_BabyFood, PhoneCall, CallDate) %>% arrange(BabySubjectID, CallDate) %>% group_by(BabySubjectID) %>% 
  filter(Q42_BabyFood %in% c(1, 4)) %>% dplyr::slice(1)%>%dplyr::mutate(FirstBabyFood=case_when(PhoneCall=="10 wk call"~"At 10 wk call", 
                                                                                  PhoneCall=="6 mo call"~"At 6 mo call", 
                                                                                  PhoneCall%in% c("12 mo call", "18 mo call", "24 mo call")~"At or after 12 mo call"))%>%
  select(BabySubjectID, FirstBabyFood)
with(BabyFood, table(FirstBabyFood))

Diet<-Diet %>% mutate_at(vars("Q38b_LiquidFormula", 'Q38c_PowderFormula', 'Q38d_CowAnimalMilk', 'Q38i_Juice', 'Q38l_Soda', 'Q43j_EatDesserts'), funs(new=factor(case_when(. %in% c(9, -9999, 1)~"Never or once", 
                                                                                                                                                                          .==2~"Every few days", .==3~"Once per day", .==4~"Several times per day"), 
                                                                                                                                                                levels=c("Never or once", "Every few days", "Once per day", "Several times per day"))))
#SSB
#variables to use "Q38g_FlavWater" "Q38g_FlavWaterSweet"    "Q38h_SportsDrink"      
# "Q38h_SportsSweet"       "Q38i_Juice" "Q38ii_JuiceWat"         "Q38j_JuiceDrink"  
#"Q38j_JuiceDrinksweet"   "Q38k_PowderMix" "Q38k_PowderMixSweet"    "Q38l_Soda"            
#"Q38l_SodaSweet"         "Q38m_Coffee" "Q38m_CoffeeSweet"       "Q38n_Tea"               
#"Q38n_TeaSweet"          "Q38o_MealDrink" "Q38o_MealDrinkSweet"    "Q38p_EnergyDrink" 
#"Q38p_EnergyDrinkSweet"
with(Diet, xtabs(~Q38g_FlavWater+Q38g_FlavWaterSweet))
with(Diet, xtabs(~Q38n_Tea+Q38n_TeaSweet))
with(Diet, xtabs(~Q38l_Soda+Q38l_SodaSweet))
#Take beveragesthat can be sugar sweetened. If they are sugar sweetened, sum across all SSB weekly frequencies
#if not SS (or missing SS info), treat frequency as 0. then z score standardize 
#SSBScore includes only beverages with SS info. 
#SSBScore2 includes these beverages and juice and watered juice frequency. 

SSB<-Diet[, c(1, 2, 3, 4, 12:15, 18:31)]
var_freq<-colnames(SSB)[seq(5, 22, by=2)]
var_sweet<-colnames(SSB)[seq(6, 22, by=2)]
var_new<-paste(var_freq, '_new', sep='')
for(i in seq(1:9)){
  SSB<-SSB %>% dplyr::mutate(!!sym(var_new[i]):=case_when(!!sym(var_sweet[i])==1~as.numeric(as.character(!!sym(var_freq[i]))), 
                                                          !!sym(var_sweet[i])!=1~0)) 
}
#check
with(SSB, xtabs(~Q38l_Soda+Q38l_Soda_new+Q38l_SodaSweet))
with(SSB, xtabs(~Q38n_Tea+Q38n_Tea_new+Q38n_TeaSweet))
with(SSB, xtabs(~Q38g_FlavWater+Q38g_FlavWater_new+Q38g_FlavWaterSweet))
#make z score standardizations of summed variables. 
SSB$SSBScore<-as.vector(scale(SSB %>% select(contains("new"))%>%rowSums(.)))
SSB$SSBScore2<-as.vector(scale(left_join(SSB, Diet %>% select(BabySubjectID, PhCall, Q38i_Juice, Q38ii_JuiceWat)%>%mutate_at(vars(contains("Juice")), funs(new=case_when(.%in%c(-9999, 1,9)~0,!(.%in%c(-9999,1,9))~as.numeric(as.character(.))))))%>%
                                 select(contains("new"))%>%rowSums(.)))
#check correlation
#ggplot(SSB, aes(x=SSBScore, y=SSBScore2, color=as.factor(Diet$Q38i_Juice)))+geom_point()  

#join together
Diet2<-left_join(left_join(Diet, SSB %>% select(BabySubjectID, PhCall, SSBScore, SSBScore2)), MetaVisit %>% select(BabySubjectID, Case)) 
CallData<-left_join(left_join(BF3, OH), left_join(left_join(Diet, BabyFood), SSB%>% select(BabySubjectID, PhCall, SSBScore, SSBScore2)))
save(CallData, file=file.path(data.out, 'Metadata', 'CallDataUSCUR.Rdata'))

