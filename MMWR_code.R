library(haven)
library(ggplot2)
library(epiR)
library(tidyr)
library(corrplot)
library(Hmisc)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load("baseline4.rdata")
# load("daily4.rdata")

baseline=read_sas("baseline.sas7bdat")
daily=read_sas("daily.sas7bdat")

#calculating saliva dichotomous result (at least 2 genes amplify)
NandO=NandS=SandO=rep(0,nrow(daily))
NandO[which(daily$saliva_ct_n>0&daily$saliva_ct_orf1>0)]=1
NandS[which(daily$saliva_ct_n>0&daily$saliva_ct_s>0)]=1
SandO[which(daily$saliva_ct_s>0&daily$saliva_ct_orf1>0)]=1
daily$saliva_result=0;daily$saliva_result[which(NandO==1|NandS==1|SandO==1)]=1

#taking average of saliva Ct
daily$saliva_ave=apply(daily[,c(21,22,24)],1,mean)
daily$saliva_ave[which(NandO==1&daily$saliva_ct_s==0)]=apply(daily[which(NandO==1&daily$saliva_ct_s==0),c(22,24)],1,mean)
daily$saliva_ave[which(SandO==1&daily$saliva_ct_n==0)]=apply(daily[which(SandO==1&daily$saliva_ct_n==0),c(22,24)],1,mean)
daily$saliva_ave[which(NandS==1&daily$saliva_ct_orf1==0)]=apply(daily[which(NandS==1&daily$saliva_ct_orf1==0),c(22,24)],1,mean)

daily$saliva_ave[which(daily$saliva_ave==0)]=NA

#Clean results
daily$nasal_swab_result[which(daily$nasal_swab_result==77)]=NA

daily$nasal_swab_antigen=factor(daily$nasal_swab_antigen,levels = c(1,0))
daily$nasal_swab_result=factor(daily$nasal_swab_result,levels = c(1,0))
daily$saliva_result=factor(daily$saliva_result,levels = c(1,0))

#######combine variables
baseline$gender=NA
baseline$gender[which(baseline$gender_man==1)]="Male"
baseline$gender[which(baseline$gender_woman==1)]="Female"
baseline$gender[which(baseline$gender_non_bin==1)]="NonBinary"
baseline$gender[which(baseline$gender_other==1)]="Other"
baseline$gender=factor(baseline$gender)
baseline$race=NA
baseline$race[which(baseline$race_amind==1)]="AmInd"
baseline$race[which(baseline$race_white==1)]="White"
baseline$race[which(baseline$race_black==1)]="Black"
baseline$race[which(baseline$race_asian==1)]="Asian"
baseline$race[which(baseline$race_pacisl==1)]="PacIs"
baseline$race[which(baseline$race_other==1)]="Other"
baseline$race=factor(baseline$race)
baseline$ethnicity=NA
baseline$ethnicity[which(baseline$hispanic_no==1)]="NonHispanic"
baseline$ethnicity[which(baseline$hispanic_mex==1|baseline$hispanic_pr==1|baseline$hispanic_cub==1|baseline$hispanic_other==1)]="Hispanic"
baseline$ethnicity=factor(baseline$ethnicity)

###combine datasets
dsmatch=match(daily$user_id,baseline$user_id)
democol=c(2:17,35:38,40:44)
alldata=cbind(daily,baseline[dsmatch,democol])
alldata$daily_dt=as.Date(alldata$symptoms_dt)
alldata$base_date=as.Date(alldata$base_date)
alldata$study_day=alldata$daily_dt-alldata$base_date

#at any time point, is any sample positive?
for(i in c(5,7:18,43:56,62:66)){alldata[,i]=factor(alldata[,i])}
alldata$parallel_result=apply(cbind(alldata$nasal_swab_result,alldata$saliva_result,alldata$nasal_swab_antigen,alldata$nasal_swab_virus),1,max,na.rm=T)-1
alldata$parallel_result[which(alldata$parallel_result==-Inf)]=NA
alldata$parallel_result=factor(alldata$parallel_result)

#combining saliva, antigen, and nasal results
alldata$combined=NA
alldata$combined[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==1)]="all"
alldata$combined[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==1)]="nasal+antigen"
alldata$combined[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==0)]="saliva+antigen"
alldata$combined[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==1)]="nasal+saliva"
alldata$combined[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==1)]="nasal"
alldata$combined[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==0)]="antigen"
alldata$combined[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==0)]="saliva"
alldata$combined[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==0)]="none"
alldata$combined=factor(alldata$combined)

#combining saliva and antigen results
alldata$combo2=NA
alldata$combo2[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==1)]="saliva+antigen"
alldata$combo2[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1)]="antigen"
alldata$combo2[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0)]="saliva"
alldata$combo2[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0)]="none"
alldata$combo2=factor(alldata$combo2)

#combining culture, saliva, antigen, and nasal results
alldata$combo3=NA
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==1)]="all"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==1)]="nasal+antigen+culture"
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==1)]="saliva+antigen+culture"
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==1)]="nasal+saliva+culture"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==1)]="nasal+culture"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==1)]="antigen+culture"
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==1)]="saliva+culture"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==1)]="culture"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==0)]="nasal+antigen"
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==0)]="saliva+antigen"
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==0)]="nasal+saliva"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==1&alldata$nasal_swab_virus==0)]="nasal"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==1&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==0)]="antigen"
alldata$combo3[which(alldata$saliva_result==1&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==0)]="saliva"
alldata$combo3[which(alldata$saliva_result==0&alldata$nasal_swab_antigen==0&alldata$nasal_swab_result==0&alldata$nasal_swab_virus==0)]="none"
alldata$combo3=factor(alldata$combo3)

#calculating times for alignment
uid=unique(alldata$user_id)
fullresult=firstantigen=max_saliva=max_nasal=viral=nonviral=symptoms=rep(NA,length(uid))
for(i in 1:length(uid)){
  pars=alldata$parallel_result[which(alldata$user_id==uid[i])]
  ag=alldata$nasal_swab_antigen[which(alldata$user_id==uid[i])]
  day=alldata$counter[which(alldata$user_id==uid[i])]
  if(length(day)>1){
    #ever positive
    fullresult[i]=ifelse(max(as.numeric(pars),na.rm=T)==2,1,0)
    #antigen
    firstantigen[i]=day[min(which(ag==1))]
    #saliva peak
    salivas=alldata$saliva_ave[which(alldata$user_id==uid[i])]
    if(min(salivas,na.rm=T)!=Inf){max_saliva[i]=day[which(salivas==min(salivas,na.rm=T))]}
    #nasal peak
    nasals=alldata$nasal_swab_ct_1[which(alldata$user_id==uid[i])]
    if(min(nasals,na.rm = T)!=Inf){
      max_nasal[i]=as.numeric(day[which(nasals==min(nasals,na.rm=T))])}
    #symptom start
    nosymptoms=alldata$no_symptoms[which(alldata$user_id==uid[i])]
    if(is.na(nosymptoms[1])){symptoms[i]=day[min(which(is.na(nosymptoms[2:length(nosymptoms)])))]
    }else{symptoms[i]=day[1]}
    #viral culture start
    culture=alldata$nasal_swab_virus[which(alldata$user_id==uid[i])]
    viral[i]=day[min(which(culture==1))]
    nonviral[i]=day[max(which(culture==1))]
  }}
alldata$userpos=fullresult[match(alldata$user_id,uid)]
alldata$firstantigenpos=firstantigen[match(alldata$user_id,uid)]
alldata$max_saliva=max_saliva[match(alldata$user_id,uid)]
alldata$max_nasal=max_nasal[match(alldata$user_id,uid)]
alldata$culture_start=viral[match(alldata$user_id,uid)]
alldata$culture_end=nonviral[match(alldata$user_id,uid)]
alldata$days_since_first_antigen_pos=as.numeric(alldata$counter-alldata$firstantigenpos)
alldata$days_since_peak_saliva=as.numeric(alldata$counter)-alldata$max_saliva
alldata$days_since_peak_nasal=as.numeric(alldata$counter)-alldata$max_nasal
alldata$days_since_culture_start=as.numeric(alldata$counter)-alldata$culture_start
alldata$days_since_culture_end=alldata$culture_end-as.numeric(alldata$counter)
alldata$days_since_symptom_start[which(is.na(alldata$days_since_symptom_start))]=alldata$counter[which(is.na(alldata$days_since_symptom_start))]-symptoms[match(alldata$user_id[which(is.na(alldata$days_since_symptom_start))],uid)]

alldata$symptoms_vs_saliva=alldata$days_since_peak_saliva-alldata$days_since_symptom_start
alldata$symptoms_vs_nasal=alldata$days_since_peak_nasal-alldata$days_since_symptom_start
alldata$symptoms_vs_antigen=alldata$days_since_first_antigen_pos-alldata$days_since_symptom_start
alldata$presymptomatic=0;alldata$presymptomatic[which(alldata$days_since_symptom_start<0)]=1
alldata$presymptomatic=factor(alldata$presymptomatic)

#count symptoms
symptoms=alldata[,7:17]
for(i in 1:ncol(symptoms)){symptoms[,i]=as.numeric(symptoms[,i])}
alldata$symptom_count=apply(symptoms,1,sum,na.rm=T)

uid=unique(alldata$user_id)
alldata$increment=match(alldata$user_id,uid)
alldata$age=as.numeric(as.character(alldata$age))

save(alldata,file="CombinedCOVIDdetectData.Rdata")
write.csv(alldata,file="cleaned_data.csv")


#######Dataset for VOI
alldata=alldata[-which(alldata$user_id==435772&alldata$counter==3),]
alldata=alldata[-which(alldata$user_id==460545&alldata$counter==15),] #15?

alldata$user_id=factor(alldata$user_id)
wdata=pivot_wider(alldata[which(is.na(alldata$study_day)==F),],id_cols = user_id,names_from = study_day,values_from = c(saliva_result,nasal_swab_antigen,nasal_swab_result))
wdata=as.data.frame(wdata)
wdata$user_result=apply(wdata[,c(3:16,18:31,33:46)],1,max,na.rm=T)
wdata$user_result=factor(wdata$user_result,levels = c(1,0))
pdata=wdata[which(wdata$user_result==1),] #select only positives, using nasal PCR
save(pdata,file="positive_wide_data.Rdata")
save(wdata,file="wide_data.Rdata")

infdata=alldata
infdata$nasal_swab_result[which(alldata$days_since_culture_end<=0)]=infdata$saliva_result[which(alldata$days_since_culture_end<=0)]=infdata$nasal_swab_antigen[which(alldata$days_since_culture_end<=0)]=NA
winfdata=pivot_wider(infdata[which(is.na(infdata$study_day)==F),],id_cols = user_id,names_from = study_day,values_from = c(saliva_result,nasal_swab_antigen,nasal_swab_result))
winfdata=as.data.frame(winfdata)
pwinfdata=winfdata[which(wdata$user_result==1),]
save(pwinfdata,file="positive_wide_infectious_data.Rdata")

uin=unique(alldata$user_id)
puin=uin[2:31]
predata=alldata[which(alldata$user_id%in%puin),]
save(predata,file="PreprintData.Rdata")


#Table 1
demograph=predata[match(puin,alldata$user_id),c(58:67)]
demograph$height_in_inches=demograph$height_ft*12+demograph$height_in
library(tableone)
t1=CreateTableOne(data = demograph[,c(1,11,5,7,9,10)])
write.csv(print(t1,showAllLevels = T),file="MMWR_Table1.csv")
table(demograph$race,demograph$ethnicity)
table(demograph$race,demograph$ethnicity)/nrow(demograph)

#Figure 1 and Table 2
spdata=predata[which(predata$user_id%in%pdata$user_id),] #select only nasal positive ever
spdata=spdata[-which(spdata$counter==0),]
#### Sensitivity of test and Ct by days since culture start
sdays=seq(from=min(spdata$days_since_culture_start,na.rm = T),to=max(spdata$days_since_culture_start,na.rm = T),by=1)
sympt_sens=matrix(NA,nrow=length(sdays),ncol=19);rownames(sympt_sens)=sdays
colnames(sympt_sens)=c(paste(rep(c("Ag Se","Saliva Se","Nasal Se","Nasal Ct"),each=3),c(rep(c(""," LCI"," UCI"),3),c(" Mean"," Min"," Max")),sep=""),
                       paste(rep(c("Ag","Saliva","Nasal"),each=2),rep(c("Npos","N"),3)),"N_culture")

for(d in 1:length(sdays)){
  dats=spdata[which(spdata$days_since_culture_start==sdays[d]),]
  if(nrow(dats)>5){
    sympt_sens[d,1:3]=binconf(length(which(dats$nasal_swab_antigen==1)),length(which(dats$counter>0))) 
    sympt_sens[d,4:6]=binconf(length(which(dats$saliva_result==1)),nrow(dats)) 
    sympt_sens[d,7:9]=binconf(length(which(dats$nasal_swab_result==1)),length(which(dats$counter>0))) 
    sympt_sens[d,10]=mean(dats$nasal_swab_ct_1,na.rm = T)
    sympt_sens[d,11:12]=range(dats$nasal_swab_ct_1,na.rm = T)
    sympt_sens[d,13]=length(which(dats$nasal_swab_antigen==1))
    sympt_sens[d,15]=length(which(dats$saliva_result==1))
    sympt_sens[d,17]=length(which(dats$nasal_swab_result==1))
    sympt_sens[d,14]=length(which(dats$counter>0))
    sympt_sens[d,16]=nrow(dats)
    sympt_sens[d,18]=length(which(dats$counter>0))
    sympt_sens[d,19]=length(which(dats$nasal_swab_virus==1))
  }}
sympt_sens=as.data.frame(sympt_sens)
sympt_sens$day=sdays
names(sympt_sens)[10:12]=c("NasalCt","CtMin","CtMax")
sympt_sens$CtMin[which(sympt_sens$CtMin==Inf)]=sympt_sens$CtMax[which(sympt_sens$CtMax==-Inf)]=NA
sympt_sens[,1:12]=round(sympt_sens[,1:12],3)
write.csv(sympt_sens,file="MMWR_Table2.csv")
sy_l=pivot_longer(sympt_sens,c(1,4,7),names_to = "Test",values_to = "Sensitivity")
sy_lL=pivot_longer(sympt_sens,c(2,5,8),names_to = "LCI")
sy_lU=pivot_longer(sympt_sens,c(3,6,9),names_to = "UCI")
sy_l$LCI=sy_lL$value
sy_l$UCI=sy_lU$value

sy_l$Test=factor(sy_l$Test);levels(sy_l$Test)=c("Antigen","Nasal RTqPCR","Saliva RTqPCR")

tiff("MMWR_Figure1.tif",width = 800,height = 900,res=100)
ggplot(sy_l[-which(is.na(sy_l$Sensitivity)),],aes(x=day,y=Sensitivity,ymin=LCI,ymax=UCI))+
  geom_ribbon(alpha=0.1)+geom_line()+facet_grid(Test~.)+
  xlab("Days Since First Viral Culture")+theme(text=element_text(size=25),strip.text.y = element_text(size = 25))+
  ylab("Daily Sensitivity")+scale_y_continuous(breaks = seq(0,1,.2),limits = c(0,1))
dev.off()

###Figure 2
pvdata=alldata[which(alldata$days_since_culture_start<0&alldata$user_id%in%pdata$user_id&alldata$counter>0),]
pvs=binom.test(length(which(pvdata$saliva_result==1)),nrow(pvdata))
pva=binom.test(length(which(pvdata$nasal_swab_antigen==1)),nrow(pvdata))
pvn=binom.test(length(which(pvdata$nasal_swab_result==1)),nrow(pvdata))
vdata=alldata[which(alldata$nasal_swab_virus==1&alldata$user_id%in%pdata$user_id&is.na(alldata$nasal_swab_virus)==F),]
vs=binom.test(length(which(vdata$saliva_result==1)),nrow(vdata))
va=binom.test(length(which(vdata$nasal_swab_antigen==1)),nrow(vdata))
vn=binom.test(length(which(vdata$nasal_swab_result==1)),nrow(vdata))
nvdata=alldata[which(alldata$nasal_swab_virus==0&alldata$user_id%in%pdata$user_id&alldata$days_since_culture_start<=7),]
nvdata=nvdata[which(nvdata$days_since_culture_start>0|is.na(nvdata$days_since_culture_start)),]
nvs=binom.test(length(which(nvdata$saliva_result==1)),nrow(nvdata))
nva=binom.test(length(which(nvdata$nasal_swab_antigen==1)),nrow(nvdata))
nvn=binom.test(length(which(nvdata$nasal_swab_result==1)),nrow(nvdata))
shedse=data.frame(Sensitivity=c(pvs$estimate,pva$estimate,pvn$estimate,vs$estimate,va$estimate,vn$estimate,nvs$estimate,nva$estimate,nvn$estimate),
                  LCI=c(pvs$conf.int[1],pva$conf.int[1],pvn$conf.int[1],vs$conf.int[1],va$conf.int[1],vn$conf.int[1],nvs$conf.int[1],nva$conf.int[1],nvn$conf.int[1]),
                  UCI=c(pvs$conf.int[2],pva$conf.int[2],pvn$conf.int[2],vs$conf.int[2],va$conf.int[2],vn$conf.int[2],nvs$conf.int[2],nva$conf.int[2],nvn$conf.int[2]),
                  Test=rep(c("Saliva","Antigen","Nasal"),3),
                  ViralCulture=rep(c("pre-positive","positive","post-positive"),each=3))
shedse$ViralCulture=factor(shedse$ViralCulture,levels=c("pre-positive","positive","post-positive"))
levels(shedse$Test)=c("Antigen","Nasal\nRTqPCR","Saliva\nRTqPCR")
tiff("MMWR_Figure2.tif",width = 800,height = 800,res=100)
ggplot(shedse,aes(x=ViralCulture,y=Sensitivity,ymin=LCI,ymax=UCI,color=Test))+
  geom_pointrange(position=position_dodge(width=.5))+scale_color_manual(values = cbbPalette)+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "bottom",text=element_text(size=25))+
  ylab("Status Sensitivity")+scale_y_continuous(breaks = seq(0,1,.2),limits = c(0,1))
dev.off()

#Figure 3 and Table 3
load("positive_wide_data.Rdata")
load("positive_wide_infectious_data.Rdata")
###Sampling frequencies vs overall result
upreprint=unique(alldata$user_id[which(is.na(alldata$nasal_swab_virus)==F)])
se_f=matrix(NA,nrow=7,ncol = 13)
colnames(se_f)=c(paste(rep(c("Saliva","Antigen","Nasal"),each=3),rep(c("Sensitivity","LCI","UCI"),3)),"N",paste(c("Saliva","Antigen","Nasal"),"Npos"))
rownames(se_f)=c("q1d","q2d","q3d","q4d","q5d","q6d","q7d")
se_i=se_f
pdata=as.matrix(pdata)
pwinfdata=as.matrix(pwinfdata)
ur=pdata[,47]
mpreprint=which(pdata[,1]%in%upreprint)
pdata=pdata[mpreprint,]
pwinfdata=pwinfdata[which(pwinfdata[,1]%in%upreprint),]

#daily
s_res1=apply(pdata[,2:16],1,max,na.rm=T)
a_res1=apply(pdata[,18:31],1,max,na.rm=T)
n_res1=apply(pdata[,33:46],1,max,na.rm=T)
sft=binom.test(length(which(s_res1==1)),length(s_res1))
se_f[1,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res1==1)),length(a_res1))
se_f[1,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res1==1)),length(n_res1))
se_f[1,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_f[1,10]=sft$parameter
se_f[1,11]=sft$statistic
se_f[1,12]=aft$statistic
se_f[1,13]=nft$statistic
s_inf1=apply(pwinfdata[,2:16],1,max,na.rm=T)
a_inf1=apply(pwinfdata[,18:31],1,max,na.rm=T)
n_inf1=apply(pwinfdata[,33:46],1,max,na.rm=T)
sft=binom.test(length(which(s_inf1==1)),length(s_inf1))
se_i[1,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf1==1)),length(a_inf1))
se_i[1,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf1==1)),length(n_inf1))
se_i[1,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[1,10]=sft$parameter
se_i[1,11]=sft$statistic
se_i[1,12]=aft$statistic
se_i[1,13]=nft$statistic

#q2d
pdata2s=rbind(pdata[,seq(from=3,to=16,by=2)],pdata[,seq(from=4,to=16,by=2)])  
pdata2a=rbind(pdata[,seq(from=18,to=31,by=2)],pdata[,seq(from=19,to=31,by=2)])  
pdata2n=rbind(pdata[,seq(from=33,to=46,by=2)],pdata[,seq(from=34,to=46,by=2)])  
s_res2=apply(pdata2s,1,max,na.rm=T)
a_res2=apply(pdata2a,1,max,na.rm=T)
n_res2=apply(pdata2n,1,max,na.rm=T)
sft=binom.test(length(which(s_res2==1)),length(s_res2))
se_f[2,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res2==1)),length(a_res2))
se_f[2,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res2==1)),length(n_res2))
se_f[2,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_f[2,10]=sft$parameter
se_f[2,11]=sft$statistic
se_f[2,12]=aft$statistic
se_f[2,13]=nft$statistic
pwinfdata2s=rbind(pwinfdata[,seq(from=3,to=16,by=2)],pwinfdata[,seq(from=4,to=16,by=2)])  
pwinfdata2a=rbind(pwinfdata[,seq(from=18,to=31,by=2)],pwinfdata[,seq(from=19,to=31,by=2)])  
pwinfdata2n=rbind(pwinfdata[,seq(from=33,to=46,by=2)],pwinfdata[,seq(from=34,to=46,by=2)])  
s_inf2=apply(pwinfdata2s,1,max,na.rm=T)
a_inf2=apply(pwinfdata2a,1,max,na.rm=T)
n_inf2=apply(pwinfdata2n,1,max,na.rm=T)
sft=binom.test(length(which(s_inf2==1)),length(s_inf2))
se_i[2,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf2==1)),length(a_inf2))
se_i[2,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf2==1)),length(n_inf2))
se_i[2,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[2,10]=sft$parameter
se_i[2,11]=sft$statistic
se_i[2,12]=aft$statistic
se_i[2,13]=nft$statistic

#q3d
pdata3s=rbind(pdata[,seq(from=3,to=16,by=3)],pdata[,seq(from=4,to=16,by=3)],cbind(pdata[,seq(from=5,to=16,by=3)],NA))  
pdata3a=rbind(pdata[,seq(from=18,to=31,by=3)],pdata[,seq(from=19,to=31,by=3)],cbind(pdata[,seq(from=20,to=31,by=3)],NA))  
pdata3n=rbind(pdata[,seq(from=33,to=46,by=3)],pdata[,seq(from=34,to=46,by=3)],cbind(pdata[,seq(from=35,to=46,by=3)],NA))  
s_res3=apply(pdata3s,1,max,na.rm=T)
a_res3=apply(pdata3a,1,max,na.rm=T)
n_res3=apply(pdata3n,1,max,na.rm=T)
sft=binom.test(length(which(s_res3==1)),length(s_res3))
se_f[3,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res3==1)),length(a_res3))
se_f[3,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res3==1)),length(n_res3))
se_f[3,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_f[3,10]=sft$parameter
se_f[3,11]=sft$statistic
se_f[3,12]=aft$statistic
se_f[3,13]=nft$statistic
pwinfdata3s=rbind(pwinfdata[,seq(from=3,to=16,by=3)],pwinfdata[,seq(from=4,to=16,by=3)],cbind(pwinfdata[,seq(from=5,to=16,by=3)],NA))  
pwinfdata3a=rbind(pwinfdata[,seq(from=18,to=31,by=3)],pwinfdata[,seq(from=19,to=31,by=3)],cbind(pwinfdata[,seq(from=20,to=31,by=3)],NA))  
pwinfdata3n=rbind(pwinfdata[,seq(from=33,to=46,by=3)],pwinfdata[,seq(from=34,to=46,by=3)],cbind(pwinfdata[,seq(from=35,to=46,by=3)],NA))  
s_inf3=apply(pwinfdata3s,1,max,na.rm=T)
a_inf3=apply(pwinfdata3a,1,max,na.rm=T)
n_inf3=apply(pwinfdata3n,1,max,na.rm=T)
sft=binom.test(length(which(s_inf3==1)),length(s_inf3))
se_i[3,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf3==1)),length(a_inf3))
se_i[3,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf3==1)),length(n_inf3))
se_i[3,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[3,10]=sft$parameter
se_i[3,11]=sft$statistic
se_i[3,12]=aft$statistic
se_i[3,13]=nft$statistic

#q4d
pdata4s=rbind(pdata[,seq(from=3,to=16,by=4)],pdata[,seq(from=4,to=16,by=4)],cbind(pdata[,seq(from=5,to=16,by=4)],NA),cbind(pdata[,seq(from=6,to=16,by=4)],NA))  
pdata4a=rbind(pdata[,seq(from=18,to=31,by=4)],pdata[,seq(from=19,to=31,by=4)],cbind(pdata[,seq(from=20,to=31,by=4)],NA),cbind(pdata[,seq(from=21,to=31,by=4)],NA))  
pdata4n=rbind(pdata[,seq(from=33,to=46,by=4)],pdata[,seq(from=34,to=46,by=4)],cbind(pdata[,seq(from=35,to=46,by=4)],NA),cbind(pdata[,seq(from=36,to=46,by=4)],NA))  
s_res4=apply(pdata4s,1,max,na.rm=T)
a_res4=apply(pdata4a,1,max,na.rm=T)
n_res4=apply(pdata4n,1,max,na.rm=T)
sft=binom.test(length(which(s_res4==1)),length(s_res4))
se_f[4,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res4==1)),length(a_res4))
se_f[4,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res4==1)),length(n_res4))
se_f[4,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_f[4,10]=sft$parameter
se_f[4,11]=sft$statistic
se_f[4,12]=aft$statistic
se_f[4,13]=nft$statistic
pwinfdata4s=rbind(pwinfdata[,seq(from=3,to=16,by=4)],pwinfdata[,seq(from=4,to=16,by=4)],cbind(pwinfdata[,seq(from=5,to=16,by=4)],NA),cbind(pwinfdata[,seq(from=6,to=16,by=4)],NA))  
pwinfdata4a=rbind(pwinfdata[,seq(from=18,to=31,by=4)],pwinfdata[,seq(from=19,to=31,by=4)],cbind(pwinfdata[,seq(from=20,to=31,by=4)],NA),cbind(pwinfdata[,seq(from=21,to=31,by=4)],NA))  
pwinfdata4n=rbind(pwinfdata[,seq(from=33,to=46,by=4)],pwinfdata[,seq(from=34,to=46,by=4)],cbind(pwinfdata[,seq(from=35,to=46,by=4)],NA),cbind(pwinfdata[,seq(from=36,to=46,by=4)],NA))  
s_inf4=apply(pwinfdata4s,1,max,na.rm=T)
a_inf4=apply(pwinfdata4a,1,max,na.rm=T)
n_inf4=apply(pwinfdata4n,1,max,na.rm=T)
sft=binom.test(length(which(s_inf4==1)),length(s_inf4))
se_i[4,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf4==1)),length(a_inf4))
se_i[4,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf4==1)),length(n_inf4))
se_i[4,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[4,10]=sft$parameter
se_i[4,11]=sft$statistic
se_i[4,12]=aft$statistic
se_i[4,13]=nft$statistic

#q5d
pdata5s=rbind(pdata[,seq(from=3,to=16,by=5)],pdata[,seq(from=4,to=16,by=5)],
              pdata[,seq(from=5,to=16,by=5)],pdata[,seq(from=6,to=16,by=5)],
              cbind(pdata[,seq(from=7,to=16,by=5)],NA))  
pdata5a=rbind(pdata[,seq(from=18,to=31,by=5)],pdata[,seq(from=19,to=31,by=5)],
              pdata[,seq(from=20,to=31,by=5)],pdata[,seq(from=21,to=31,by=5)],
              cbind(pdata[,seq(from=22,to=31,by=5)],NA))  
pdata5n=rbind(pdata[,seq(from=33,to=46,by=5)],pdata[,seq(from=34,to=46,by=5)],
              pdata[,seq(from=35,to=46,by=5)],pdata[,seq(from=36,to=46,by=5)],
              cbind(pdata[,seq(from=37,to=46,by=5)],NA))  
s_res5=apply(pdata5s,1,max,na.rm=T)
a_res5=apply(pdata5a,1,max,na.rm=T)
n_res5=apply(pdata5n,1,max,na.rm=T)
sft=binom.test(length(which(s_res5==1)),length(s_res5))
se_f[5,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res5==1)),length(a_res5))
se_f[5,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res5==1)),length(n_res5))
se_f[5,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_f[5,10]=sft$parameter
se_f[5,11]=sft$statistic
se_f[5,12]=aft$statistic
se_f[5,13]=nft$statistic
pwinfdata5s=rbind(pwinfdata[,seq(from=3,to=16,by=5)],pwinfdata[,seq(from=4,to=16,by=5)],
                  pwinfdata[,seq(from=5,to=16,by=5)],pwinfdata[,seq(from=6,to=16,by=5)],
                  cbind(pwinfdata[,seq(from=7,to=16,by=5)],NA))  
pwinfdata5a=rbind(pwinfdata[,seq(from=18,to=31,by=5)],pwinfdata[,seq(from=19,to=31,by=5)],
                  pwinfdata[,seq(from=20,to=31,by=5)],pwinfdata[,seq(from=21,to=31,by=5)],
                  cbind(pwinfdata[,seq(from=22,to=31,by=5)],NA))  
pwinfdata5n=rbind(pwinfdata[,seq(from=33,to=46,by=5)],pwinfdata[,seq(from=34,to=46,by=5)],
                  pwinfdata[,seq(from=35,to=46,by=5)],pwinfdata[,seq(from=36,to=46,by=5)],
                  cbind(pwinfdata[,seq(from=37,to=46,by=5)],NA))  
s_inf5=apply(pwinfdata5s,1,max,na.rm=T)
a_inf5=apply(pwinfdata5a,1,max,na.rm=T)
n_inf5=apply(pwinfdata5n,1,max,na.rm=T)
sft=binom.test(length(which(s_inf5==1)),length(s_inf5))
se_i[5,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf5==1)),length(a_inf5))
se_i[5,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf5==1)),length(n_inf5))
se_i[5,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[5,10]=sft$parameter
se_i[5,11]=sft$statistic
se_i[5,12]=aft$statistic
se_i[5,13]=nft$statistic

#q6d
pdata6s=rbind(pdata[,seq(from=3,to=16,by=6)],pdata[,seq(from=4,to=16,by=6)],
              cbind(pdata[,seq(from=5,to=16,by=6)],NA),cbind(pdata[,seq(from=6,to=16,by=6)],NA),
              cbind(pdata[,seq(from=7,to=16,by=6)],NA),cbind(pdata[,seq(from=8,to=16,by=6)],NA))  
pdata6a=rbind(pdata[,seq(from=18,to=31,by=6)],pdata[,seq(from=19,to=31,by=6)],
              cbind(pdata[,seq(from=20,to=31,by=6)],NA),cbind(pdata[,seq(from=21,to=31,by=6)],NA),
              cbind(pdata[,seq(from=22,to=31,by=6)],NA),cbind(pdata[,seq(from=23,to=31,by=6)],NA))  
pdata6n=rbind(pdata[,seq(from=33,to=46,by=6)],pdata[,seq(from=34,to=46,by=6)],
              cbind(pdata[,seq(from=35,to=46,by=6)],NA),cbind(pdata[,seq(from=36,to=46,by=6)],NA),
              cbind(pdata[,seq(from=37,to=46,by=6)],NA),cbind(pdata[,seq(from=38,to=46,by=6)],NA))  
s_res6=apply(pdata6s,1,max,na.rm=T)
a_res6=apply(pdata6a,1,max,na.rm=T)
n_res6=apply(pdata6n,1,max,na.rm=T)
sft=binom.test(length(which(s_res6==1)),length(s_res6))
se_f[6,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res6==1)),length(a_res6))
se_f[6,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res6==1)),length(n_res6))
se_f[6,7:9]=c(nft$estimate,nft$conf.int[1:2])
s_t6=a_t6=n_t6=rep(NA,nrow(pdata6s))
se_f[6,10]=sft$parameter
se_f[6,11]=sft$statistic
se_f[6,12]=aft$statistic
se_f[6,13]=nft$statistic
pwinfdata6s=rbind(pwinfdata[,seq(from=3,to=16,by=6)],pwinfdata[,seq(from=4,to=16,by=6)],
                  cbind(pwinfdata[,seq(from=5,to=16,by=6)],NA),cbind(pwinfdata[,seq(from=6,to=16,by=6)],NA),
                  cbind(pwinfdata[,seq(from=7,to=16,by=6)],NA),cbind(pwinfdata[,seq(from=8,to=16,by=6)],NA))  
pwinfdata6a=rbind(pwinfdata[,seq(from=18,to=31,by=6)],pwinfdata[,seq(from=19,to=31,by=6)],
                  cbind(pwinfdata[,seq(from=20,to=31,by=6)],NA),cbind(pwinfdata[,seq(from=21,to=31,by=6)],NA),
                  cbind(pwinfdata[,seq(from=22,to=31,by=6)],NA),cbind(pwinfdata[,seq(from=23,to=31,by=6)],NA))  
pwinfdata6n=rbind(pwinfdata[,seq(from=33,to=46,by=6)],pwinfdata[,seq(from=34,to=46,by=6)],
                  cbind(pwinfdata[,seq(from=35,to=46,by=6)],NA),cbind(pwinfdata[,seq(from=36,to=46,by=6)],NA),
                  cbind(pwinfdata[,seq(from=37,to=46,by=6)],NA),cbind(pwinfdata[,seq(from=38,to=46,by=6)],NA))  
s_inf6=apply(pwinfdata6s,1,max,na.rm=T)
a_inf6=apply(pwinfdata6a,1,max,na.rm=T)
n_inf6=apply(pwinfdata6n,1,max,na.rm=T)
sft=binom.test(length(which(s_inf6==1)),length(s_inf6))
se_i[6,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf6==1)),length(a_inf6))
se_i[6,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf6==1)),length(n_inf6))
se_i[6,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[6,10]=sft$parameter
se_i[6,11]=sft$statistic
se_i[6,12]=aft$statistic
se_i[6,13]=nft$statistic

#q7d
pdata7s=rbind(pdata[,seq(from=3,to=16,by=7)],pdata[,seq(from=4,to=16,by=7)],
              pdata[,seq(from=5,to=16,by=7)],pdata[,seq(from=7,to=16,by=7)],
              pdata[,seq(from=7,to=16,by=7)],pdata[,seq(from=8,to=16,by=7)],
              pdata[,seq(from=9,to=16,by=7)])  
pdata7a=rbind(pdata[,seq(from=18,to=31,by=7)],pdata[,seq(from=19,to=31,by=7)],
              pdata[,seq(from=20,to=31,by=7)],pdata[,seq(from=21,to=31,by=7)],
              pdata[,seq(from=22,to=31,by=7)],pdata[,seq(from=23,to=31,by=7)],
              pdata[,seq(from=24,to=31,by=7)])  
pdata7n=rbind(pdata[,seq(from=33,to=46,by=7)],pdata[,seq(from=34,to=46,by=7)],
              pdata[,seq(from=35,to=46,by=7)],pdata[,seq(from=36,to=46,by=7)],
              pdata[,seq(from=37,to=46,by=7)],pdata[,seq(from=38,to=46,by=7)],
              pdata[,seq(from=39,to=46,by=7)])  
s_res7=apply(pdata7s,1,max,na.rm=T)
a_res7=apply(pdata7a,1,max,na.rm=T)
n_res7=apply(pdata7n,1,max,na.rm=T)
sft=binom.test(length(which(s_res7==1)),length(s_res7))
se_f[7,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_res7==1)),length(a_res7))
se_f[7,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_res7==1)),length(n_res7))
se_f[7,7:9]=c(nft$estimate,nft$conf.int[1:2])
s_t7=a_t7=n_t7=rep(NA,nrow(pdata7s))
se_f[7,10]=sft$parameter
se_f[7,11]=sft$statistic
se_f[7,12]=aft$statistic
se_f[7,13]=nft$statistic
pwinfdata7s=rbind(pwinfdata[,seq(from=3,to=16,by=7)],pwinfdata[,seq(from=4,to=16,by=7)],
                  pwinfdata[,seq(from=5,to=16,by=7)],pwinfdata[,seq(from=7,to=16,by=7)],
                  pwinfdata[,seq(from=7,to=16,by=7)],pwinfdata[,seq(from=8,to=16,by=7)],
                  pwinfdata[,seq(from=9,to=16,by=7)])  
pwinfdata7a=rbind(pwinfdata[,seq(from=18,to=31,by=7)],pwinfdata[,seq(from=19,to=31,by=7)],
                  pwinfdata[,seq(from=20,to=31,by=7)],pwinfdata[,seq(from=21,to=31,by=7)],
                  pwinfdata[,seq(from=22,to=31,by=7)],pwinfdata[,seq(from=23,to=31,by=7)],
                  pwinfdata[,seq(from=24,to=31,by=7)])  
pwinfdata7n=rbind(pwinfdata[,seq(from=33,to=46,by=7)],pwinfdata[,seq(from=34,to=46,by=7)],
                  pwinfdata[,seq(from=35,to=46,by=7)],pwinfdata[,seq(from=36,to=46,by=7)],
                  pwinfdata[,seq(from=37,to=46,by=7)],pwinfdata[,seq(from=38,to=46,by=7)],
                  pwinfdata[,seq(from=39,to=46,by=7)])  
s_inf7=apply(pwinfdata7s,1,max,na.rm=T)
a_inf7=apply(pwinfdata7a,1,max,na.rm=T)
n_inf7=apply(pwinfdata7n,1,max,na.rm=T)
sft=binom.test(length(which(s_inf7==1)),length(s_inf7))
se_i[7,1:3]=c(sft$estimate,sft$conf.int[1:2])
aft=binom.test(length(which(a_inf7==1)),length(a_inf7))
se_i[7,4:6]=c(aft$estimate,aft$conf.int[1:2])
nft=binom.test(length(which(n_inf7==1)),length(n_inf7))
se_i[7,7:9]=c(nft$estimate,nft$conf.int[1:2])
se_i[7,10]=sft$parameter
se_i[7,11]=sft$statistic
se_i[7,12]=aft$statistic
se_i[7,13]=nft$statistic

se_f[,1:9]=round(se_f[,1:9],3)
se_i[,1:9]=round(se_i[,1:9],3)
se_t=cbind(se_f[,c(10,4)],se_i[,4],se_f[,12],se_i[,12],se_f[,1],se_i[,1],se_f[,11],se_i[,11],se_f[,7],se_i[,7],se_f[,13],se_i[,13])
colnames(se_t)=c("N","AgSe","AgISe","AgN","AgIN","SSe","SISe","SN","SIN","NSe","NISe","NN","NIN")
write.csv(se_t,file="MMWR_Table3.csv")

se_fl=rbind(se_f[,1:3],se_f[,4:6],se_f[,7:9]);colnames(se_fl)=c("Sensitivity","LCI","UCI")
se_fl=as.data.frame(se_fl)
se_fl$frequency=rep(rownames(se_f),3);se_fl$frequency=factor(se_fl$frequency)
se_fl$test=rep(c("Saliva","Antigen","Nasal"),each=nrow(se_f));se_fl$test=factor(se_fl$test)

se_il=rbind(se_i[,1:3],se_i[,4:6],se_i[,7:9]);colnames(se_il)=c("Sensitivity","LCI","UCI")
se_il=as.data.frame(se_il)
se_il$frequency=rep(rownames(se_i),3);se_il$frequency=factor(se_il$frequency)
se_il$test=rep(c("Saliva","Antigen","Nasal"),each=nrow(se_i));se_il$test=factor(se_il$test)

se_all=rbind(se_fl,se_il)
se_all$Outcome=rep(c("Detection of Infection","Detection before or\nwhile viral culture positive"),each=nrow(se_fl))
se_all$Outcome=factor(se_all$Outcome,levels = levels(factor(se_all$Outcome))[c(2,1)])
levels(se_all$frequency)=c("Daily","Every\nOther\nDay","Every\nThird\nDay","Every\nFourth\nDay","Every\nFifth\nDay","Every\nSixth\nDay","Weekly")
levels(se_all$test)=c("Antigen","Nasal\nRTqPCR","Saliva\nRTqPCR")
tiff("MMWR_Figure3.tif",width = 800,height = 1200,res = 150)
ggplot(se_all,aes(x=frequency,y=Sensitivity,ymin=LCI,ymax=UCI,color=test))+geom_pointrange(position = position_dodge(width = .5))+
  theme(text=element_text(size=15),legend.title = element_blank(),legend.position = "top")+
  scale_color_manual(values = cbbPalette)+
  facet_grid(Outcome~.)+xlab(NULL)+ylab("Protocol Sensitivity")+scale_y_continuous(breaks = seq(0,1,.2),limits = c(0,1))
dev.off()




