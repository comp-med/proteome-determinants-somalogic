## Set working directory
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Julia/proteome_determinants/bin/")

options(stringsAsFactors = F)

rm(list=ls())

## load packages needed 
library(readstata13)
library(RNOmni)
library(caret)
library(glmnet)
library(doMC)
library(tidyr)
library(dplyr)
library(data.table)

## read in parameter file 
params <- read.delim("mrc_seqid_param_file.txt", sep = "\ ", he=F) 
head(params)

# fs.mat <- read.delim(paste0("../data_output/feature_selection/",params$V1[1],"_FS_determinants.txt"), sep = "\t", header = T)
# fs.mat$var <- row.names(fs.mat)
# fs.mat <- fs.mat[,c("var","select")]
# fs.mat$var[grep("cis",fs.mat$var)] <- "cis.score"
# fs.mat$var[grep("trans",fs.mat$var)] <- "trans.score"
# colnames(fs.mat) <- c("var",params$V1[1])
# 
# for (i in params$V1[2:length(params$V1)]){
#   tmp <- read.delim(paste0("../data_output/feature_selection/",i,"_FS_determinants.txt"), sep = "\t", header = T)
#   tmp$var <- row.names(tmp)
#   tmp <- tmp[,c("var","select")]
#   tmp$var[grep("cis",tmp$var)] <- "cis.score"
#   tmp$var[grep("trans",tmp$var)] <- "trans.score"
# 
#   colnames(tmp) <- c("var", i)
#   fs.mat <- merge(fs.mat, tmp, by = "var", all = T)
# }

# write.table(fs.mat, file = "../data_output/Featrue_selection_proteome_determinants.txt", sep = "\t", row.names = F)
fs.mat <- read.delim("../data_output/Featrue_selection_proteome_determinants.txt", sep = "\t", he=T)


fs.gather <- gather(fs.mat, prot, select, 2:ncol(fs.mat))
head(fs.gather)
fs.gather$select <- ifelse(fs.gather$var=="cis.score"&is.na(fs.gather$select),0,fs.gather$select)
fs.gather$select <- ifelse(fs.gather$var=="trans.score"&is.na(fs.gather$select),0,fs.gather$select)
fs.gather$select <- fs.gather

## generate an order
prot.order <- fs.gather%>%
  group_by(var)%>%
  summarise(mean=mean(select, na.rm=T))
fs.gather <- merge(fs.gather, prot.order, by = "var")
library(ggplot2)
library(forcats)

fs.gather$group <- ifelse(fs.gather$var%in%c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                             "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                             "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass"),"Anthropometric",
                          ifelse(fs.gather$var%in%c("alkphos0","alt0","liver_score","Bilirubin0"),"Liver",
                                  ifelse(fs.gather$var%in%c("chol0","hdl0","ldl0","triglyceride0"),"Lipids",
                                         ifelse(fs.gather$var%in%c("A1c_pc_all","glucose0","Insulin"),"Glycaemic",
                                                ifelse(fs.gather$var%in%c("TSH","FT3","FT4"),"Thyroid",
                                                       ifelse(fs.gather$var%in%c("CAD.GRS","DBP.GRS","SBP.GRS","eGFR.GRS",
                                                                                 "T2D.GRS_scale","FG.GRS_scale","h2PG.GRS_scale","IR.GRS_scale",
                                                                                 "BMI.GRS_scale","HIP.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale"),"Disease-GRS",
                                                              ifelse(fs.gather$var=="cis.score","Cis-pQTL_score",
                                                                     ifelse(fs.gather$var=="trans.score","Trans-pQTL_score",
                                                                            ifelse(fs.gather$var%in%c("Alcohol","g_smoke_current","g_smoke_ever","PAEE","g_vitamincresult",
                                                                                                      "tMDS_2","DASH_AccordanceScore","Energy_kcal"),"Dietary_lifestyle",
                                                                                   ifelse(fs.gather$var%in%c("BPDia1","BPSys1","eGFR"),"Cardio",
                                                                                          ifelse(fs.gather$var%in%c("anti_htn","anti_lipid"),"Medication",
                                                                                                 ifelse(fs.gather$var=="g_hs_crp","Inflammatory",
                                                                                                        ifelse(fs.gather$var=="age","Age",
                                                                                                               ifelse(fs.gather$var=="sex","Sex",
                                                                                                                      ifelse(fs.gather$var%in%c("TestSite.2","TestSite.3",
                                                                                                                                                "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                                                                                                                                "PC1.20","PC1.5","PC1.005"),"Techincal",NA)))))))))))))))


x.labs <- list("bmi"="BMI","whr"="WHR","subcutaneousfat_qc1"="Subcutaneous fat","visceralfatg_qc1"="Visceral fat","peripheralfat_qc1"="Peripheral fat",
               "AD04i_iDEXA_Total_bone_mass"="Bone mass","AD05i_iDEXA_Total_fat_mass"="Total fat mass","AD06i_iDEXA_Total_lean_mass"="Lean mass",
               "AD14i_iDEXA_arms_fat_mass"="Arms fat mass","AD38i_iDEXA_legs_fat_mass"="Legs fat mass",
               "alkphos0"="ALP","alt0"="ALT","Bilirubin0"="Bilirubin","chol0"="Total cholesterol","g_hs_crp"="CRP",
               "g_vitamincresult"="Vitamin C","hdl0"="HDL","ldl0"="LDL","triglyceride0"="Triglycerides",
               "BPDia1"="Diastolic BP","BPSys1"="Systolic BP","eGFR"="eGFR",
               "tMDS_2"="MDS","DASH_AccordanceScore"="DASH","Energy_kcal"="Total energy intake",
               "A1c_pc_all"="HbA1c","glucose0"="Glucose","Insulin"="Insulin","liver_score"="Liver score",
               "CAD.GRS"="CAD-GRS","DBP.GRS"="DBP-GRS","SBP.GRS"="SBP-GRS","eGFR.GRS"="eGFR-GRS",
               "T2D.GRS_scale"="T2D-GRS","FG.GRS_scale"="FG-GRS","h2PG.GRS_scale"="2hPG-GRS","IR.GRS_scale"="FI-GRS",
               "BMI.GRS_scale"="BMI-GRS","HIP.GRS_scale"="Hip-GRS","WAIST.GRS_scale"="Waist-GRS","WHR.GRS_scale"="WHR-GRS",
               "cis.score"="cis-pQTL score","trans.score"="trans-pQTL score",
               "Alcohol"="Alcohol intake","g_smoke_current"="Current smoker","g_smoke_ever"="Ever smoker",
               "PAEE"="PAEE","anti_htn"="Anti-hypertensive medication","anti_lipid"="Lipid lowering medication",
               "TSH"="TSH","FT3"="T3","FT4"="T4",
               "age"="Age","sex"="Sex",
               "TestSite.2"="Test site 2","TestSite.3"="Test site 3",
               "AppDate_Attended_ROUNDED"="Appointment date","AppDate.spring"="Appointment date spring",
               "AppDate.summer"="Appointment date summer","AppDate.winter"="Appointment date winter",
               "PC1.20"="PC1 20% proteins","PC1.5"="PC1 5% proteins","PC1.005"="PC1 0.005% proteins")

var.expl <- read.delim("../data_output/variance_explained_all_proteins.txt", sep ="\t", header = T)
var.gather <- gather(var.expl, prot, var.ex, 2:ncol(var.expl))
head(var.gather)

fs.count <- var.gather %>%
  filter(is.na(var.ex)!=T)%>%
  group_by(prot)%>%
  count()

fs.gather <- var.gather %>%
  filter(is.na(var.ex)!=T)%>%
  group_by(var)%>%
  count()
fs.gather$group <- ifelse(fs.gather$var%in%c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                             "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                             "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass"),"Anthropometric",
                          ifelse(fs.gather$var%in%c("alkphos0","alt0","liver_score","Bilirubin0"),"Liver",
                                 ifelse(fs.gather$var%in%c("chol0","hdl0","ldl0","triglyceride0"),"Lipids",
                                        ifelse(fs.gather$var%in%c("A1c_pc_all","glucose0","Insulin"),"Glycaemic",
                                               ifelse(fs.gather$var%in%c("TSH","FT3","FT4"),"Thyroid",
                                                      ifelse(fs.gather$var%in%c("CAD.GRS","DBP.GRS","SBP.GRS","eGFR.GRS",
                                                                                "T2D.GRS_scale","FG.GRS_scale","h2PG.GRS_scale","IR.GRS_scale",
                                                                                "BMI.GRS_scale","HIP.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale"),"Disease-GRS",
                                                             ifelse(fs.gather$var=="cis.score","Cis-pQTL_score",
                                                                    ifelse(fs.gather$var=="trans.score","Trans-pQTL_score",
                                                                           ifelse(fs.gather$var%in%c("Alcohol","g_smoke_current","g_smoke_ever","PAEE","g_vitamincresult",
                                                                                       "tMDS_2","DASH_AccordanceScore","Energy_kcal"),"Dietary_lifestyle",
                                                                                  ifelse(fs.gather$var%in%c("BPDia1","BPSys1","eGFR"),"Cardio",
                                                                                         ifelse(fs.gather$var%in%c("anti_htn","anti_lipid"),"Medication",
                                                                                                ifelse(fs.gather$var=="g_hs_crp","Inflammatory",
                                                                                                       ifelse(fs.gather$var=="age","Age",
                                                                                                              ifelse(fs.gather$var=="sex","Sex",
                                                                                                                     ifelse(fs.gather$var%in%c("TestSite.2","TestSite.3",
                                                                                                                                 "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                                                                                                                 "PC1.20","PC1.5","PC1.005"),"Techincal",NA)))))))))))))))


var.gather$group <- ifelse(var.gather$var%in%c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                             "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                             "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass"),"Anthropomteric",
                          ifelse(var.gather$var%in%c("alkphos0","alt0","liver_score","Bilirubin0"),"Liver",
                                 ifelse(var.gather$var%in%c("chol0","hdl0","ldl0","triglyceride0"),"Lipids",
                                        ifelse(var.gather$var%in%c("A1c_pc_all","glucose0","Insulin"),"Glycaemic",
                                               ifelse(var.gather$var%in%c("TSH","FT3","FT4"),"Thyroid",
                                                      ifelse(var.gather$var%in%c("CAD.GRS","DBP.GRS","SBP.GRS","eGFR.GRS",
                                                                                "T2D.GRS_scale","FG.GRS_scale","h2PG.GRS_scale","IR.GRS_scale",
                                                                                "BMI.GRS_scale","HIP.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale"),"Disease-GRS",
                                                             ifelse(var.gather$var=="cis.score","Cis-pQTL_score",
                                                                    ifelse(var.gather$var=="trans.score","Trans-pQTL_score",
                                                                           ifelse(var.gather$var%in%c("Alcohol","g_smoke_current",
                                                                                                      "g_smoke_ever","PAEE","g_vitamincresult",
                                                                                                      "tMDS_2","DASH_AccordanceScore"
                                                                                                      ,"Energy_kcal"),"Dietary_lifestyle",
                                                                                  ifelse(var.gather$var%in%c("BPDia1","BPSys1","eGFR"),"Cardio",
                                                                                         ifelse(var.gather$var%in%c("anti_htn","anti_lipid"),"Medication",
                                                                                                ifelse(var.gather$var=="g_hs_crp","Inflammatory",
                                                                                                       ifelse(var.gather$var=="age","Age",
                                                                                                              ifelse(var.gather$var=="sex","Sex",
                                                                                                                     ifelse(var.gather$var%in%c("TestSite.2","TestSite.3",
                                                                                                                                 "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                                                                                                                 "PC1.20","PC1.5","PC1.005"),"Techincal",NA)))))))))))))))



var.gather <- merge(var.gather, fs.gather[,1:2], by = "var")
library(ggbreak)
var.gather$var.ex <- var.gather$var.ex*100
v.plot<- var.gather %>% 
  filter(!var%in%c("noise.1","noise.2","noise.3")&group!="Techincal") %>% 
  mutate(var.o= fct_reorder(var,var.ex, median, na.rm=T)) %>% 
  ggplot(aes(x=var.o, y=sqrt(var.ex), colour=group))+
  geom_boxplot()+
  # scale_y_break(c(0.125,0.8))+
  coord_flip()+ggtitle("a.")+
  scale_x_discrete(labels=x.labs)+theme_bw()+
  scale_y_continuous(breaks = sqrt(c(1,2,5,10,25,50,75)), labels=c(1,2,5,10,25,50,75))+
  ylab("Explained variance (%)")+
  labs(color="Determinant group")+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.position = "none",
    plot.title = element_text(face="bold"))+
  # facet_wrap(.~group, scales="free")+
  scale_colour_manual(values = c("Age"="#D7B5A6" ,"Sex"="#9D7660",
                                 "Anthropomteric"="#B07AA1","Lipids"="#D37295",
                                 "Cardio"="#E15759","Medication"="#FF9D9A",
                                 "Dietary_lifestyle"="#F1CE63",
                                 "Disease-GRS"="#86BCB6","Cis-pQTL_score"="#4E79A7","Trans-pQTL_score"="#A0CBE8",
                                 "Glycaemic"="#59A14F",
                                 "Inflammatory"="#F28E2B","Thyroid"="#FFBE7D",
                                 "Liver"="#D4A6C8"))

ord.var.m <- var.gather %>% 
  filter(!var%in%c("noise.1","noise.2","noise.3")&group!="Techincal") %>% 
  mutate(var.o= fct_reorder(var,var.ex, median, na.rm=T)) %>% 
  select(var.o)
a <- fs.gather %>% 
  filter(!var%in%c("noise.1","noise.2","noise.3")&group!="Techincal") %>% 
  mutate(var.o = factor(var, levels = levels(ord.var.m$var.o))) %>% 
  ggplot(aes(x=var.o, y=n, fill=group))+
  geom_bar(stat = "identity")+
  coord_flip()+ggtitle("b.")+
  scale_x_discrete(labels=x.labs)+theme_bw()+
  ylab("Number of proteins selected for")+
  labs(fill="Variable group")+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), legend.position = c(0.8,0.2), 
        plot.title = element_text(face = "bold"), legend.key.size = unit(0.09,"cm"))+
  scale_fill_manual(values = c("Age"="#D7B5A6" ,"Sex"="#9D7660",
                               "Anthropometric"="#B07AA1","Lipids"="#D37295",
                               "Cardio"="#E15759","Medication"="#FF9D9A",
                               "Dietary_lifestyle"="#F1CE63",
                               "Disease-GRS"="#86BCB6","Cis-pQTL_score"="#4E79A7","Trans-pQTL_score"="#A0CBE8",
                               "Glycaemic"="#59A14F",
                               "Inflammatory"="#F28E2B","Thyroid"="#FFBE7D",
                               "Liver"="#D4A6C8"
  ))


library(gridExtra)
var.res.all <- subset(var.expl, var!="n"&var!="Residuals"&var!="TestSite.2"&var!="TestSite.3"&var!="AppDate_Attended_ROUNDED"&var!="AppDate.spring"&
                        var!="AppDate.summer"&var!="AppDate.winter"&var!="PC1.20"&var!="PC1.5"&var!="PC1.005")
var.res.gather <- gather(var.res.all, prot, exp.var, 2:4980)

## order variables  
var.res.gather$var <- factor(var.res.gather$var,
                             levels = c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                        "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                        "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass",
                                        "alkphos0","alt0","Bilirubin0","chol0","g_hs_crp","g_vitamincresult","hdl0","ldl0","triglyceride0",
                                        "BPDia1","BPSys1","eGFR",
                                        "tMDS_2","DASH_AccordanceScore","Energy_kcal",
                                        "A1c_pc_all","glucose0","Insulin","liver_score",
                                        "CAD.GRS","DBP.GRS","SBP.GRS","eGFR.GRS",
                                        "T2D.GRS_scale","FG.GRS_scale","h2PG.GRS_scale","IR.GRS_scale",
                                        "BMI.GRS_scale","HIP.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale",
                                        "cis.score","trans.score",
                                        "Alcohol","g_smoke_current","g_smoke_ever","PAEE","anti_htn","anti_lipid",
                                        "TSH","FT3","FT4",
                                        "age","sex"#,
                                        # "TestSite.2","TestSite.3",
                                        # "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                        # "PC1.20","PC1.5","PC1.005"
                             ))
library(forcats)
order.prots <- var.res.gather%>%
  group_by(prot)%>%
  summarise(total =sum(exp.var, na.rm = T)*100 )%>%
  mutate(prot.ord=fct_reorder(prot,desc(total)))
var.res.gather$prot <- factor(var.res.gather$prot, levels = levels(order.prots$prot.ord)) 
## plot stacked bar plot
library(ggplot2)

## collapsed barplot
var.res.gather$group <- ifelse(var.res.gather$var%in%c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                                       "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                                       "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass"),"Anthropometric",
                               ifelse(var.res.gather$var%in%c("alkphos0","alt0","liver_score","Bilirubin0"),"Liver",
                                      ifelse(var.res.gather$var%in%c("chol0","hdl0","ldl0","triglyceride0"),"Lipids",
                                             ifelse(var.res.gather$var%in%c("A1c_pc_all","glucose0","Insulin"),"Glycaemic",
                                                    ifelse(var.res.gather$var%in%c("TSH","FT3","FT4"),"Thyroid",
                                                           ifelse(var.res.gather$var%in%c("BMI.GRS_scale","CAD.GRS","DBP.GRS","eGFR.GRS","FG.GRS_scale","h2PG.GRS_scale","HIP.GRS_scale",
                                                                                      "IR.GRS_scale","SBP.GRS","T2D.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale"),"Disease-GRS",
                                                                  ifelse(var.res.gather$var=="cis.score","Cis-pQTL_score",
                                                                         ifelse(var.res.gather$var=="trans.score","Trans-pQTL_score",
                                                                                ifelse(var.res.gather$var%in%c("Alcohol","g_smoke_current",
                                                                                                               "g_smoke_ever","PAEE",
                                                                                                               "g_vitamincresult","tMDS_2","DASH_AccordanceScore","Energy_kcal"),"Dietary_lifestyle",
                                                                                       ifelse(var.res.gather$var%in%c("BPDia1","BPSys1","eGFR"),"Cardio",
                                                                                              ifelse(var.res.gather$var%in%c("anti_htn","anti_lipid"),"Medication",
                                                                                                     ifelse(var.res.gather$var=="g_hs_crp","Inflammatory",
                                                                                                            ifelse(var.res.gather$var=="age","Age",
                                                                                                                   ifelse(var.res.gather$var=="sex","Sex",
                                                                                                                          ifelse(var.res.gather$var%in%c("TestSite.2","TestSite.3",
                                                                                                                                           "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                                                                                                                           "PC1.20","PC1.5","PC1.005"),"Techincal",NA)))))))))))))))



var.res.gather$group <- factor(var.res.gather$group, levels = rev(c("Cis-pQTL_score","Trans-pQTL_score","Disease-GRS","Age","Sex",
                                                                "Anthropometric","Glycaemic","Inflammatory","Lipids","Liver", "Cardio","Thyroid","Dietary_lifestyle","Medication")))

prots.10 <- var.res.gather %>%
  group_by(prot) %>%
  summarise(total = sum(exp.var, na.rm=T)) %>%
  filter(total>0.10)

c <- ggplot(subset(var.res.gather, prot%in%prots.10$prot), aes(x = prot, y = exp.var*100, fill=group))+
  theme_bw()+
  geom_bar(stat = "identity", position = position_stack())+
  ggtitle("b.")+
  xlab("Protein")+
  # scale_y_break(c(20,22), 0.25)+
  theme(plot.title = element_text(face = "bold"), #axis.title.x = element_blank(),
        legend.text = element_text(size = 10),axis.title.x = element_blank(),
        legend.key.size = unit(0.7,"line"), axis.ticks.x = element_blank(),#axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = c(0.75,0.65),axis.text.x = element_blank())+
  #   scale_x_discrete(labels=x.labs)+
  scale_fill_manual(values = list("Age"="#D7B5A6" ,"Sex"="#9D7660",
                                  "Anthropometric"="#B07AA1","Lipids"="#D37295",
                                  "Cardio"="#E15759","Medication"="#FF9D9A",
                                  "Dietary_lifestyle"="#F1CE63",
                                  "Disease-GRS"="#86BCB6","Cis-pQTL_score"="#4E79A7","Trans-pQTL_score"="#A0CBE8",
                                  "Glycaemic"="#59A14F",
                                  "Inflammatory"="#F28E2B","Thyroid"="#FFBE7D",
                                  "Liver"="#D4A6C8"))+
  labs(fill="Variable group")+guides(fill = guide_legend(ncol = 2))+
  ylab("Explained variance (%)")+scale_y_continuous(limits = c(0,100),expand = c(0,0))
 
c2 <- var.res.gather %>% 
  group_by(prot) %>% 
  mutate(freq = 100*(exp.var / sum(exp.var, na.rm = T))) %>% 
  filter(prot%in%prots.10$prot) %>% 
  ggplot(aes(x = prot, y = freq, fill=group))+
  theme_bw()+
  geom_bar(stat = "identity", position = position_stack())+
  # ggtitle("a.")+
  xlab("Protein")+
  # scale_y_break(c(20,22), 0.25)+
  theme(plot.title = element_text(face = "bold"), #axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.7,"line"), axis.ticks.x = element_blank(),#axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",axis.text.x = element_blank())+
  #   scale_x_discrete(labels=x.labs)+
  scale_fill_manual(values = list("Age"="#D7B5A6" ,"Sex"="#9D7660",
                                  "Anthropometric"="#B07AA1","Lipids"="#D37295",
                                  "Cardio"="#E15759","Medication"="#FF9D9A",
                                  "Dietary_lifestyle"="#F1CE63",
                                  "Disease-GRS"="#86BCB6","Cis-pQTL_score"="#4E79A7","Trans-pQTL_score"="#A0CBE8",
                                  "Glycaemic"="#59A14F",
                                  "Inflammatory"="#F28E2B","Thyroid"="#FFBE7D",
                                  "Liver"="#D4A6C8"))+
  # labs(fill="Explanatory variable group")+
  ylab("Proportional contribution (%)")+scale_y_continuous(limits = c(0,100),expand = c(0,0))
l.pa <- grid.arrange(c,c2, nrow=2, heights = c(0.4,0.25))

## including technical variation
c2.t <- var.expl  %>% 
  filter( var!="n"&var!="Residuals") %>% 
  gather( prot, exp.var, 2:4980) %>%
  mutate(group=ifelse(var%in%c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                              "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                              "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass"),"Anthropometric",
                      ifelse(var%in%c("alkphos0","alt0","liver_score","Bilirubin0"),"Liver",
                             ifelse(var%in%c("chol0","hdl0","ldl0","triglyceride0"),"Lipids",
                                    ifelse(var%in%c("A1c_pc_all","glucose0","Insulin"),"Glycaemic",
                                           ifelse(var%in%c("TSH","FT3","FT4"),"Thyroid",
                                                  ifelse(var%in%c("BMI.GRS_scale","CAD.GRS","DBP.GRS","eGFR.GRS","FG.GRS_scale","h2PG.GRS_scale","HIP.GRS_scale",
                                                                                 "IR.GRS_scale","SBP.GRS","T2D.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale"),"Disease-GRS",
                                                         ifelse(var=="cis.score","Cis-pQTL_score",
                                                                ifelse(var=="trans.score","Trans-pQTL_score",
                                                                       ifelse(var%in%c("Alcohol","g_smoke_current",
                                                                                                      "g_smoke_ever","PAEE",
                                                                                                      "g_vitamincresult","tMDS_2","DASH_AccordanceScore","Energy_kcal"),"Dietary_lifestyle",
                                                                              ifelse(var%in%c("BPDia1","BPSys1","eGFR"),"Cardio",
                                                                                     ifelse(var%in%c("anti_htn","anti_lipid"),"Medication",
                                                                                            ifelse(var=="g_hs_crp","Inflammatory",
                                                                                                   ifelse(var=="age","Age",
                                                                                                          ifelse(var=="sex","Sex",
                                                                                                                 ifelse(var%in%c("TestSite.2","TestSite.3",
                                                                                                                                 "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                                                                                                                 "PC1.20","PC1.5","PC1.005"),"Technical",NA))))))))))))))),
         group= factor(group, levels = rev(c("Cis-pQTL_score","Trans-pQTL_score","Disease-GRS","Age","Sex",
                                                     "Anthropometric","Glycaemic","Inflammatory","Lipids","Liver", "Cardio","Thyroid","Dietary_lifestyle","Medication","Technical"))))%>% 
  group_by(prot) %>% 
  mutate(freq = 100*(exp.var / sum(exp.var, na.rm = T))) %>% 
  # filter(prot%in%as.character(prots.10$prot)) %>% 
  ggplot(aes(x = prot, y = freq, fill=group))+
  theme_bw()+
  geom_bar(stat = "identity", position = position_stack())+
  # ggtitle("a.")+
  xlab("Protein")+
  # scale_y_break(c(20,22), 0.25)+
  theme(plot.title = element_text(face = "bold"), #axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.7,"line"), axis.ticks.x = element_blank(),#axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",axis.text.x = element_blank())+
  #   scale_x_discrete(labels=x.labs)+
  scale_fill_manual(values = list("Age"="#D7B5A6" ,"Sex"="#9D7660",
                                  "Anthropometric"="#B07AA1","Lipids"="#D37295",
                                  "Cardio"="#E15759","Medication"="#FF9D9A",
                                  "Dietary_lifestyle"="#F1CE63",
                                  "Disease-GRS"="#86BCB6","Cis-pQTL_score"="#4E79A7","Trans-pQTL_score"="#A0CBE8",
                                  "Glycaemic"="#59A14F",
                                  "Inflammatory"="#F28E2B","Thyroid"="#FFBE7D",
                                  "Liver"="#D4A6C8","Technical"="darkgrey"))+
  # labs(fill="Explanatory variable group")+
  ylab("Proportional contribution (%)")+scale_y_continuous(limits = c(0,100),expand = c(0,0))
l.pa <- grid.arrange(c,c2.t, nrow=2, heights = c(0.4,0.25))


########################################################
####  UMAP exploration using proteome determinants  ####
########################################################

##check optimal number of cluster
# library(factoextra)
# y <- var.expl[-which(var.expl$var%in%c("TestSite.2","TestSite.3",
#                                        "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
#                                        "PC1.20","PC1.5","PC1.005")),]
# ## replace missing by 0
# y <- sapply(params$V1, function(x){
#   ifelse(is.na(y[,x]),0,y[,x])
# })
# 
# fviz_nbclust(y, FUNcluster = kmeans, method = "gap_stat",nboot = 50, k.max = 50)

library(umap)
## UMAP exlucding technical variation 
u.labels.s <- sapply(params$V1, function(x){
  ii <- which.max(var.res.all[,x])
  deter <- var.res.all$var[ii]
  return(deter)
})

u.data.s <- var.res.all
row.names(u.data.s) <- u.data.s$var
u.data.s <- as.data.frame(t(u.data.s[,2:ncol(u.data.s)]))
head(u.data.s)

u.data.s[,1:ncol(u.data.s)] <- sapply(1:ncol(u.data.s), function(x){
  ifelse(is.na(u.data.s[,x]),0,u.data.s[,x])
})


## set configureation to run umap
costum.config <- umap.defaults
costum.config$random_state <- 10
costum.config$metric <- "pearson"
costum.config$n_epochs <- 1000
costum.config$input <- "data"
costum.config$init <- "random"
# umap.list <- list()
# 
# for (i in 5:25){
#   print(i)
#   costum.config$n_neighbors=i
#   
#   umap.deter.pheno <- umap(u.data.s, config = costum.config)
#   df.pheno <- data.frame(x=umap.deter.pheno$layout[,1], y=umap.deter.pheno$layout[,2], determinant=u.labels.s)
#   df.pheno$determinant <- factor(df.pheno$determinant, levels = c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
#                                                                   "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
#                                                                   "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass",
#                                                                   "alkphos0","alt0","Bilirubin0","chol0","g_hs_crp","g_vitamincresult","hdl0","ldl0","triglyceride0",
#                                                                   "BPSys1","eGFR",
#                                                                   "DASH_AccordanceScore",
#                                                                   "A1c_pc_all","glucose0","Insulin","liver_score",
#                                                                   "DBP.GRS","HIP.GRS_scale",
#                                                                   "cis.score","trans.score",
#                                                                   "Alcohol","g_smoke_current","g_smoke_ever","anti_htn","anti_lipid",
#                                                                   "TSH","FT4",
#                                                                   "age","sex"))
#   umap.list[[i]] <- ggplot(df.pheno, aes(x,y,colour=determinant))+
#     geom_point()+theme_bw()+ggtitle(i)+
#     scale_color_manual(values = list("bmi"=cols[1],"whr"=cols[2],"subcutaneousfat_qc1"=cols[3],"visceralfatg_qc1"=cols[4],"peripheralfat_qc1"=cols[5],
#                                      "AD04i_iDEXA_Total_bone_mass"=cols[6],"AD05i_iDEXA_Total_fat_mass"=cols[7],"AD06i_iDEXA_Total_lean_mass"=cols[8],
#                                      "AD14i_iDEXA_arms_fat_mass"=cols[9],"AD38i_iDEXA_legs_fat_mass"=cols[10],
#                                      "alkphos0"=cols[11],"alt0"=cols[12],"Bilirubin0"=cols[13],"chol0"=cols[14],"g_hs_crp"=cols[15],
#                                      "g_vitamincresult"=cols[16],"hdl0"=cols[17],"ldl0"=cols[18],"triglyceride0"=cols[19],
#                                      "BPDia1"=cols[20],"BPSys1"=cols[21],"eGFR"=cols[22],
#                                      "tMDS_2"=cols[23],"DASH_AccordanceScore"=cols[24],"Energy_kcal"=cols[25],
#                                      "A1c_pc_all"=cols[26],"glucose0"=cols[27],"Insulin"=cols[28],"liver_score"=cols[29],
#                                      "CAD.GRS"=cols[30],"DBP.GRS"=cols[31],"SBP.GRS"=cols[32],"eGFR.GRS"=cols[33],
#                                      "T2D.GRS_scale"=cols[34],"FG.GRS_scale"=cols[35],"h2PG.GRS_scale"=cols[36],"IR.GRS_scale"=cols[37],
#                                      "BMI.GRS_scale"=cols[38],"HIP.GRS_scale"=cols[39],"WAIST.GRS_scale"=cols[40],"WHR.GRS_scale"=cols[41],
#                                      "cis.score"=cols[42],"trans.score"=cols[43],
#                                      "Alcohol"=cols[44],"g_smoke_current"=cols[45],"g_smoke_ever"=cols[46],
#                                      "PAEE"=cols[47],"anti_htn"=cols[48],"anti_lipid"=cols[49],
#                                      "TSH"=cols[50],"FT3"=cols[51],"FT4"=cols[52],
#                                      "age"=cols[53],"sex"=cols[54],
#                                      "TestSite.2"=cols[55],"TestSite.3"=cols[56],
#                                      "AppDate_Attended_ROUNDED"=cols[57],"AppDate.spring"=cols[58],"AppDate.summer"=cols[59],"AppDate.winter"=cols[60],
#                                      "PC1.20"=cols[61],"PC1.5"=cols[62],"PC1.005"=cols[63]))+
#     theme(legend.position="none"
#       # legend.position = c(0.27,0.83),legend.text = element_text(size = 5),legend.title = element_blank(),
#           # legend.key.size = unit(0.7,"line")
#       )+
#     xlab("UMAP 1")+ylab("UMAP 2")+guides(color=guide_legend(ncol = 3))
# }
# 
# do.call(grid.arrange, umap.list[5:25])


## plot final configurations
costum.config$n_neighbors=13

umap.deter.pheno <- umap(u.data.s, config = costum.config)

df.pheno <- data.frame(x=umap.deter.pheno$layout[,1], y=umap.deter.pheno$layout[,2], determinant=u.labels.s)
df.pheno$determinant <- factor(df.pheno$determinant, levels = c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                                                "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                                                "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass",
                                                                "alkphos0","alt0","Bilirubin0","chol0","g_hs_crp","g_vitamincresult","hdl0","ldl0","triglyceride0",
                                                                "BPSys1","eGFR",
                                                                "DASH_AccordanceScore",
                                                                "A1c_pc_all","glucose0","Insulin","liver_score",
                                                                "DBP.GRS","HIP.GRS_scale",
                                                                "cis.score","trans.score",
                                                                "Alcohol","g_smoke_current","g_smoke_ever","anti_htn","anti_lipid",
                                                                "TSH","FT4",
                                                                "age","sex"))

## now with high-level determinants
df.pheno$group <- ifelse(df.pheno$determinant%in%c("bmi","whr","subcutaneousfat_qc1","visceralfatg_qc1","peripheralfat_qc1",
                                                   "AD04i_iDEXA_Total_bone_mass","AD05i_iDEXA_Total_fat_mass","AD06i_iDEXA_Total_lean_mass",
                                                   "AD14i_iDEXA_arms_fat_mass","AD38i_iDEXA_legs_fat_mass"),"Anthropometric",
                         ifelse(df.pheno$determinant%in%c("alkphos0","alt0","liver_score","Bilirubin0"),"Liver",
                                ifelse(df.pheno$determinant%in%c("chol0","hdl0","ldl0","triglyceride0"),"Lipids",
                                       ifelse(df.pheno$determinant%in%c("A1c_pc_all","glucose0","Insulin"),"Glycaemic",
                                              ifelse(df.pheno$determinant%in%c("TSH","FT3","FT4"),"Thyroid",
                                                     ifelse(df.pheno$determinant%in%c("CAD.GRS","DBP.GRS","SBP.GRS","eGFR.GRS",
                                                                                "T2D.GRS_scale","FG.GRS_scale","h2PG.GRS_scale","IR.GRS_scale",
                                                                                "BMI.GRS_scale","HIP.GRS_scale","WAIST.GRS_scale","WHR.GRS_scale"),"Disease-GRS",
                                                            ifelse(df.pheno$determinant=="cis.score","Cis-pQTL_score",
                                                                   ifelse(df.pheno$determinant=="trans.score","Trans-pQTL_score",
                                                                          ifelse(df.pheno$determinant%in%c("Alcohol","g_smoke_current",
                                                                                                           "g_smoke_ever","PAEE",
                                                                                                           "g_vitamincresult","tMDS_2","DASH_AccordanceScore",
                                                                                                           "Energy_kcal"),"Dietary_lifestyle",
                                                                                 ifelse(df.pheno$determinant%in%c("BPDia1","BPSys1","eGFR"),"Cardio",
                                                                                        ifelse(df.pheno$determinant%in%c("anti_htn","anti_lipid"),"Medication",
                                                                                               ifelse(df.pheno$determinant=="g_hs_crp","Inflammatory",
                                                                                                      ifelse(df.pheno$determinant=="age","Age",
                                                                                                             ifelse(df.pheno$determinant=="sex","Sex",
                                                                                                                    ifelse(df.pheno$determinant%in%c("TestSite.2","TestSite.3",
                                                                                                                                       "AppDate_Attended_ROUNDED","AppDate.spring","AppDate.summer","AppDate.winter",
                                                                                                                                       "PC1.20","PC1.5","PC1.005"),"Techincal",NA)))))))))))))))



u.plot <- ggplot(df.pheno, aes(x,y,colour=group))+
  geom_point(shape=16)+theme_bw()+ggtitle("a.")+
  scale_color_manual(values = list("Age"="#D7B5A6" ,"Sex"="#9D7660",
                                   "Anthropometric"="#B07AA1","Lipids"="#D37295",
                                   "Cardio"="#E15759","Medication"="#FF9D9A",
                                   "Dietary_lifestyle"="#F1CE63",
                                   "Disease-GRS"="#86BCB6","Cis-pQTL_score"="#4E79A7","Trans-pQTL_score"="#A0CBE8",
                                   "Glycaemic"="#59A14F",
                                   "Inflammatory"="#F28E2B","Thyroid"="#FFBE7D",
                                   "Liver"="#D4A6C8"))+
  theme(legend.position = "none",plot.title = element_text(face = "bold"),
        legend.text = element_text(size = 7),legend.title = element_text(size = 8),
        legend.key.size = unit(0.7,"line"))+
  xlab("UMAP 1")+ylab("UMAP 2")+
  guides(color=guide_legend(ncol = 3, title = "Explanatory variable group"))+
  labs(color="Explanatory variable group")
# p.top <- grid.arrange(v2.1, v2.2,ncol=2)

# pdf("../figures/FIG2_VExplSectors_umap.pdf", width = 9, height=9)
# grid.arrange( p.top, u.plot, nrow=2, heights=c(0.35,0.4))
# dev.off()
pdf("../figures/FIG3_FS_21032023.pdf", width = 9, height=6)
grid.arrange(v.plot,a, nrow=1)
dev.off()
# pdf("../figures/FIG2_FS_VExpl_umap.pdf", width = 9, height=9)
# lo.p <- grid.arrange( a, v.plot,nrow=1)
# dev.off()
pdf("../figures/FIG2_VExpl_freq_umap_Tehcnical.pdf", width = 9, height=12)
grid.arrange( u.plot, l.pa,nrow=2)
dev.off()
