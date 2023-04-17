#Variable processing#
#data$Pathological.findings<-factor(data$Pathological.findings,levels = c(0,1),labels = c("Benign","Malignant"))
data$Proportion.of.cystic.components<-data$Proportion.of.cystic.components
data$Echogenicity<-factor(data$Echogenicity,levels = c(0,1),labels = c("Non-hypoechoic","Hypoechoic"))
data$Shape<-factor(data$Shape,levels = c(0,1),labels = c("Wider than tall","Taller than wide"))
data$Ultrasound.Margin<-factor(data$Ultrasound.Margin,levels = c(0,1),labels = c("Regular","Irregular"))
data$Large.comet.tail.artifacts<-factor(data$Large.comet.tail.artifacts,levels = c(0,1),labels = c("No","Yes"))
data$Macrocalcifications<-factor(data$Macrocalcifications,levels = c(0,1),labels = c("No","Yes"))
data$Peripheral..rim..calcifications<-factor(data$Peripheral..rim..calcifications,levels = c(0,1),labels = c("No","Yes"))
data$Punctate.echogenic.foci<-factor(data$Punctate.echogenic.foci,levels = c(0,1),labels = c("No","Yes"))
data$Score<-data$Score
data$TIRADS<-factor(data$TIRADS,levels = c(1,2,3,4,5),labels = c("1","2","3","4","5"))
data$Diameter<-data$Diameter
data$Management.of.ACR.TI.RADS<-factor(data$Management.of.ACR.TI.RADS,levels = c(0,1,2),labels = c("No FNA or follow-up","Follow-up","FNA"))
data$Posterior.echogenicity.of.nodules<-factor(data$Posterior.echogenicity.of.nodules,levels = c(0,1,2),labels = c("No change","Enhancement","Attenuation or shadow"))
data$Acoustic.halo<-factor(data$Acoustic.halo,levels = c(0,1,2),labels = c("No","Partial","Complete"))
data$Elastography<-factor(data$Elastography,levels = c(0,1),labels = c("Soft","Stiff"))
data$Blood.volume<-factor(data$Blood.volume,levels = c(0,1),labels = c("Other","Low"))
data$Enhancement.uniformity<-factor(data$Enhancement.uniformity,levels = c(0,1),labels = c("Homogeneous","Heterogeneous"))
data$Enhancement.margin<-factor(data$Enhancement.margin,levels = c(0,1),labels = c("Regular","Irregular"))
data$Ring.enhancement<-factor(data$Ring.enhancement,levels = c(0,1,2),labels = c("No","Hypo","Hyper"))
data$Change.of.nodule.size<-factor(data$Change.of.nodule.size,levels = c(0,1,2),labels = c("Decreased","Fixed","Increased"))
data$Age<-data$Age
data$Gender<-factor(data$Gender,levels = c(0,1),labels = c("Famale","Male"))
data$BMI<-data$BMI
data$SBP<-data$SBP
data$DBP<-data$DBP
data$Thyroid.diffuse.lesions<-factor(data$Thyroid.diffuse.lesions,levels = c(0,1),labels = c("No","Yes"))
#Univariable analysis#
uni_glm_model<-function(x){
  FML<-as.formula(paste0("Pathological.findings==1~",x))
  glm1<-glm(FML,data = data,family = binomial)
  glm2<-summary(glm1)
  OR<-round(exp(coef(glm1)),2)
  SE<-round(glm2$coefficients[,2],3)  
  CI2.5<-round(exp(coef(glm1)-1.96*SE),2)
  CI97.5<-round(exp(coef(glm1)+1.96*SE),2)
  CI<-paste0(CI2.5,'-',CI97.5)
  B<-round(glm2$coefficients[,1],3)
  Z<-round(glm2$coefficients[,3],3)
  P<-round(glm2$coefficients[,4],4)
  uni_glm_model<-data.frame('characteristics'=x,
                            'B'=B,
                            'SE'=SE,
                            'OR'=OR,
                            'CI'=CI,
                            'Z' =Z,
                            'P'=P)[-1,]
  
  return(uni_glm_model)
}
variable.names<-colnames(data)[c(3,5,7,9,10,12,14:26,34,35,37:44)]  
variable.names
uni_glm<-lapply(variable.names,uni_glm_model)
uni_glm
library(plyr)
uni_glm<-ldply(uni_glm,data.frame)
uni_glm
View(uni_glm)
write.csv(uni_glm, "uni.csv")
uni_glm1 <- uni_glm[uni_glm$P<= 0.05,]
uni_glm1
uni_glm$characteristics[uni_glm$P<= 0.05]
write.csv(uni_glm1, "p5.csv")
#Multivariable analysis and model construction#
fml<- as.formula(paste0('Pathological.findings==1~',paste0(uni_glm$characteristics[uni_glm$P<0.05],collapse = '+'))) 
fml
fmlm<-as.formula(Pathological.findings == 1 ~ Elastography + Proportion.of.cystic.components + Echogenicity + Shape 
                 + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                 + Posterior.echogenicity.of.nodules +  Acoustic.halo + Blood.volume 
                 + Enhancement.margin + Ring.enhancement + Age)
fmlCEUS<-as.formula(Pathological.findings == 1 ~  Proportion.of.cystic.components + Echogenicity + Shape 
                    + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                    + Posterior.echogenicity.of.nodules +  Acoustic.halo + Blood.volume 
                    + Enhancement.margin + Ring.enhancement + Age)
fmlSE<-as.formula(Pathological.findings == 1 ~ Elastography + Proportion.of.cystic.components + Echogenicity + Shape 
                  + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                  + Posterior.echogenicity.of.nodules +  Acoustic.halo + Age)
fmlUS<-as.formula(Pathological.findings == 1 ~  Proportion.of.cystic.components + Echogenicity + Shape 
                  + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                  + Posterior.echogenicity.of.nodules +  Acoustic.halo +  Age)
fml0<-as.formula(Pathological.findings == 1 ~  Management.of.ACR.TI.RADS)
#forward#
fmodel<-glm(Pathological.findings~1,data = data,family=binomial)
fmodel
Fmodelm<-step(fmodel,scope=list(upper=~ Elastography + Proportion.of.cystic.components + Echogenicity + Shape 
                                + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                                + Posterior.echogenicity.of.nodules +  Acoustic.halo + Blood.volume 
                                + Enhancement.margin + Ring.enhancement + Age,
                                lower=~1),data = data,family=binomial,direction ="forward")
Fmodelm
summary(Fmodelm)
FmodelCEUS<-step(fmodel,scope=list(upper=~  Proportion.of.cystic.components + Echogenicity + Shape 
                                   + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                                   + Posterior.echogenicity.of.nodules +  Acoustic.halo + Blood.volume 
                                   + Enhancement.margin + Ring.enhancement + Age,
                                   lower=~1),data = data,family=binomial,direction ="forward")
FmodelCEUS
summary(FmodelCEUS)
FmodelSE<-step(fmodel,scope=list(upper=~ Elastography + Proportion.of.cystic.components + Echogenicity + Shape 
                                 + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                                 + Posterior.echogenicity.of.nodules +  Acoustic.halo  
                                 + Enhancement.margin + Ring.enhancement + Age,
                                 lower=~1),data = data,family=binomial,direction ="forward")
FmodelSE
summary(FmodelSE)
FmodelUS<-step(fmodel,scope=list(upper=~ Proportion.of.cystic.components + Echogenicity + Shape 
                                 + Ultrasound.Margin + Macrocalcifications + Punctate.echogenic.foci + Diameter
                                 + Posterior.echogenicity.of.nodules +  Acoustic.halo 
                                 + Enhancement.margin + Ring.enhancement + Age,
                                 lower=~1),data = data,family=binomial,direction ="forward")
FmodelUS
summary(FmodelUS)
Fmodel0<-step(fmodel,scope=list(upper=~ Management.of.ACR.TI.RADS,lower=~1),data = data,family=binomial,direction ="forward")
Fmodel0
summary(Fmodel0)
cbind(coef=coef(Fmodelm),confint(Fmodelm))
cbind(coef=coef(FmodelCEUS),confint(FmodelCEUS))
cbind(coef=coef(FmodelSE),confint(FmodelSE))
cbind(coef=coef(FmodelUS),confint(FmodelUS))
cbind(coef=coef(Fmodel0),confint(Fmodel0))
exp(cbind(coef=coef(Fmodelm),confint(Fmodelm)))
exp(cbind(coef=coef(FmodelCEUS),confint(FmodelCEUS)))
exp(cbind(coef=coef(FmodelSE),confint(FmodelSE)))
exp(cbind(coef=coef(FmodelUS),confint(FmodelUS)))
exp(cbind(coef=coef(Fmodel0),confint(Fmodel0)))
AIC(Fmodel0,Fmodelm,FmodelCEUS,FmodelSE,FmodelUS)
anova(Fmodel0,FmodelUS,test = "Chisq")
anova(Fmodel0,FmodelSE,test = "Chisq")
anova(Fmodel0,FmodelCEUS,test = "Chisq")
anova(Fmodel0,Fmodelm,test = "Chisq")
anova(FmodelUS,FmodelSE,test = "Chisq")
anova(FmodelUS,FmodelCEUS,test = "Chisq")
anova(FmodelUS,Fmodelm,test = "Chisq")
anova(FmodelSE,FmodelCEUS,test = "Chisq")
anova(FmodelSE,Fmodelm,test = "Chisq")
anova(FmodelCEUS,Fmodelm,test = "Chisq")
glmm<-summary(Fmodelm)
glmm
glmCEUS<-summary(FmodelCEUS)
glmCEUS
glmSE<-summary(FmodelSE)
glmSE
glmUS<-summary(FmodelUS)
glmUS
glm0<-summary(Fmodel0)
glm0
glmm$coefficients
OR<-round(exp(glmm$coefficients[,1]),2)
OR
SE<-round(glmm$coefficients[,2],3)
CI2.5<-round(exp(coef(Fmodelm)-1.96*SE),2)
CI97.5<-round(exp(coef(Fmodelm)+1.96*SE),2)
CI<-paste0(CI2.5,'-',CI97.5)
B<-round(glmm$coefficients[,1],3)
Z<-round(glmm$coefficients[,3],3)
P<-round(glmm$coefficients[,4],3)
mlogitm<-data.frame( 
  'B'=B,
  'SE'=SE,
  'OR'=OR,
  'CI'=CI,
  'Z' =Z,
  'P'=P)[-1,]   
mlogitm
View(mlogitm)
names(data)
fml
fmlm
multinamesm<-as.character(colnames(data)[c(3,7,9,10,12,18,34,40)])
multinamesm
mlogitm<-data.frame('characteristics'=multinamesm,mlogitm)
mlogitm
View(mlogitm)
write.csv(mlogitm, "multim.csv")
finalm<-merge.data.frame(uni_glm,mlogitm,by='characteristics',all = T,sort = T)
finalm
View(finalm)
write.csv(finalm, "finalm.csv")
glmCEUS$coefficients
OR<-round(exp(glmCEUS$coefficients[,1]),2)
OR
SE<-round(glmCEUS$coefficients[,2],3)
CI2.5<-round(exp(coef(FmodelCEUS)-1.96*SE),2)
CI97.5<-round(exp(coef(FmodelCEUS)+1.96*SE),2)
CI<-paste0(CI2.5,'-',CI97.5)
B<-round(glmCEUS$coefficients[,1],3)
Z<-round(glmCEUS$coefficients[,3],3)
P<-round(glmCEUS$coefficients[,4],3)
mlogitm<-data.frame( 
  'B'=B,
  'SE'=SE,
  'OR'=OR,
  'CI'=CI,
  'Z' =Z,
  'P'=P)[-1,]   
mlogitm
View(mlogitm)
names(data)
fml
fmlCEUS
multinamesCEUS<-as.character(colnames(data)[c(3,7,9,10,12,18,40)])
multinamesCEUS
mlogitm<-data.frame('characteristics'=multinamesCEUS,mlogitm)
mlogitm
View(mlogitm)
write.csv(mlogitm, "multiCEUS.csv")
finalCEUS<-merge.data.frame(uni_glm,mlogitm,by='characteristics',all = T,sort = T)
finalCEUS
View(finalCEUS)
write.csv(finalCEUS, "finalCEUS.csv")
glmSE$coefficients
OR<-round(exp(glmSE$coefficients[,1]),2)
OR
SE<-round(glmSE$coefficients[,2],3)
CI2.5<-round(exp(coef(FmodelSE)-1.96*SE),2)
CI97.5<-round(exp(coef(FmodelSE)+1.96*SE),2)
CI<-paste0(CI2.5,'-',CI97.5)
B<-round(glmSE$coefficients[,1],3)
Z<-round(glmSE$coefficients[,3],3)
P<-round(glmSE$coefficients[,4],3)
mlogitm<-data.frame( 
  'B'=B,
  'SE'=SE,
  'OR'=OR,
  'CI'=CI,
  'Z' =Z,
  'P'=P)[-1,]   
mlogitm
View(mlogitm)
names(data)
fml
fmlSE
multinamesSE<-as.character(colnames(data)[c(3,7,9,10,12,18,34)])
multinamesSE
mlogitm<-data.frame('characteristics'=multinamesSE,mlogitm)
mlogitm
View(mlogitm)
write.csv(mlogitm, "multiSE.csv")
finalSE<-merge.data.frame(uni_glm,mlogitm,by='characteristics',all = T,sort = T)
finalSE
View(finalSE)
write.csv(finalSE, "finalSE.csv")
glmUS$coefficients
OR<-round(exp(glmUS$coefficients[,1]),2)
OR
SE<-round(glmUS$coefficients[,2],3)
CI2.5<-round(exp(coef(FmodelUS)-1.96*SE),2)
CI97.5<-round(exp(coef(FmodelUS)+1.96*SE),2)
CI<-paste0(CI2.5,'-',CI97.5)
B<-round(glmUS$coefficients[,1],3)
Z<-round(glmUS$coefficients[,3],3)
P<-round(glmUS$coefficients[,4],3)
mlogitm<-data.frame( 
  'B'=B,
  'SE'=SE,
  'OR'=OR,
  'CI'=CI,
  'Z' =Z,
  'P'=P)[-1,]   
mlogitm
View(mlogitm)
names(data)
fml
fmlUS
multinamesUS<-as.character(colnames(data)[c(3,7,9,10,12,18)])
multinamesUS
mlogitm<-data.frame('characteristics'=multinamesUS,mlogitm)
mlogitm
View(mlogitm)
write.csv(mlogitm, "multiUS.csv")
finalUS<-merge.data.frame(uni_glm,mlogitm,by='characteristics',all = T,sort = T)
finalUS
View(finalUS)
write.csv(finalUS, "finalUS.csv")