require(MASS)
library(CASdatasets)
library(tidyverse)
library(gamlss.tr)
library(copula)
library(sampling)
library(AUC)
library(fitdistrplus)
library(gamlss)
library(cplm)

rm(list=ls())

gen.trun(par = c(0),
         family = PO,
         type = "left") ## generate the truncated poisson distribution functions
ggplot2::theme_set(theme_bw())
ggplot2::theme_update(text = element_text(size = 20))

# frequency data

data(ausprivauto0405) 
dat1 <- ausprivauto0405[ausprivauto0405$VehValue!=0,]
# dat1 <- dat1[dat1$Exposure > 0.2, ]
dat1 <- dat1[dat1$Exposure > 0.7, ]

aggregate(dat1$VehBody,by=list(dat1$VehBody),length)
dat1$VehBody<-as.character(dat1$VehBody)
dat1$VehBody[dat1$VehBody=="Panel van"]<-"Van"
dat1$VehBody[dat1$VehBody=="Station wagon"]<-"Wagon"
dat1$VehBody[dat1$VehBody=="Motorized caravan"]<-"Caravan"
dat1$VehBody[dat1$VehBody=="Bus"]<-"Others"
dat1$VehBody[dat1$VehBody=="Convertible"]<-"Others"
dat1$VehBody[dat1$VehBody=="Roadster"]<-"Others"
dat1$VehBody<-as.factor(dat1$VehBody)
dat1$VehBody <- relevel(dat1$VehBody, ref = 'Sedan')

sum(dat1$ClaimOcc)/nrow(dat1);sum(dat1$ClaimNb!=0)/nrow(dat1)
dat1$VehAge<-as.character(dat1$VehAge)
dat1$VehAge[dat1$VehAge=="old cars"]<-"old"
dat1$VehAge[dat1$VehAge=="young cars"]<-"young"
dat1$VehAge[dat1$VehAge=="oldest cars"]<-"oldest"
dat1$VehAge[dat1$VehAge=="youngest cars"]<-"youngest"
dat1$VehAge<-as.factor(dat1$VehAge)
dat1$VehAge <- relevel(dat1$VehAge, ref = 'old')

dat1$DrivAge<-as.character(dat1$DrivAge)
dat1$DrivAge[dat1$DrivAge=="young people"]<-"young"
dat1$DrivAge[dat1$DrivAge=="older work. people"]<-"older_work."
dat1$DrivAge[dat1$DrivAge=="oldest people"]<-"oldest"
dat1$DrivAge[dat1$DrivAge=="working people"]<-"working"
dat1$DrivAge[dat1$DrivAge=="old people"]<-"old"
dat1$DrivAge[dat1$DrivAge=="youngest people"]<-"youngest"
dat1$DrivAge<-as.factor(dat1$DrivAge)
dat1$DrivAge <- relevel(dat1$DrivAge, ref = 'older_work.')

dat1$Gender <- relevel(dat1$Gender, ref = "Female")

dat1$LnVehValue<-log(dat1$VehValue)
hist(dat1$VehValue)
hist(dat1$LnVehValue)

# dat1$ClaimNb_rc <- pmin(dat1$ClaimNb,2) #right censored at 2

# data split

aggregate(dat1$ClaimNb,by=list(dat1$ClaimNb),length)
dat1$ind<-NA
dat1$ind[dat1$ClaimNb==0]<-rep(1:5,13000)[1:sum(dat1$ClaimNb==0)]
dat1$ind[dat1$ClaimNb==1]<-rep(1:5,13000)[1:sum(dat1$ClaimNb==1)]
dat1$ind[dat1$ClaimNb==2]<-rep(1:5,13000)[1:sum(dat1$ClaimNb==2)]
dat1$ind[dat1$ClaimNb>2]<-rep(1:5,13000)[1:sum(dat1$ClaimNb>2)]
aggregate(dat1$Exposure,by=list(dat1$ind),length)
aggregate(dat1$Exposure,by=list(dat1$ind),sum)
aggregate(dat1$ClaimNb,by=list(dat1$ind),sum)
aggregate(dat1$ClaimNb,by=list(dat1$ind),sum)$x/aggregate(dat1$Exposure,by=list(dat1$ind),sum)$x

# design matrix
dat1_mat<-model.matrix(~., data=dat1)[,-1]
str(dat1)
4*11*2*6
dat1_mat2<-as.data.frame(model.matrix(~ Exposure + ClaimNb + ind + LnVehValue*VehAge*VehBody*Gender*DrivAge, data=dat1)[,-1])
dim(dat1_mat2)
names(dat1_mat2)[1:50]
cat_exposure<-dat1_mat2$Exposure*(dat1_mat2[,5:ncol(dat1_mat2)]!=0)
cat_exposure_sum<-apply(cat_exposure,2,sum)
hist(cat_exposure_sum)
sum(cat_exposure_sum>1000)
ind_e<-names(cat_exposure_sum)[which(cat_exposure_sum>1000)]
dat1_mat2<-dat1_mat2[,c(1:4,which(names(dat1_mat2)%in%ind_e))]
names(dat1_mat2)
#LnVehValue_interact<-dat1_mat2$LnVehValue*dat1_mat2[,c(5:ncol(dat1_mat2))]
#dat1_mat2<-data.frame(dat1_mat2,LnVehValue_interact)

head(dat1_mat)
dat1_mat<-data.frame(dat1_mat)

# severity data

dat2 <- dat1 %>% 
  filter(ClaimOcc==1) %>% 
  mutate(ClaimSev = ClaimAmount/ClaimNb)
dat2$LnSev<-log(dat2$ClaimSev)
aggregate(dat2$Exposure,by=list(dat2$ind),length)
aggregate(dat2$Exposure,by=list(dat2$ind),sum)
aggregate(dat2$ClaimNb,by=list(dat2$ind),sum)
aggregate(dat2$ClaimSev,by=list(dat2$ind),mean)

# design matrix
dat2_mat<-model.matrix(~., data=dat2)[,-1]
head(dat2_mat)
dat2_mat<-data.frame(dat2_mat)

dat2_mat2<-as.data.frame(
  model.matrix(
    ~ Exposure + ClaimNb + ClaimAmount + ClaimSev + LnSev + ind + LnVehValue*VehAge*VehBody*Gender*DrivAge, data=dat2)[,-1])
dim(dat2_mat2)
names(dat2_mat2)[1:10]
cat_exposure<-dat2_mat2$Exposure*(dat2_mat2[,8:ncol(dat2_mat2)]!=0)
cat_exposure_sum<-apply(cat_exposure,2,sum)
hist(cat_exposure_sum)
sum(cat_exposure_sum>100)
ind_e<-names(cat_exposure_sum)[which(cat_exposure_sum>100)]
dat2_mat2<-dat2_mat2[,c(1:7,which(names(dat2_mat2)%in%ind_e))]
names(dat2_mat2)
#LnVehValue_interact<-dat2_mat2$LnVehValue*dat2_mat2[,c(8:ncol(dat2_mat2))]
#dat2_mat2<-data.frame(dat2_mat2,LnVehValue_interact)

# Gamma Deviance Loss

Gamma.Deviance <-function(pred, obs, wei){
  2*sum(wei*((obs/pred)-1-log(obs/pred)))/length(pred)
}

# Poisson Deviance Loss

Poisson.Deviance <- function(pred, obs){
  2*(sum(pred)-sum(obs)+sum(log((obs/pred)^(obs))))/length(pred)
}

# Minus ZTP log-likelihood

minus_ZTP_logLik <- function(pred, obs){
  -sum(dPOtr(obs, mu=pred,log = TRUE)) / length(pred)
}

