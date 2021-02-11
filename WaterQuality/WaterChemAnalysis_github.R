#Written by Hannah Aichelman
#Contact: hannahaichelman@gmail.com
#Last updated: 2/11/2021

#This script is used to plot water chemistry data associated with the time course physiology manuscript
#as well as summarize the parameters through time, which is reported in supplementary tables 1 and 2 
#and plotted in supplementary figures 1 and 2.

#set working directory
setwd("~/Documents/BU/NEU_experiment/waterchem")

#load necessary libraries
library("ggplot2")
library("segmented")
library("plotrix")
library("lubridate")
library("chron")
library("Hmisc")
library("pgirmess")
library("lsmeans")
library("MASS")
library("Rmisc")
library("dplyr")
library("tidyverse")
library("readxl")


#read in new sheet with more data
calc_data <- read_excel("compiledwaterchem_github.xlsx", sheet = "calculated")
msd_data <- read_excel("compiledwaterchem_github.xlsx", sheet = "measured")

#reorganize levels of treatment factor
calc_data$target_treat = factor(calc_data$target_treat)
levels(calc_data$target_treat)
calc_data$target_treat = factor(calc_data$target_treat,levels = c("280_28", "280_31", "400_28", "400_31", "700_28", "700_31", "2800_28", "2800_31"))
#levels(data$treat) <- c("ambient pCO2, 28°C","Xlow pCO2, 31°C","Xlow pCO2, 28°C","ambient pCO2, 31°C","next century pCO2, 28°C","next century pCO2, 31°C","extreme pCO2, 28°C","extreme pCO2, 31°C")

msd_data$target_treat = factor(msd_data$target_treat)
levels(msd_data$target_treat)
msd_data$target_treat = factor(msd_data$target_treat,levels = c("280_28", "280_31", "400_28", "400_31", "700_28", "700_31", "2800_28", "2800_31"))

msd_data = msd_data %>%
  mutate(mtemp = as.numeric(mtemp), pHm = as.numeric(pHm), sal = as.numeric(sal)) %>%
  mutate(msmt = as.factor(msmt), tank = as.factor(tank), temp = as.factor(temp), pco2 = as.factor(pco2), tp = as.factor(tp), tp_combined = as.factor(tp_combined))
  
# remove data from August 25th (measurement 2) because it is incomplete  
calc_data = calc_data %>%
  mutate(mpco2 = as.numeric(mpco2), alk = as.numeric(alk), dic = as.numeric(dic), pHc = as.numeric(pHc), 
         hco3 = as.numeric(hco3), co2 = as.numeric(co2), co3 = as.numeric(co3), sat = as.numeric(sat)) %>%
  mutate(msmt = as.factor(msmt), tank = as.factor(tank), temp = as.factor(temp), pco2 = as.factor(pco2), tp = as.factor(tp), tp_combined = as.factor(tp_combined)) %>%
  filter(msmt != "2") #remove second measurement between T0 and T30 because of incomplete data

#### WATER TEMP ####

# plot average water temperature by time point and treatment
#plot 3 tanks separately:
tempPlot <- ggplot(msd_data, aes(x = tp, y = mtemp, fill = target_treat, shape = tank)) +
  geom_point(aes(color = target_treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Mean Temperature")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(target_treat~., ncol=2)
tempPlot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs_newdata/Temp_row.pdf", tempPlot, width=6, height=8, units=c("in"), useDingbats=FALSE)

str(msd_data)
head(msd_data)

# filter out the treatments not considered here
temp_data <- msd_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

summary(temp_data)

#explore the data
T0_temp = temp_data %>%
  filter(tp == "T0")

meantemp_T0 <- summarySE(T0_temp, measurevar="mtemp", groupvars=c("target_treat"), na.rm=TRUE)

lmT0=aov(mtemp~temp*pco2, data=T0_temp)
summary(lmT0)
par(mfrow=c(2,2))
plot(lmT0)

##these data are not normal at T0
lmT0_log=aov(log(mtemp)~temp*pco2, data=T0_temp)
summary(lmT0_log)
par(mfrow=c(2,2))
plot(lmT0_log)

TukeyHSD(lmT0) #all different at T0

T30_temp_cum = temp_data %>%
  filter(tp == "T0" | tp == "T30")

meantemp_T30 <- summarySE(T30_temp_cum, measurevar="mtemp", groupvars=c("target_treat"), na.rm=TRUE)

lmT30=aov(mtemp~temp*pco2, data=T30_temp_cum)
summary(lmT30)
par(mfrow=c(2,2))
plot(lmT30) #normal data

TukeyHSD(lmT30) #all different at T30

T60_temp_cum = temp_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60")

meantemp_T60 <- summarySE(T60_temp_cum, measurevar="mtemp", groupvars=c("target_treat"), na.rm=TRUE)

lmT60=aov(mtemp~temp*pco2, data=T60_temp_cum)
summary(lmT60)
par(mfrow=c(2,2))
plot(lmT60) #normal data

TukeyHSD(lmT60) #all different at T60

T90_temp_cum = temp_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90")

meantemp_T90 <- summarySE(T90_temp_cum, measurevar="mtemp", groupvars=c("target_treat"), na.rm=TRUE)

lmT90=aov(mtemp~temp*pco2, data=T90_temp_cum)
summary(lmT90)
par(mfrow=c(2,2))
plot(lmT90) #normal data

TukeyHSD(lmT90) #all different at T30

# plot boxplots of cumulative treatments
tempPlot <- ggplot(T90_temp_cum, aes(x = target_treat, y = mtemp)) +
  geom_boxplot(fill="gray")+
  geom_jitter(shape=16, position=position_jitter(0.1), size=2)+ 
  #ggtitle("Mean Temperature")+
  xlab("Treatment")+
  ylab("Temperature (°C)")+
  ylim(27,32)+
  theme(legend.position = "none")+
  theme_bw()
tempPlot

#SummarySE, average temperature by time point
#Throughout, use SummarySE to get data summaries that were used in supplementary tables
meanTemp <- summarySE(msd_data, measurevar="mtemp", groupvars=c("target_treat","tp"), na.rm=TRUE)
summary(meanTemp)

meanTemp <- summarySE(msd_data, measurevar="mtemp", groupvars=c("target_treat"), na.rm=TRUE)
summary(meanTemp)

#plot
meantempPlot <- ggplot(meanTemp, aes(x = tp, y = mtemp, fill = target_treat)) +
  geom_point(aes(color = target_treat), size = 5, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=mtemp-se,ymax=mtemp+se), lwd=0.4, width=0.2)+
  ggtitle("Mean Temperature")+
  xlab("Time")+
  ylab("Temperature (°C)")+
  ylim(27,32)+
  facet_wrap(target_treat~., ncol=2)+
  theme(legend.position = "none")
meantempPlot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs_newdata/MeanTemp_row.pdf", meantempPlot, width=6, height=8, units=c("in"), useDingbats=FALSE)


#### PCO2 ####

# plot average pco2 by time point and treatment
#plot 3 tanks separately:
pco2Plot <- ggplot(calc_data, aes(x = tp, y = mpco2, fill = target_treat, shape = tank)) +
  geom_point(aes(color = target_treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("pCO2")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(target_treat~., ncol=2, scales = "free")
pco2Plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs_newdata/pCO2_row.pdf", pco2Plot, width=6, height=8, units=c("in"), useDingbats=FALSE)

# filter out the treatments not considered here
pCO2_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

# summarize data to create cumulative treatments
T0_pco2 = pCO2_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,mpco2)

meanpco2_T0 <- summarySE(T0_pco2, measurevar="mpco2", groupvars=c("target_treat"), na.rm=TRUE)

lmT0_pco2=aov(mpco2~temp*pco2, data=T0_pco2)
summary(lmT0_pco2)
par(mfrow=c(2,2))
plot(lmT0_pco2)

##weirdly enough these data are normal at T0
TukeyHSD(lmT0_pco2)
##relevant things not significantly different at T0: 31:700-31:400 AND 31:700-28:280

T30_pco2 <- pCO2_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(tp,target_treat,tank,temp,pco2,mpco2)

meanpco2_T30 <- summarySE(T30_pco2, measurevar="mpco2", groupvars=c("target_treat"), na.rm=TRUE)

lmT30_pco2=aov(mpco2~temp*pco2, data=T30_pco2)
summary(lmT30_pco2)
par(mfrow=c(2,2))
plot(lmT30_pco2)
###looks weird, need to transform

lmT30_pco2_log=aov(log(mpco2)~temp*pco2, data=T30_pco2)
summary(lmT30_pco2_log)
par(mfrow=c(2,2))
plot(lmT30_pco2_log)

TukeyHSD(lmT30_pco2_log)
### = 31:700-28:400, 28:700-28:400, 31:700-28:280 not sig different

T60_pco2 <- pCO2_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,mpco2)

meanpco2_T60 <- summarySE(T60_pco2, measurevar="mpco2", groupvars=c("target_treat"), na.rm=TRUE)

lmT60_pco2=aov(mpco2~temp*pco2, data=T60_pco2)
summary(lmT60_pco2)
par(mfrow=c(2,2))
plot(lmT60_pco2)
###looks weird, need to transform

lmT60_pco2_log=aov(log(mpco2)~temp*pco2, data=T60_pco2)
summary(lmT60_pco2_log)
par(mfrow=c(2,2))
plot(lmT60_pco2_log)
#looks good

TukeyHSD(lmT60_pco2_log)
### all of these look good when the data is included from T30-T60 and from T0-T60

T90_pco2 <- pCO2_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,mpco2)

meanpco2_T90 <- summarySE(T90_pco2, measurevar="mpco2", groupvars=c("target_treat"), na.rm=TRUE)

lmT90_pco2=aov(mpco2~temp*pco2, data=T90_pco2)
summary(lmT90_pco2)
par(mfrow=c(2,2))
plot(lmT90_pco2)
###looks a little weird, need to transform

lmT90_pco2_log=aov(log(mpco2)~temp*pco2, data=T90_pco2)
summary(lmT90_pco2_log)
par(mfrow=c(2,2))
plot(lmT90_pco2_log)
#looks better

TukeyHSD(lmT90_pco2_log)
### All are different with all data T0-T90 included and T30-T90


# plot boxplots of cumulative treatments
pco2Plot <- ggplot(T90_pco2, aes(x = target_treat, y = mpco2)) +
  geom_boxplot(fill="gray", outlier.shape=NA)+
  geom_jitter(shape=16, position=position_jitter(0.1), size=1.5)+ 
  xlab("Treatment")+
  ylab("pCO2")+
  scale_y_continuous(breaks=seq(200,3800, by = 500), limits=c(200,3800))+
  theme(legend.position = "none")+
  theme_bw()
pco2Plot
#save as 4"x4" pdf
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs_newdata/pCO2/T90_pCO2.pdf", pco2Plot, width=4, height=4, units=c("in"), useDingbats=FALSE)


# Dates corresponding to time points in the experiment:
# T0 = August 10, 2015 (calc_data$msmt = 1)
# T30 = September 19, 2015 (calc_data$msmt = 3)
# T60 = October 9 2015 (calc_data$msmt = 5)
# T90 = November 10, 2015 (calc_data$msmt = 8)



#### MEASURED PH ####
# plot measured average pH by time point and treatment
#plot 3 tanks separately:
ggplot(msd_data, aes(x = tp, y = pHm, fill = target_treat)) +
  geom_point(aes(color = target_treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Measured Mean pH")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(target_treat~., nrow=4, ncol=4)

# filter out the treatments not considered here
pHm_data = msd_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_mpH <- pHm_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,pHm)

meanpHm_T0 <- summarySE(T0_mpH, measurevar="pHm", groupvars=c("target_treat"), na.rm=TRUE)

T30_mpH <- pHm_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,pHm)

meanpHm_T30 <- summarySE(T30_mpH, measurevar="pHm", groupvars=c("target_treat"), na.rm=TRUE)

T60_mpH <- pHm_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,pHm)

meanpHm_T60 <- summarySE(T60_mpH, measurevar="pHm", groupvars=c("target_treat"), na.rm=TRUE)

T90_mpH <- pHm_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,pHm)

meanpHm_T90 <- summarySE(T90_mpH, measurevar="pHm", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meanpHplot <- ggplot(meanpH, aes(x = tp, y = pHm, fill = treat)) +
  geom_point(aes(color = treat), size = 5, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= pHm-se,ymax= pHm+se), lwd=0.4,width=0.3)+
  ggtitle("Measured Mean pH")+
  facet_wrap(treat~., nrow=1, ncol=6)
meanpHplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeasuredMeanpH_row.pdf", meanpHplot, width=12, height=6, units=c("in"), useDingbats=FALSE)

#### CALCULATED PH ####
# plot calculated average pH by time point and treatment
#plot 3 tanks separately:
ggplot(data, aes(x = tp, y = pHc, fill = treat)) +
  geom_point(aes(color = treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Calculated Mean pH")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(treat~., nrow=4, ncol=4)

# filter out the treatments not considered here
pHc_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_cpH <- pHc_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,pHc)

meanpHc_T0 <- summarySE(T0_cpH, measurevar="pHc", groupvars=c("target_treat"), na.rm=TRUE)

T30_cpH <- pHc_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,pHc)

meanpHc_T30 <- summarySE(T30_cpH, measurevar="pHc", groupvars=c("target_treat"), na.rm=TRUE)

T60_cpH <- pHc_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,pHc)

meanpHc_T60 <- summarySE(T60_cpH, measurevar="pHc", groupvars=c("target_treat"), na.rm=TRUE)

T90_cpH <- pHc_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,pHc)

meanpHc_T90 <- summarySE(T90_cpH, measurevar="pHc", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meanpHplot <- ggplot(meanpHc, aes(x = tp, y = pHc, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= pHc-se,ymax= pHc+se), lwd=0.4,width=0.3)+
  ggtitle("Calculated Mean pH")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanpHplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/CalculatedMeanpH_row.pdf", meanpHplot, width=12, height=6, units=c("in"), useDingbats=FALSE)

#### SALINITY ####
# plot average salinity by time point and treatment
#plot 3 tanks separately:
ggplot(data, aes(x = tp, y = sal, fill = treat)) +
  geom_point(aes(color = treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Mean pH")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(treat~., nrow=4, ncol=4)

# filter out the treatments not considered here
sal_data = msd_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_sal <- sal_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,sal)

meansal_T0 <- summarySE(T0_sal, measurevar="sal", groupvars=c("target_treat"), na.rm=TRUE)

T30_sal <- sal_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,sal)

meansal_T30 <- summarySE(T30_sal, measurevar="sal", groupvars=c("target_treat"), na.rm=TRUE)

T60_sal <- sal_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,sal)

meansal_T60 <- summarySE(T60_sal, measurevar="sal", groupvars=c("target_treat"), na.rm=TRUE)

T90_sal <- sal_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,sal)

meansal_T90 <- summarySE(T90_sal, measurevar="sal", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meansalplot <- ggplot(meansal, aes(x = tp, y = sal, fill = treat)) +
  geom_point(aes(color = treat), size = 5, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= sal-se,ymax= sal+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Salinity")+
  xlab("Time")+
  scale_y_continuous(name ="Salinity", limits = c(28,34), breaks = seq(28,34,2))+
  facet_wrap(treat~., nrow=1, ncol=6)+
  theme(legend.position="none")
meansalplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanSalinity_row.pdf", meansalplot, width=13, height=6, units=c("in"), useDingbats=FALSE)

#### ALKALINITY ####

# filter out the treatments not considered here
alk_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_alk <- alk_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,alk)

meanalk_T0 <- summarySE(T0_alk, measurevar="alk", groupvars=c("target_treat"), na.rm=TRUE)

T30_alk <- alk_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,alk)

meanalk_T30 <- summarySE(T30_alk, measurevar="alk", groupvars=c("target_treat"), na.rm=TRUE)

T60_alk <- alk_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,alk)

meanalk_T60 <- summarySE(T60_alk, measurevar="alk", groupvars=c("target_treat"), na.rm=TRUE)

T90_alk <- alk_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,alk)

meanalk_T90 <- summarySE(T90_alk, measurevar="alk", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meanalkplot <- ggplot(meanalk, aes(x = tp, y = alk, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= alk-se,ymax= alk+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Alkalinity")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanalkplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanAlkalinity_row.pdf", meanalkplot, width=13, height=6, units=c("in"), useDingbats=FALSE)

#### DIC ####

# filter out the treatments not considered here
dic_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_dic <- dic_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,dic)

meandic_T0 <- summarySE(T0_dic, measurevar="dic", groupvars=c("target_treat"), na.rm=TRUE)

T30_dic <- dic_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,dic)

meandic_T30 <- summarySE(T30_dic, measurevar="dic", groupvars=c("target_treat"), na.rm=TRUE)

T60_dic <- dic_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,dic)

meandic_T60 <- summarySE(T60_dic, measurevar="dic", groupvars=c("target_treat"), na.rm=TRUE)

T90_dic <- dic_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,dic)

meandic_T90 <- summarySE(T90_dic, measurevar="dic", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meandicplot <- ggplot(meandic, aes(x = tp, y = dic, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= dic-se,ymax= dic+se), lwd=0.4,width=0.3)+
  ggtitle("Mean DIC")+
  facet_wrap(treat~., nrow=1, ncol=8)
meandicplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanDIC_row.pdf", meandicplot, width=13, height=6, units=c("in"), useDingbats=FALSE)

#### Bicarbonate (HCO3) ####

# filter out the treatments not considered here
hco3_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_hco3 <- hco3_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,hco3)

meanhco3_T0 <- summarySE(T0_hco3, measurevar="hco3", groupvars=c("target_treat"), na.rm=TRUE)

T30_hco3 <- hco3_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,hco3)

meanhco3_T30 <- summarySE(T30_hco3, measurevar="hco3", groupvars=c("target_treat"), na.rm=TRUE)

T60_hco3 <- hco3_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,hco3)

meanhco3_T60 <- summarySE(T60_hco3, measurevar="hco3", groupvars=c("target_treat"), na.rm=TRUE)

T90_hco3 <- hco3_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,hco3)

meanhco3_T90 <- summarySE(T90_hco3, measurevar="hco3", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meanhco3plot <- ggplot(meanhco3, aes(x = tp, y = hco3, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= hco3-se,ymax= hco3+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Bicarbonate (HCO3)")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanhco3plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanHCO3_row.pdf", meanhco3plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### Carbon Dioxide (CO2) ####

# filter out the treatments not considered here
co2_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_co2 <- co2_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,co2)

meanco2_T0 <- summarySE(T0_co2, measurevar="co2", groupvars=c("target_treat"), na.rm=TRUE)

T30_co2 <- co2_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,co2)

meanco2_T30 <- summarySE(T30_co2, measurevar="co2", groupvars=c("target_treat"), na.rm=TRUE)

T60_co2 <- co2_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,co2)

meanco2_T60 <- summarySE(T60_co2, measurevar="co2", groupvars=c("target_treat"), na.rm=TRUE)

T90_co2 <- co2_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,co2)

meanco2_T90 <- summarySE(T90_co2, measurevar="co2", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meanco2plot <- ggplot(meanco2, aes(x = tp, y = co2, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= co2-se,ymax= co2+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Carbon Dioxide (CO2)")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanco2plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanCO2_row.pdf", meanco2plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### Carbonate Ion (CO3) ####

# filter out the treatments not considered here
co3_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_co3 <- co3_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,co3)

meanco3_T0 <- summarySE(T0_co3, measurevar="co3", groupvars=c("target_treat"), na.rm=TRUE)

T30_co3 <- co3_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,co3)

meanco3_T30 <- summarySE(T30_co3, measurevar="co3", groupvars=c("target_treat"), na.rm=TRUE)

T60_co3 <- co3_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,co3)

meanco3_T60 <- summarySE(T60_co3, measurevar="co3", groupvars=c("target_treat"), na.rm=TRUE)

T90_co3 <- co3_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,co3)

meanco3_T90 <- summarySE(T90_co3, measurevar="co3", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meanco3plot <- ggplot(meanco3, aes(x = tp, y = co3, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= co3-se,ymax= co3+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Carbonate (CO3)")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanco3plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanCO3_row.pdf", meanco3plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### Aragonite Saturation State (sat) ####

# filter out the treatments not considered here
sat_data = calc_data %>%
  filter(target_treat!="280_31") %>%
  filter(target_treat!="400_28")

#SummarySE by time point
T0_sat <- sat_data %>%
  filter(tp == "T0") %>%
  select(target_treat,tank,temp,pco2,sat)

meansat_T0 <- summarySE(T0_sat, measurevar="sat", groupvars=c("target_treat"), na.rm=TRUE)

T30_sat <- sat_data %>%
  filter(tp == "T0" | tp == "T30") %>%
  select(target_treat,tank,temp,pco2,sat)

meansat_T30 <- summarySE(T30_sat, measurevar="sat", groupvars=c("target_treat"), na.rm=TRUE)

T60_sat <- sat_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60") %>%
  select(target_treat,tank,temp,pco2,sat)

meansat_T60 <- summarySE(T60_sat, measurevar="sat", groupvars=c("target_treat"), na.rm=TRUE)

T90_sat <- sat_data %>%
  filter(tp == "T0" | tp == "T30" | tp == "T60" | tp == "T90") %>%
  select(target_treat,tank,temp,pco2,sat)

meansat_T90 <- summarySE(T90_sat, measurevar="sat", groupvars=c("target_treat"), na.rm=TRUE)

#plot
meansatplot <- ggplot(meansat, aes(x = tp, y = sat, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= sat-se,ymax= sat +se), lwd=0.4,width=0.3)+
  ggtitle("Mean Saturation State")+
  facet_wrap(treat~., nrow=1, ncol=8)
meansatplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanSAT_row.pdf", meansatplot, width=13, height=6, units=c("in"), useDingbats=FALSE)
