#Written by Hannah Aichelman
#Contact: hannahaichelman@gmail.com
#Last updated: 4/22/2020

#This script is used to plot water chemistry data associated with the time course physiology manuscript
#as well as summarize the parameters through time, which is reported in supplementary tables 1 and 2.

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

#read in data
data <- read.csv("compiledwaterchem.csv")
head(data)

#reorganize levels of treatment factor
levels(data$treat)
data$treat = factor(data$treat,levels = c("358_28", "325_31", "471_28", "424_31", "674_28", "606_31", "2750_28", "2917_31"))
levels(data$treat) <- c("ambient pCO2, 28°C","Xlow pCO2, 31°C","Xlow pCO2, 28°C","ambient pCO2, 31°C","next century pCO2, 28°C","next century pCO2, 31°C","extreme pCO2, 28°C","extreme pCO2, 31°C")

#filter out the treatments not considered here
data2 <- data %>%
  filter(treat!="Xlow pCO2, 31°C") %>%
  filter(treat!="Xlow pCO2, 28°C") %>%
  drop_na()
summary(data2)

levels(data2$treat)

#### WATER TEMP ####
# plot average water temperature by time point and treatment
#plot 3 tanks separately:
quartz()
ggplot(data2, aes(x = tp, y = mtemp, fill = treat)) +
  geom_point(aes(color = treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Mean Temperature")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(treat~., nrow=4, ncol=4)

#SummarySE, average temperature by time point
#Throughout, use SummarySE to get data summaries that were used in supplementary tables
meanTemp <- summarySE(data2, measurevar="mtemp", groupvars=c("treat","tp"), na.rm=TRUE)
summary(meanTemp)

#plot
meantempPlot <- ggplot(meanTemp, aes(x = tp, y = mtemp, fill = treat)) +
  geom_point(aes(color = treat), size = 5, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=mtemp-se,ymax=mtemp+se), lwd=0.4, width=0.3)+
  ggtitle("Mean Temperature")+
  xlab("Time")+
  ylab("Temperature (°C)")+
  ylim(27,32)+
  facet_grid(treat~., nrow=1, ncol=6)+
  theme(legend.position = "none")
meantempPlot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanTemp_row.pdf", meantempPlot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### PCO2 ####
# plot average pco2 by time point and treatment
#plot 3 tanks separately:
ggplot(data, aes(x = tp, y = mpco2, fill = treat)) +
  geom_point(aes(color = treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Mean Temperature")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(treat~., nrow=4, ncol=4)

#SummarySE
meanpco2 <- summarySE(data2, measurevar="mpco2", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanpco2)

#plot
meanpco2Plot <- ggplot(meanpco2, aes(x = tp, y = mpco2, fill = treat)) +
  geom_point(aes(color = treat), size = 5, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= mpco2-se,ymax= mpco2+se), lwd=0.4,width=0.3)+
  ggtitle("Mean pCO2")+
  xlab("Time")+
  scale_y_continuous(name="pCO2 (uatm)", breaks = seq(0,3500,500))+
  #ylim()+
  facet_wrap(treat~., nrow=1, ncol=6)+
  theme(legend.position="none")
meanpco2Plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanpCO2_row.pdf", meanpco2Plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### MEASURED PH ####
# plot measured average pH by time point and treatment
#plot 3 tanks separately:
ggplot(data, aes(x = tp, y = pHm, fill = treat)) +
  geom_point(aes(color = treat), size = 2, position = position_dodge(width=0.75)) +
  ggtitle("Measured Mean pH")+
  #scale_fill_manual(values=c("F"="lightskyblue1", "N"="firebrick2"))+
  facet_wrap(treat~., nrow=4, ncol=4)

#SummarySE
meanpH <- summarySE(data2, measurevar="pHm", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanpH)

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

#SummarySE
meanpHc <- summarySE(data2, measurevar="pHc", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanpHc)
#write.csv(meanpHc, "~/Documents/BU/NEU_experiment/waterchem/MeanpHc.csv")

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

#SummarySE
meansal <- summarySE(data2, measurevar="sal", groupvars=c("treat","tp"), na.rm=TRUE)
head(meansal)

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

#SummarySE
meanalk <- summarySE(data2, measurevar="alk", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanalk)
#write.csv(meanalk, "~/Documents/BU/NEU_experiment/waterchem/Meanalk.csv")

#plot
meanalkplot <- ggplot(meanalk, aes(x = tp, y = alk, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= alk-se,ymax= alk+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Alkalinity")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanalkplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanAlkalinity_row.pdf", meanalkplot, width=13, height=6, units=c("in"), useDingbats=FALSE)

#### DIC ####

#SummarySE
meandic <- summarySE(data2, measurevar="dic", groupvars=c("treat","tp"), na.rm=TRUE)
head(meandic)
#write.csv(meandic, "~/Documents/BU/NEU_experiment/waterchem/Meandic.csv")

#plot
meandicplot <- ggplot(meandic, aes(x = tp, y = dic, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= dic-se,ymax= dic+se), lwd=0.4,width=0.3)+
  ggtitle("Mean DIC")+
  facet_wrap(treat~., nrow=1, ncol=8)
meandicplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanDIC_row.pdf", meandicplot, width=13, height=6, units=c("in"), useDingbats=FALSE)

#### Bicarbonate (HCO3) ####

#SummarySE
meanhco3 <- summarySE(data2, measurevar="hco3", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanhco3)
write.csv(meanhco3, "~/Documents/BU/NEU_experiment/waterchem/MeanHCO3.csv")

#plot
meanhco3plot <- ggplot(meanhco3, aes(x = tp, y = hco3, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= hco3-se,ymax= hco3+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Bicarbonate (HCO3)")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanhco3plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanHCO3_row.pdf", meanhco3plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### Carbon Dioxide (CO2) ####

#SummarySE
meanco2 <- summarySE(data2, measurevar="co2", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanco2)
#write.csv(meanco2, "~/Documents/BU/NEU_experiment/waterchem/MeanCO2.csv")

#plot
meanco2plot <- ggplot(meanco2, aes(x = tp, y = co2, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= co2-se,ymax= co2+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Carbon Dioxide (CO2)")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanco2plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanCO2_row.pdf", meanco2plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### Carbonate Ion (CO3) ####

#SummarySE
meanco3 <- summarySE(data2, measurevar="co3", groupvars=c("treat","tp"), na.rm=TRUE)
head(meanco3)
#write.csv(meanco3, "~/Documents/BU/NEU_experiment/waterchem/MeanCO3.csv")

#plot
meanco3plot <- ggplot(meanco3, aes(x = tp, y = co3, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= co3-se,ymax= co3+se), lwd=0.4,width=0.3)+
  ggtitle("Mean Carbonate (CO3)")+
  facet_wrap(treat~., nrow=1, ncol=8)
meanco3plot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanCO3_row.pdf", meanco3plot, width=13, height=6, units=c("in"), useDingbats=FALSE)


#### Aragonite Saturation State (sat) ####

#SummarySE
meansat <- summarySE(data2, measurevar="sat", groupvars=c("treat","tp"), na.rm=TRUE)
head(meansat)
#write.csv(meansat, "~/Documents/BU/NEU_experiment/waterchem/Meansat.csv")

#plot
meansatplot <- ggplot(meansat, aes(x = tp, y = sat, fill = treat)) +
  geom_point(aes(color = treat), size = 4, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin= sat-se,ymax= sat +se), lwd=0.4,width=0.3)+
  ggtitle("Mean Saturation State")+
  facet_wrap(treat~., nrow=1, ncol=8)
meansatplot
ggsave(file="~/Documents/BU/NEU_experiment/waterchem/figs/MeanSAT_row.pdf", meansatplot, width=13, height=6, units=c("in"), useDingbats=FALSE)
