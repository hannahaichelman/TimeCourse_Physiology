#Written by Hannah Aichelman
#Contact: hannahaichelman@gmail.com
#Last updated: 2/11/2021

#This script is used to analyze changes in coral physiology data through time and is associated with Aichelman et al. (2020) 
#The dataset includes host physiology (calcification rate, total protein, and total carbohydrate) 
#as well as symbiont physiology (symbiont density and chlorophyll a content) from two species of 
#Caribbean reef-building corals (Siderastrea siderea and Pseudodiploria strigosa) across reef zone (forereef and nearshore).
#These corals were exposed to temperature (28 and 31C) and acidification (present day, next century, and extreme pco2 conditions)
#for 95 days, and physiology measurements were taken every 30 days. 
#Throughout, SSID = Siderastrea siderea and PSTR = Pseudodiploria strigosa


#set working directory
setwd("~/Documents/BU/NEU_experiment/data_sheets")

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
library("tidyverse") #easily select physdata to plot
library("plotly") #make cool interactive plots
library("cowplot") #makes your ggplots pretty and simple
library("lme4") #for linear models/stats
library("car")
library("multcomp") #for glhtest and post-host analysis of lsmeans
library("lmerTest")
library("corrplot")
library("MuMIn")

#### Read in data sheets ####
calcdata <- read.csv("CalculatedGrowthRates_final.csv")
physdata <- read.csv("NEU_cor_nub_2018_cleanedforR.csv")

# add Carbohydrate Data in to physiology data
carbdata <- read.csv("NEU_Carbohydrates_Final.csv")
str(carbdata)

carb_sum = carbdata %>%
  group_by(id) %>%
  summarise(carbs_mg_mL = mean(calculated_carbs_mg_mL))
str(carb_sum)
head(carb_sum)

# write.csv(carb_sum, "carb_summary.csv") #used this to look at high standard deviation measurements

#combine carb data with remaining physiology, join by coral ID
physdata <- left_join(physdata, carb_sum, by="id")
head(physdata)

#correct carbohydrate data for surface area of fragments
physdata <- physdata %>%
  mutate(carbs_mg_cm2 = carbs_mg_mL/sur_area)


#### Looking for outliers ####

# Looking for outliers with calcification data first
# Notes include references to photos and datasheets to confirm if the outlier was reasonable or not

# SIDERASTREA SIDEREA
# > which.min(calcdataSid$Growthmg.cm2.day)
# [1] 105
# > calcdataSid[105,]
        # id     tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 147 SNSY27 400_31_2  400   31   S  N   Y 9.945446  T30        -10.02831    400_31

#removing this ID after looking at original calculation spreadsheet, T0 and T30 weights were very different and there is no T0 bw measurement on Colleen's original BW spreadsheet...maybe weight was added in mistakenly somewhere? Regardless, taking it out.

# > which(calcdataSid$Growthmg.cm2.day < -1)
# [1] 189 206 222
# > calcdataSid[189,]
        # id      tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 460 SNSX30 2800_28_1 2800   28   S  N   X 9.014287  T60           -1.182   2800_28 
# DON'T SEE ANY REASON ON DAVIES_CALCRATES SPREADSHEET TO TOSS OUT, but weird because looks healthy in photos at T60 and T90

# > calcdataSid[206,]
        # id      tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 494 SNSY31 2800_28_2 2800   28   S  N   Y 7.565368  T60           -1.057   2800_28
# DON'T SEE ANY REASON ON DAVIES_CALCRATES SPREADSHEET TO TOSS OUT, again strange because it looks healthy in the T60 photo

# > calcdataSid[222,]
        # id      tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 526 SNSZ42 2800_31_3 2800   31   S  N   Z 8.058396  T60           -1.417   2800_31
# DON'T SEE ANY REASON ON DAVIES_CALCRATES SPREADSHEET TO TOSS OUT, again strange because it looks healthy in the T60 photo

# > which(calcdataSid$Growthmg.cm2.day > 3)
# [1] 71 85
# > calcdataSid[71,]
       # id     tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 106 SNSX2 280_31_1  280   31   S  N   X 6.087735  T30         3.811296    280_31
# DON'T SEE ANY REASON ON DAVIES_CALCRATES SPREADSHEET OR PHOTOS TO TOSS OUT
	
# > calcdataSid[85,]
       # id     tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 123 SNSX4 280_31_1  280   31   S  N   X  9.53357  T30         3.628373    280_31
# DON'T SEE ANY REASON ON DAVIES_CALCRATES SPREADSHEET OR PHOTOS TO TOSS OUT

# PSEUDODIPLORIA STRIGOSA
# > which(calcdataDip$Growthmg.cm2.day < -2.5)
# [1] 174 176 281
# > calcdataDip[174,]
        # id     tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 676 SNPZ39 280_31_2  280   31   P  N   Z 8.038008  T60            -4.27    280_31
# SNPZ39 is very bleached, looks like it loses tissue between T30 and T60 (no T0 weight or photo for this dude either). But even if dead I don't know that this is a reasonable rate...
#remove this individual

# > calcdataDip[176,]
        # id     tank pco2 temp spp rz gen sur_area time Growthmg.cm2.day treatment
# 678 SNPZ43 400_28_3  400   28   P  N   Z 7.014957  T60           -4.319    400_28

# > calcdataDip[281,]
         # id     tank pco2 temp spp rz gen sur_area   time Growthmg.cm2.day treatment
# 1362 SNPZ43 400_28_3  400   28   P  N   Z 7.014957 T0_T60           -3.468    400_28
# SNP43 is very bleached and looks like it is losing tissue at T30 and T60, so makes sense that it has very negative growth


#### Data Analysis ####

#### Calcification ####
head(calcdata)
dim(calcdata)
str(calcdata)

#create combined factors for treatment (pco2 and temperature) and genotype/rz, both will be used for stats
calcdata <-cbind(calcdata[0:10],"treatment"=paste(calcdata$pco2, calcdata$temp, sep="_"))
calcdata <-cbind(calcdata[0:11],"gen_rz"=paste(calcdata$gen, calcdata$rz, sep="_"))

#make pco2 and temperature and gen_rz factors
calcdata$pco2 <- as.factor(calcdata$pco2)
calcdata$temp <- as.factor(calcdata$temp)
calcdata$gen_rz <- as.factor(calcdata$gen_rz)
calcdata$rz <- as.factor(calcdata$rz)

#check treatment names
calcdata$treatment = factor(calcdata$treatment,levels = c("280_28", "280_31", "400_28", "400_31", "700_28", "700_31", "2800_28", "2800_31"))
levels(calcdata$treatment)

#rename levels based on actual averages of pco2
levels(calcdata$treatment) <- c("298uatm, 28.0°C","325uatm, 31°C","471uatm, 28°C","388uatm, 31.1°C","663uatm, 28.0°C","662uatm, 31.0°C","2973uatm, 28.1°C","3245uatm, 30.7°C")

#re-name pco2 levels for stats
levels(calcdata$pco2)
#re-name levels of factors by category instead of number
levels(calcdata$pco2) <- c("ambient","ambient","next century","extreme")

#check timepoint names
calcdata$time = factor(calcdata$time,levels = c("T30", "T60", "T90", "T0_T60", "T0_T90"))
levels(calcdata$time)


### Create sub-sets of physdata to plot
# the filters can be changed by time to look at either time course or one specific time point
# the treatment filters reflect the decision to remove treatments that had problems in the experiment

# SSID 
calcdataSid <- calcdata %>%
	select(id, rz, spp, gen, gen_rz, time, pco2, temp, treatment, Growthmg.cm2.day) %>%
	filter(spp=="S") %>%
	filter(time=="T30" | time=="T60" | time=="T90") %>%
	filter(treatment!="325uatm, 31°C") %>%
	filter(treatment!="471uatm, 28°C") %>%
	filter(id != "SNSY27" ) %>% #see above for outlier analysis
	drop_na()
summary(calcdataSid)

#PSTR 
calcdataDip <- calcdata %>%
	select(id, rz, spp, gen, gen_rz, time, pco2, temp, treatment, Growthmg.cm2.day) %>%
	filter(spp=="P")%>%
  filter(time=="T30" | time=="T60" | time=="T90") %>%
  filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
	filter(id != "SNPZ39")%>% #see above for outlier analysis
	drop_na()
summary(calcdataDip)

## check normality and homogeneity of variance
# PSTR
shapiro.test(calcdataDip$Growthmg.cm2.day) #Dips are normal
# data:  calcdataDip$Growthmg.cm2.day
# W = 0.98575, p-value = 0.07569

# qq plot
ggqqplot(calcdataDip$Growthmg.cm2.day, main = "qq Plot PSTR Growth")
# density plot
ggdensity(calcdataDip$Growthmg.cm2.day, main = "Density Plot PSTR Growth")

# SSID
shapiro.test(calcdataSid$Growthmg.cm2.day+2) #Sids are not normal
# data:  calcdataSid$Growthmg.cm2.day
# W = 0.97634, p-value = 0.001618
shapiro.test(log(calcdataSid$Growthmg.cm2.day+2)) #Sids are still not normal
# data:  log(calcdataSid$Growthmg.cm2.day + 2)
# W = 0.85095, p-value = 3.399e-13
shapiro.test(sqrt(calcdataSid$Growthmg.cm2.day+2)) #Sids are still not normal
# data:  sqrt(calcdataSid$Growthmg.cm2.day + 2)
# W = 0.93299, p-value = 4.589e-08
sidcalc_cube <- calcdataSid$Growthmg.cm2.day^(1/3) #Sids are still not normal
shapiro.test(sidcalc_cube)
# data:  sidcalc_cube
# W = 0.93374, p-value = 1.635e-07

#qq plot
ggqqplot(calcdataSid$Growthmg.cm2.day, main = "qq Plot SSID Growth")

#density plot
ggdensity(calcdataSid$Growthmg.cm2.day, main = "Density Plot SSID Growth - sqrt")

#none of the transformations work for SSID calcification, so going to have to use gamlss for models (see below)

#### Calcification Model Selection ####
#large datasets do well for backwards model selection, but small datasets works well for forward model selection
#don't include interactons that you dont have hypotheses for, or that you don't have as main effects 

#Using a forward model selection method, following Dr. Chris Schmitt's (BU) notes from here:
#https://fuzzyatelin.github.io/bioanth-stats/module-16/module-16.html

#SSID model selection using gamlss (this section was completed with help from Sara Smith-Wuitchik)
library(goft)
library(fitdistrplus)
library(gamlss)

sidcalcmodeldata <- calcdataSid %>%
  mutate(newgrowth = Growthmg.cm2.day + 1.5)  #add 1.5 to create all positive numbers for distribution fit

# Make a Cullen and Frey plot to help decide which distributions to test
descdist(sidcalcmodeldata$newgrowth, discrete = F) # lognormal or Weibull or gamma
fit.weibull <- fitdist(sidcalcmodeldata$newgrowth, "weibull")
fit.lognorm <- fitdist(sidcalcmodeldata$newgrowth, "lnorm")
fit.gamma <- fitdist(sidcalcmodeldata$newgrowth, "gamma")
# check fit of model visually
plot(fit.weibull) # best fit of the three
plot(fit.lognorm)
plot(fit.gamma)
# check AIC to confirm visual fit
fit.weibull$aic # lowest aic value
fit.lognorm$aic
fit.gamma$aic


# model selection for gamlss
# write out whatever your full model is with all possible variable/interactions that may be relevant
m1 <- gamlss(newgrowth ~ time_relevel*pco2_relevel*temp_relevel*rz + random(gen_rz), family = WEI(), method = RS(), data = sidcalcmodeldata)
step.full <- stepGAIC(m1, direction = c("forward"), trace = T)
# get final model
formula(step.full)
# get pvalues, check interactions, pull estimates
step.full$anova
# get estimates, std error, p-values, interactions
summary(step.full)
# visualize
plot(step.full)

#re-order factor levels for SSID calcification to get other comparisons output from gamlss
#then re-run above code
levels(sidcalcmodeldata$time_relevel)
sidcalcmodeldata$time_relevel = factor(sidcalcmodeldata$time,levels = c("T90","T60","T30", "T0_T60","T0_T90"))

levels(sidcalcmodeldata$pco2)
sidcalcmodeldata$pco2_relevel = factor(sidcalcmodeldata$pco2,levels = c("extreme", "next century", "ambient"))

levels(sidcalcmodeldata$temp)
sidcalcmodeldata$temp_relevel = factor(sidcalcmodeldata$temp,levels = c("31", "28"))

#PSTR null model 
m0 <- lmer(Growthmg.cm2.day ~ 1 + (1|gen_rz), data=calcdataDip)
summary(m0)
AIC(m0)
#[1] 435.6938

#test adding parameters and check AIC
m1 <- lmer(Growthmg.cm2.day ~ time+temp+pco2+rz + time:temp + temp:pco2 + temp:rz + (1|gen_rz), data=calcdataDip)
summary(m1)
anova(m1)
AIC(m1)
#[1] 409.7194 - time + (1|gen_rz)
#[1] 394.4723 - time + temp + (1|gen_rz)
#[1] 390.0929 - time + temp + pco2 + (1|gen_rz)
#[1] 387.669 - time + temp + pco2 + rz +(1|gen_rz)
#[1] 372.3691 - time+temp+pco2+rz + time:temp + (1|gen_rz)
#[1] 364.5652 - time+temp+pco2+rz + time:temp + temp:pco2 + (1|gen_rz)
#[1] 363.3046 - time+temp+pco2+rz + time:temp + temp:pco2 + temp:rz + (1|gen_rz)

#post hoc tests using Tukey-adjusted comparisons
lsmeans(m1, pairwise~temp:rz, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "PSTR Calcification Residuals")

#to look at pairwise of all factors, the output from lsmeans is easier to interpret
m1 <- lmer(Growthmg.cm2.day ~ time*temp*pco2*rz + (1|gen_rz), data=calcdataDip)
summary(glht(m1, lsm(pairwise ~ rz|time|pco2|temp)))

#summarySE to group data for figures 
Sidgrowth.means <- summarySE(calcdataSid, measurevar="Growthmg.cm2.day", groupvars=c("rz","time","treatment"))
summary(Sidgrowth.means)
Dipgrowth.means <- summarySE(calcdataDip, measurevar="Growthmg.cm2.day", groupvars=c("rz","time","treatment"))
summary(Dipgrowth.means)

#SSID calcification plot  
sidgrowth <- ggplot(Sidgrowth.means, aes(x = time, y = Growthmg.cm2.day, fill = rz)) +
  geom_errorbar(aes(x = time, ymax = Growthmg.cm2.day+se, ymin = Growthmg.cm2.day-se), width = .2, position = position_dodge(width=0.2)) +
  geom_point(aes(color=rz), position = position_dodge(width = 0.2), size=3)+
  #ggtitle("S. siderea - Calcification Rate")+
  ylab(bquote("Calcification rate (mg" ~cm^-2~~day^-1~")"))+
  xlab("Time")+
  ylim(-1,2.5)+
  theme_cowplot()+
  scale_color_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  scale_fill_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  facet_wrap(treatment~., nrow=3, ncol=2)+
  #geom_smooth(aes(group=rz, color=rz, fill=rz), position = position_dodge(width=0.2))+
  stat_smooth(method="lm",se=FALSE, formula = y~poly(x,2),size=1, aes(group=rz, color=rz))+
  geom_hline(yintercept=0, linetype="dashed", color="gray", size=1)
  #theme(legend.position="right", legend.title=)
  #scale_fill_discrete(name = "Reef Zone", labels = c("FR", "NS"))
sidgrowth
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/SSID_calc_final.pdf", sidgrowth, width=6, height=5, units=c("in"), useDingbats=FALSE)

#PSTR calcification plot
dipgrowth <- ggplot(Dipgrowth.means, aes(x = time, y = Growthmg.cm2.day, fill = rz)) +
  geom_errorbar(aes(x = time, ymax = Growthmg.cm2.day+se, ymin = Growthmg.cm2.day-se), width = .2, position = position_dodge(width=0.2)) +
  geom_point(aes(color=rz), position = position_dodge(width = 0.2), size=3)+
  #ggtitle("P. strigosa - Calcification Rate")+
  ylab(bquote("Calcification rate (mg" ~cm^-2~~day^-1~")"))+
  xlab("Time")+
  ylim(-1,2.5)+
  theme_cowplot()+
  scale_color_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  scale_fill_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  facet_wrap(treatment~., nrow=3, ncol=2)+
  #geom_smooth(aes(group=rz, color=rz, fill=rz), position = position_dodge(width=0.2))+
  stat_smooth(method="lm",se=FALSE, formula = y~poly(x,2),size=1, aes(group=rz, color=rz))+
  geom_hline(yintercept=0, linetype="dashed", color="gray", size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
dipgrowth
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/PSTR_calc_final.pdf", dipgrowth, width=6, height=5, units=c("in"), useDingbats=FALSE)


#### Remaining host and symbiont physiology data ####
head(physdata)
dim(physdata)
str(physdata)

#make treatment factor and rz_gen (genotype) factor for later stats
physdata <-cbind(physdata[0:20],"treatment"=paste(physdata$pco2, physdata$temp, sep="_"))
physdata <-cbind(physdata[0:21],"gen_rz"=paste(physdata$gen, physdata$rz, sep="_"))

#relevel and rename factors
levels(physdata$time)
physdata$time = factor(physdata$time, levels = c("T0", "T30", "T60", "T90"))

levels(physdata$treatment)
physdata$treatment = factor(physdata$treatment,levels = c("280_28", "280_31", "400_28", "400_31", "700_28", "700_31", "2800_28", "2800_31"))
levels(physdata$treatment) <- c("358uatm, 28°C","325uatm, 31°C","471uatm, 28°C","424uatm, 31°C","674uatm, 28°C","606uatm, 31°C","2750uatm, 28°C","2917uatm, 31°C")

physdata$pco2 = factor(physdata$pco2)
physdata$temp = factor(physdata$temp)
physdata$gen_rz = factor(physdata$gen_rz)
physdata$gen = factor(physdata$gen)
physdata$rz = factor(physdata$rz)

#re-name pco2 levels for stats
levels(physdata$pco2)
levels(physdata$pco2) <- c("ambient","ambient","next century","extreme")

#remove genotypes A and B (not included in this experiment)
physdata <- subset(physdata, gen != "A")
physdata <- subset(physdata, gen != "B")

#### Symbiont density ####

class(physdata$zoox_cm2)
#remove outliers up top in order to plot all physdata together
#treatment filters are the same as the calcification data above
#SSID
sidzoox <- physdata %>%
	select(id, rz, pco2, gen_rz, temp, spp, gen, time, treatment, zoox_cm2) %>%
	filter(spp=="S")%>%
  mutate(zoox_cm2_div = zoox_cm2/1000000) %>% #makes y-axes easier to read later
	filter(id != "SFPZ23")%>%
	filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
  drop_na()
summary(sidzoox)
str(sidzoox)

#PSTR
dipzoox <- physdata %>%
	select(id, rz, gen_rz, pco2, temp, spp, gen, time, treatment, zoox_cm2) %>%
	filter(spp=="P")%>%
  mutate(zoox_cm2_div = zoox_cm2/1000000) %>% #makes y-axes easier to read later
	filter(id != "SNPZ29")%>%
	filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
  mutate(zoox_cm2_cube = zoox_cm2^(1/3)) %>%
  drop_na()
summary(dipzoox)
str(dipzoox)

# check normality and homogeneity of variance
# PSTR
shapiro.test(dipzoox$zoox_cm2) #Dips are not normal
#data:  dipzoox$zoox_cm2
# W = 0.77064, p-value = 4.66e-12
shapiro.test(log(dipzoox$zoox_cm2)) #Dips are not normal
# data:  log(dipzoox$zoox_cm2)
# W = 0.90388, p-value = 5.483e-07
shapiro.test(sqrt(dipzoox$zoox_cm2)) #Dips are not normal
# data:  sqrt(dipzoox$zoox_cm2)
# W = 0.94493, p-value = 0.0001423
dipzoox_cube <- dipzoox$zoox_cm2^(1/3) # Cube root transformatin works
shapiro.test(dipzoox_cube)
#data:  dipzoox_cube
#W = 0.97862, p-value = 0.06494

#qq plot
ggqqplot(dipzoox_cube, main = "qq Plot cube root(PSTR Syms)")

#density plot
ggdensity(dipzoox_cube, main = "Density Plot cube root(PSTR Syms)")

# SSID
shapiro.test(sidzoox$zoox_cm2) #Sids are not normal
# data:  sidzoox$zoox_cm2
# W = 0.89999, p-value = 4.91e-08
shapiro.test(log(sidzoox$zoox_cm2)) #Log transformation works
# data:  log(sidzoox$zoox_cm2)
# W = 0.98361, p-value = 0.1053

#qq plot
ggqqplot(log(sidzoox$zoox_cm2), main = "qq Plot log(SSID Syms)")
#density plot
ggdensity(log(sidzoox$zoox_cm2), main = "Density Plot log(SSID Syms)")

#### Symbiont density model selection ####
#SSID null model
m0 <- lmer(log(zoox_cm2) ~ 1 + (1|gen_rz), data=sidzoox)
summary(m0)
AIC(m0)
#[1] 298.83

#test adding parameters and check AIC
m1 <- lmer(log(zoox_cm2) ~ time + temp + pco2 + time:pco2 + (1|gen_rz), data=sidzoox)
summary(m1)
anova(m1)
AIC(m1)
#[1] 293.82 - time + (1|gen_rz)
#[1] 289.63 - time + temp + (1|gen_rz)
#[1] 292.13 - time + temp + pco2 + (1|gen_rz)
#[1] 288.84 - time+temp+pco2 + time:pco2 + (1|gen_rz)
#doing the interaction of temp*time*pco2 reduced aic but produced lots of insignificant factors in the anova 

#post-hoc tests for significant parameters
lsmeans(m1, pairwise~pco2*time, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "log(SSID Sym Density) Residuals")

#PSTR null model
m0 <- lmer(zoox_cm2_cube ~ 1 + (1|gen_rz), data=dipzoox)
summary(m0)
AIC(m0)
#[1] 1152.774

#test adding parameters and check AIC
m1 <- lmer(log(zoox_cm2) ~ temp + time + (1|gen_rz), data=dipzoox)
summary(m1)
anova(m1)
AIC(m1)
#[1] 360 - temp + (1|gen_rz)
#[1] 344.9 - time + temp + (1|gen_rz)
#adding any more variables reduced the aic but they were not significant

#post-hoc tests for significant parameters
lsmeans(m1, pairwise~time, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "log(PSTR) Sym Density Residuals")

#Summarize data to plot
Sidzoox.means <- summarySE(sidzoox, measurevar="zoox_cm2_div", groupvars=c("time","temp","pco2"))
Dipzoox.means <- summarySE(dipzoox, measurevar="zoox_cm2_div", groupvars=c("time","temp","pco2"))

#Plot
#SSID plot
p.sidzoox <- ggplot(Sidzoox.means, aes(x = time, y = zoox_cm2_div, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + #un-comment if you want to look at all individuals as points
  geom_errorbar(aes(x = time, ymax = zoox_cm2_div+se, ymin = zoox_cm2_div-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  #ggtitle("SSID-Symbiont Density")+
  ylim(0,8)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  xlab("Time")+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.sidzoox
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/SSID_zoox_combinedpco2.pdf", p.sidzoox, width=6, height=5, units=c("in"), useDingbats=FALSE)

#PSTR plot
p.dipzoox <- ggplot(Dipzoox.means, aes(x = time, y = zoox_cm2_div, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + 
  geom_errorbar(aes(x = time, ymax = zoox_cm2_div+se, ymin = zoox_cm2_div-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  #ggtitle("PSTR-Symbiont Density")+
  ylim(0,8)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  xlab("Time")+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.dipzoox
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/PSTR_zoox_combinedpco2.pdf", p.dipzoox, width=6, height=5, units=c("in"), useDingbats=FALSE)



#### Chlorophyll-a Density #### 
class(physdata$chlA_ug_cm2)

#SSID data
sidchl <- physdata %>%
	dplyr::select(id, rz, spp, pco2, temp, gen, gen_rz, time, treatment, chlA_ug_cm2) %>%
	filter(spp=="S")%>%
	filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
  mutate(chlA_cubert = chlA_ug_cm2^(1/3)) %>%
  drop_na()
summary(sidchl)

#PSTR data
dipchl <- physdata %>%
	dplyr::select(id, rz, spp, pco2, temp, gen, gen_rz, time, treatment, chlA_ug_cm2) %>%
	filter(spp=="P")%>%
	filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
  mutate(chlA_sqrt = sqrt(chlA_ug_cm2)) %>%
  drop_na()
summary(dipchl)

# check normality and homogeneity of variance
# PSTR
shapiro.test(dipchl$chlA_ug_cm2) #Dips are not normal
# data:  dipchl$chlA_ug_cm2
# W = 0.95256, p-value = 0.000432
shapiro.test(log(dipchl$chlA_ug_cm2)) #Dips are not normal
# data:  log(dipchl$chlA_ug_cm2)
# W = 0.9447, p-value = 0.0001191
shapiro.test(sqrt(dipchl$chlA_ug_cm2)) #square root transformation works
# data:  sqrt(dipchl$chlA_ug_cm2)
# W = 0.99006, p-value = 0.5652

#qq plot
ggqqplot(sqrt(dipchl$chlA_ug_cm2), main = "qq Plot sqrt(PSTR ChlA)")

#density plot
ggdensity(sqrt(dipchl$chlA_ug_cm2), main = "Density Plot sqrt(PSTR ChlA)")

# SSID
shapiro.test(sidchl$chlA_ug_cm2) #Sids are not normal
# data:  sidchl$chlA_ug_cm2
# W = 0.89695, p-value = 2.573e-08
shapiro.test(log(sidchl$chlA_ug_cm2)) #Sids are not normal
# data:  log(sidchl$chlA_ug_cm2)
# W = 0.84395, p-value = 8.7e-11
shapiro.test(sqrt(sidchl$chlA_ug_cm2)) #Sids are not normal
# data:  sqrt(sidchl$chlA_ug_cm2)
# W = 0.97856, p-value = 0.02856
sidchl_cube <- sidchl$chlA_ug_cm2^(1/3) #cube root transformation works
shapiro.test(sidchl_cube)
# data:  sidchl_cube
# W = 0.98287, p-value = 0.08131

#qq plot
ggqqplot(sidchl_cube, main = "qq Plot cube root(SSID ChlA)")

#density plot
ggdensity(sidchl_cube, main = "Density Plot cube root(SSID ChlA)")

#### Chl-a Model Selection ####
#SSID null model  
m0 <- lmer(chlA_cubert ~ 1 + (1|gen_rz), data=sidchl)
summary(m0)
AIC(m0)
#[1] 178.1507

#test adding parameters and check AIC
m1 <- lmer(log(chlA_ug_cm2) ~ time + pco2 + (1|gen_rz), data=sidchl)
summary(m1)
anova(m1)
AIC(m1)
#[1] 303.7958 - time + (1|gen_rz)
#[1] 298.8724 - time + pco2 + (1|gen_rz)
#adding in time:pco2 really reduces aic value but p = 0.06 for that interaction, so don't include

#post-hoc tests for significant parameters
lsmeans(m1, pairwise~time, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "SSID ChlA Residuals")

#PSTR null model
m0 <- lmer(chlA_sqrt ~ 1 + (1|gen_rz), data=dipchl)
summary(m0)
AIC(m0)
#[1] 343.8846

#test adding parameters and check AIC
m1 <- lmer(chlA_sqrt ~ temp + pco2 + (1|gen_rz), data=dipchl)
summary(m1)
anova(m1)
AIC(m1)
#[1] 315.4296 - temp + (1|gen_rz)
#[1] 315.0349 - temp + pco2 + (1|gen_rz)
#adding in temp:pco2 really reduces aic value but p = 0.78 for that interaction, so not including

#post-hoc tests for significant parameters
lsmeans(m1, pairwise~pco2, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "log(PSTR ChlA) Residuals")


#Summarize data to plot
Sidchl.means <- summarySE(sidchl, measurevar="chlA_ug_cm2", groupvars=c("time","temp","pco2"))
Dipchl.means <- summarySE(dipchl, measurevar="chlA_ug_cm2", groupvars=c("time","temp","pco2"))

#Plot
#SSID
p.sidchl <- ggplot(Sidchl.means, aes(x = time, y = chlA_ug_cm2, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + #un-comment if you want to look at all individuals as points
  geom_errorbar(aes(x = time, ymax = chlA_ug_cm2+se, ymin = chlA_ug_cm2-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  #ggtitle("SSID-Symbiont Density")+
  ylim(0,32)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Chla ("*mu~'g' ~cm^-2~')'))+
  xlab("Time")+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.sidchl
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/SSID_chl_combinedpco2.pdf", p.sidchl, width=6, height=5, units=c("in"), useDingbats=FALSE)

#PSTR
p.dipchl <- ggplot(Dipchl.means, aes(x = time, y = chlA_ug_cm2, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + 
  geom_errorbar(aes(x = time, ymax = chlA_ug_cm2+se, ymin = chlA_ug_cm2-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  #ggtitle("PSTR-Symbiont Density")+
  ylim(0,32)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Chla ("*mu~'g' ~cm^-2~')'))+
  xlab("Time")+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.dipchl
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/PSTR_chl_combinedpco2.pdf", p.dipchl, width=6, height=5, units=c("in"), useDingbats=FALSE)



#### Protein ####
class(physdata$prot_mg_cm2)

#SSID
sidprot <- physdata %>%
	dplyr::select(id, rz, gen_rz, spp, pco2, temp, gen, time, treatment, prot_mg_cm2) %>%
	filter(spp=="S")%>%
	filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
  drop_na()
summary(sidprot)

#PSTR
dipprot <- physdata %>%
	dplyr::select(id, rz, gen_rz, spp, pco2, temp, gen, time, treatment, prot_mg_cm2) %>%
	filter(spp=="P")%>%
	filter(treatment!="325uatm, 31°C")%>%
	filter(treatment!="471uatm, 28°C")%>%
  drop_na()
summary(dipprot)

# check normality and homogeneity of variance
# PSTR
shapiro.test(dipprot$prot_mg_cm2) #Dips are normal
# data:  dipprot$prot_mg_cm2
# W = 0.98625, p-value = 0.2812

#qq plot
ggqqplot(dipprot$prot_mg_cm2, main = "qq Plot PSTR Protein")
#density plot
ggdensity(dipprot$prot_mg_cm2, main = "Density Plot PSTR Protein")

# SSID
shapiro.test(sidprot$prot_mg_cm2) #Sids are normal
# data:  sidprot$prot_mg_cm2
# W = 0.98122, p-value = 0.05448

#qq plot
ggqqplot(sidprot$prot_mg_cm2, main = "qq Plot SSID Protein")
#density plot
ggdensity(sidprot$prot_mg_cm2, main = "Density Plot SSID Protein")

#### Protein Model Selection ####
#SSID null model
m0 <- lmer(prot_mg_cm2 ~ 1 + (1|gen_rz), data=sidprot)
summary(m0)
AIC(m0)
#[1] -2.54317

#test adding parameters and check AIC
m1 <- lmer(prot_mg_cm2 ~ time + temp + (1|gen_rz), data=sidprot)
summary(m1)
anova(m1)
AIC(m1)
#[1] -10.17901 - time + (1|gen_rz)
#[1] -10.27894 - time + temp + (1|gen_rz)
#adding in temp doesn't really change aic value but p = 0.008 for that interaction

#post-hoc tests for significant parameters
lsmeans(m1, pairwise~time, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "SSID Protein Residuals")

#PSTR null model
m0 <- lmer(prot_mg_cm2 ~ 1 + (1|gen_rz), data=dipprot)
summary(m0)
AIC(m0)
#[1] 29.44885

#test adding parameters and check AIC
m1 <- lmer(prot_mg_cm2 ~ temp + (1|gen_rz), data=dipprot)
summary(m1)
anova(m1)
AIC(m1)
#[1] 10.14369 - temp + (1|gen_rz)

#post-hoc tests for significant parameters
lsmeans(m1, pairwise~temp, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "PSTR Protein Residuals")

#SummarySE for figures
Sidprot.means <- summarySE(sidprot, measurevar="prot_mg_cm2", groupvars=c("time","temp","pco2"))
Dipprot.means <- summarySE(dipprot, measurevar="prot_mg_cm2", groupvars=c("time","temp","pco2"))

#Plot
#SSID
p.sidprot <- ggplot(Sidprot.means, aes(x = time, y = prot_mg_cm2, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + #un-comment if you want to look at all individuals as points
  geom_errorbar(aes(x = time, ymax = prot_mg_cm2+se, ymin = prot_mg_cm2-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  ylim(0,1)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  xlab('')+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.sidprot
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/SSID_prot_combinedpco2.pdf", p.sidprot, width=6, height=5, units=c("in"), useDingbats=FALSE)

#PSTR
p.dipprot <- ggplot(Dipprot.means, aes(x = time, y = prot_mg_cm2, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + 
  geom_errorbar(aes(x = time, ymax = prot_mg_cm2+se, ymin = prot_mg_cm2-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  ylim(0,1)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  xlab("")+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.dipprot
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/PSTR_prot_combinedpco2.pdf", p.dipprot, width=6, height=5, units=c("in"), useDingbats=FALSE)


#Now to look at protein separated by reef zone through time:
#summarySE for figures 
Sidprot.means <- summarySE(sidprot, measurevar="prot_mg_cm2", groupvars=c("rz", "time","treatment"))
summary(Sidprot.means)
Dipprot.means <- summarySE(dipprot, measurevar="prot_mg_cm2", groupvars=c("rz","time","treatment"))
summary(Dipprot.means)

#Plot
#SSID
sidprot <- ggplot(Sidprot.means, aes(x = time, y = prot_mg_cm2, fill = rz)) +
  geom_errorbar(aes(x = time, ymax = prot_mg_cm2+se, ymin = prot_mg_cm2-se), width = .2, position = position_dodge(width=0.2)) +
  geom_point(aes(color=rz), position = position_dodge(width = 0.2), size=3)+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  xlab("Time")+
  ylim(-.5,1.5)+
  theme_cowplot()+
  scale_color_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  scale_fill_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  facet_wrap(treatment~., nrow=3, ncol=2)+
  #geom_smooth(aes(group=rz, color=rz, fill=rz), position = position_dodge(width=0.2))+
  stat_smooth(method="lm",se=FALSE, formula = y~poly(x,2),size=1, aes(group=rz, color=rz))+
  geom_hline(yintercept=0, linetype="dashed", color="gray", size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
sidprot
ggsave(file="~/Documents/BU/NEU_experiment/figs/prot/SSID_prot_rz.pdf", sidprot, width=6, height=5, units=c("in"), useDingbats=FALSE)

#PSTR
dipprot <- ggplot(Dipprot.means, aes(x = time, y = prot_mg_cm2, fill = rz)) +
  geom_errorbar(aes(x = time, ymax = prot_mg_cm2+se, ymin = prot_mg_cm2-se), width = .2, position = position_dodge(width=0.2)) +
  geom_point(aes(color=rz), position = position_dodge(width = 0.2), size=3)+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  xlab("Time")+
  ylim(-.5,1.5)+
  theme_cowplot()+
  scale_color_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  scale_fill_manual(values=c("F"="slategray4", "N"="black"), name = "Reef Zone", labels = c("FR", "NS"))+
  facet_wrap(treatment~., nrow=3, ncol=2)+
  #geom_smooth(aes(group=rz, color=rz, fill=rz), position = position_dodge(width=0.2))+
  stat_smooth(method="lm",se=FALSE, formula = y~poly(x,2),size=1, aes(group=rz, color=rz))+
  geom_hline(yintercept=0, linetype="dashed", color="gray", size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
dipprot
ggsave(file="~/Documents/BU/NEU_experiment/figs/prot/PSTR_prot_rz.pdf", dipprot, width=6, height=5, units=c("in"), useDingbats=FALSE)



#### Carbohydrate ####
class(physdata$carbs_mg_mL)

#SSID
sidcarb <- physdata %>%
  dplyr::select(id, rz, gen_rz, spp, pco2, temp, gen, time, treatment, carbs_mg_mL) %>%
  filter(spp=="S")%>%
  filter(treatment!="325uatm, 31°C")%>%
  filter(treatment!="471uatm, 28°C")%>%
  mutate(carbs_sqrt = sqrt(carbs_mg_mL)) %>%
  drop_na()
summary(sidcarb)

#PSTR
dipcarb <- physdata %>%
  dplyr::select(id, rz, gen_rz, spp, pco2, temp, gen, time, treatment, carbs_mg_mL) %>%
  filter(spp=="P")%>%
  filter(treatment!="325uatm, 31°C")%>%
  filter(treatment!="471uatm, 28°C")%>%
  mutate(carbs_sqrt = sqrt(carbs_mg_mL)) %>%
  drop_na()
summary(dipcarb)

# check for normality
# PSTR
shapiro.test(dipcarb$carbs_mg_mL) #Dips are not normal
# data:  dipcarb$carbs_mg_mL
# W = 0.97327, p-value = 0.02198
shapiro.test(log(dipcarb$carbs_mg_mL)) #Dips are not normal
# data:  log(dipcarb$carbs_mg_mL)
# W = 0.88099, p-value = 4.382e-08
shapiro.test(sqrt(dipcarb$carbs_mg_mL)) #square root transformation works
# data:  sqrt(dipcarb$carbs_mg_mL)
# W = 0.99072, p-value = 0.6363

#qq plot
ggqqplot(sqrt(dipcarb$carbs_mg_mL), main = "qq Plot sqrt(PSTR Carbs)")

#density plot
ggdensity(sqrt(dipcarb$carbs_mg_mL), main = "Density Plot sqrt(PSTR Carbs)")

shapiro.test(sidcarb$carbs_mg_mL) #Sids are not normal
# data:  sidcarb$carbs_mg_mL
# W = 0.97654, p-value = 0.01828
shapiro.test(log(sidcarb$carbs_mg_mL)) #Sids are not normal
# data:  log(sidcarb$carbs_mg_mL)
# W = 0.93582, p-value = 6.46e-06
shapiro.test(sqrt(sidcarb$carbs_mg_mL)) #square root transformation works
# data:  sqrt(sidcarb$carbs_mg_mL)
# W = 0.98721, p-value = 0.2353

#qq plot
ggqqplot(sqrt(sidcarb$carbs_mg_mL), main = "qq Plot sqrt(SSID Carbs)")

#density plot
ggdensity(sqrt(sidcarb$carbs_mg_mL), main = "Density Plot sqrt(SSID Carbs)")

#### Carb Model Selection ####
#SSID null model 
m0 <- lmer(carbs_sqrt ~ 1 + (1|gen_rz), data=sidcarb)
summary(m0)
AIC(m0)
#[1] 171.4

#test adding parameters and check AIC
m1 <- lmer(carbs_sqrt ~ temp + (1|gen_rz), data=sidcarb)
summary(m1)
anova(m1)
AIC(m1)
#[1] 168.2067 - temp + (1|gen_rz)

#post hoc tests
lsmeans(m1, pairwise~temp, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "sqrt(SSID Carbohydrate) Residuals")

#PSTR null model
m0 <- lmer(carbs_sqrt ~ 1 + (1|gen_rz), data=dipcarb)
summary(m0)
AIC(m0)
#[1] 160.1132

#test adding parameters and check AIC
m1 <- lmer(carbs_sqrt ~ temp + time + time:temp + (1|gen_rz), data=dipcarb)
summary(m1)
anova(m1)
AIC(m1)
#[1] 144.5942 - temp + (1|gen_rz)
#[1] 142.6836 - temp + time + (1|gen_rz)
#[1] 142.4132 - temp + time + time:temp + (1|gen_rz)

#post hoc tests
lsmeans(m1, pairwise~time, adjust="tukey")

#look at homogeneity of variance
par(mfrow=c(2,2))
plot(m1, main = "log(PSTR Carbohydrate) Residuals")

#summarySE for plotting
Sidcarb.means <- summarySE(sidcarb, measurevar="carbs_mg_mL", groupvars=c("time","temp","pco2"))
Dipcarb.means <- summarySE(dipcarb, measurevar="carbs_mg_mL", groupvars=c("time","temp","pco2"))

#Plot
#SSID
p.sidcarb <- ggplot(Sidcarb.means, aes(x = time, y = carbs_mg_mL, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + #un-comment if you want to look at all individuals as points
  geom_errorbar(aes(x = time, ymax = carbs_mg_mL+se, ymin = carbs_mg_mL-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  ylim(0,8)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Total Carbohydrate (mg" ~mL^-1~')'))+
  xlab('')+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.sidcarb
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/SSID_carb_combinedpco2.pdf", p.sidcarb, width=6, height=5, units=c("in"), useDingbats=FALSE)

p.dipcarb <- ggplot(Dipcarb.means, aes(x = time, y = carbs_mg_mL, shape=pco2)) +
  geom_point(aes(color=pco2),size = 3.5, position = position_dodge(width=0.2)) + 
  geom_errorbar(aes(x = time, ymax = carbs_mg_mL+se, ymin = carbs_mg_mL-se, color=pco2), width = .2, position = position_dodge(width=0.2)) +
  ylim(0,8)+
  scale_shape_manual(values=c("ambient"=16, "ambient"=16,"next century"=17, "extreme"=15), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  scale_fill_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+ #green = ambient pco2, orange = 700, purple = 2800
  scale_color_manual(values=c("ambient"="#1B9E77", "ambient"="#1B9E77","next century"="#D95F02", "extreme"="#7570B3"), name = "pCO2 (uatm)", labels = c("ambient", "next century","extreme"))+
  ylab(bquote("Total Carbohydrate (mg" ~mL^-1~')'))+
  xlab("")+
  facet_wrap(temp~., nrow=1, ncol=2)+
  geom_smooth(aes(group=pco2, color=pco2),method="lm",se=FALSE, size=1)+
  theme(legend.title = element_text(size = 12), legend.text=element_text(size = 10), legend.position = "right", legend.background = element_rect(color = "black", linetype = "solid", size = .5))
p.dipcarb
ggsave(file="~/Documents/BU/NEU_experiment/figs/final/PSTR_carb_combinedpco2.pdf", p.dipcarb, width=6, height=5, units=c("in"), useDingbats=FALSE)


#### Correlation Matrices ####

# Join All Data
head(physdata)
head(calcdata)

all_data <- physdata %>%
  left_join(calcdata, by = c("id","time","spp","treatment","rz","pco2","temp", "tank", "gen", "sur_area","gen_rz"))%>%
  #filter(!is.na(Growthmg.cm2.day)) %>% #took this out to keep T0 measurements for other physiologies
  filter(treatment!="325uatm, 31°C") %>%
  filter(treatment!="471uatm, 28°C") %>%
  filter(id !="SFPZ23") %>%
  filter(id !="SNPZ29") %>%
  filter(id !="SNSY27") %>%
  filter(id !="SNPZ39") %>%
  filter(gen != "A")%>%
  filter(gen != "B")

summary(all_data)
head(all_data)
dim(all_data)


#write.csv(all_data, "AllPhysData_combined.csv")

#separate out species and time points
time30sid <- all_data %>%
  filter(time=="T30") %>%
  filter(spp=="S") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time30dip <- all_data %>%
  filter(time=="T30") %>%
  filter(spp=="P") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time60sid <- all_data %>%
  filter(time=="T60") %>%
  filter(spp=="S") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2)%>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time60dip <- all_data %>%
  filter(time=="T60") %>%
  filter(spp=="P") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time90sid <- all_data %>%
  filter(time=="T90") %>%
  filter(spp=="S") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time90dip <- all_data %>%
  filter(time=="T90") %>%
  filter(spp=="P") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()


#source the rquery.cormat function, will need it below
source("http://www.sthda.com/upload/rquery_cormat.r")

#Function to calculate p-value of correlations for corrplot()
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#PSTR correlations
#PSTR time 30
#compute correlation matrix
time30dip_m <- cor(time30dip) #create a correlation matrix of data
time30dip_p.mat <- cor.mtest(time30dip) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
time30dip_table <- rquery.cormat(time30dip, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/NEU_experiment/Correlation_Matrices/PSTR_t30_corrmat.pdf", width=6, height=6, useDingbats=FALSE)
corrplot(time30dip_m, method = "circle", 
         type = "upper", 
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = time30dip_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
) 
dev.off()

#PSTR time 60
time60dip_m <- cor(time60dip) #create a correlation matrix of data
time60dip_p.mat <- cor.mtest(time60dip) #calculate significance of correlations

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
time60dip_table <- rquery.cormat(time60dip, type="flatten", graph=FALSE)

pdf(file = "/Users/hannahaichelman/Documents/BU/NEU_experiment/Correlation_Matrices/PSTR_t60_corrmat.pdf", width=6, height=6, useDingbats=FALSE)
corrplot(time60dip_m, method = "circle", 
         type = "upper", 
         order="alphabet",
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = time60dip_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
) 
dev.off()

#PSTR time 90
time90dip_m <- cor(time90dip) #create a correlation matrix of data
time90dip_p.mat <- cor.mtest(time90dip) #calculate significance of correlations

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
time90dip_table <- rquery.cormat(time90dip, type="flatten", graph=FALSE)

pdf(file = "/Users/hannahaichelman/Documents/BU/NEU_experiment/Correlation_Matrices/PSTR_t90_corrmat.pdf", width=6, height=6, useDingbats=FALSE)
corrplot(time90dip_m, method = "circle", 
         type = "upper", 
         order="alphabet",
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = time90dip_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         ) 
dev.off()

#SSID correlations
#SSID time 30
#compute correlation matrix
time30sid_m <- cor(time30sid)
time30sid_p.mat <- cor.mtest(time30sid) #calculate significance of correlations

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
time30sid_table <- rquery.cormat(time30sid, type="flatten", graph=FALSE)

#plot correlation matrix
pdf(file = "/Users/hannahaichelman/Documents/BU/NEU_experiment/Correlation_Matrices/SSID_t30_corrmat.pdf", width=6, height=6, useDingbats=FALSE)
corrplot(time30sid_m, method = "circle", 
         type = "upper", 
         order="alphabet",
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = time30sid_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
) 
dev.off()

#SSID time 60
time60sid_m <- cor(time60sid)
time60sid_p.mat <- cor.mtest(time60sid) #calculate significance of correlations

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
time60sid_table <- rquery.cormat(time60sid, type="flatten", graph=FALSE)

pdf(file = "/Users/hannahaichelman/Documents/BU/NEU_experiment/Correlation_Matrices/SSID_t60_corrmat.pdf", width=6, height=6, useDingbats=FALSE)
corrplot(time60sid_m, method = "circle", 
         type = "upper", 
         order="alphabet",
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = time60sid_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
) 
dev.off()

#SSID time 90
time90sid_m <- cor(time90sid)
time90sid_p.mat <- cor.mtest(time90sid) #calculate significance of correlations

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
time90sid_table <- rquery.cormat(time90sid, type="flatten", graph=FALSE)

pdf(file = "/Users/hannahaichelman/Documents/BU/NEU_experiment/Correlation_Matrices/SSID_t90_corrmat.pdf", width=6, height=6, useDingbats=FALSE)
corrplot(time90sid_m, method = "circle", 
         type = "upper", 
         order="alphabet",
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = time90sid_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
) 
dev.off()

#### Linear regression/Correlation analysis ####

#separate out species and time points
#slightly different for the linear regressions compared to above, we need to keep more info
time30sid <- all_data %>%
  filter(time=="T30") %>%
  filter(spp=="S") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2, temp, rz, pco2, gen_rz) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time30dip <- all_data %>%
  filter(time=="T30") %>%
  filter(spp=="P") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2, temp, rz, pco2, gen_rz) %>%
  mutate(zoox_cm2_div = zoox_cm2/1000000) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2_div) %>%
  drop_na()

time60sid <- all_data %>%
  filter(time=="T60") %>%
  filter(spp=="S") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2, temp, pco2, gen_rz)%>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time60dip <- all_data %>%
  filter(time=="T60") %>%
  filter(spp=="P") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2, temp, pco2, rz, gen_rz) %>%
  mutate(zoox_cm2_div = zoox_cm2/1000000) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2_div) %>%
  drop_na()

time90sid <- all_data %>%
  filter(time=="T90") %>%
  filter(spp=="S") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2, temp, pco2, gen_rz) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2) %>%
  drop_na()

time90dip <- all_data %>%
  filter(time=="T90") %>%
  filter(spp=="P") %>%
  select(Growthmg.cm2.day,carbs_mg_cm2,prot_mg_cm2,chlA_ug_cm2,zoox_cm2, temp, pco2, rz, gen_rz) %>%
  mutate(zoox_cm2_div = zoox_cm2/1000000) %>%
  rename(calcification = Growthmg.cm2.day, carbs = carbs_mg_cm2, protein = prot_mg_cm2, Chla = chlA_ug_cm2, syms = zoox_cm2_div) %>%
  drop_na()

time90dip_28 <- time90dip %>%
  filter(temp=="28")
time90dip_31 <- time90dip %>%
  filter(temp=="31")


##Calcification ~ protein regression

#PSTR plot, color by temperature treatment and shape is reef zone
#just change the time of the object plotting as needed
plot <- ggplot(time30dip, aes(x=protein,y=calcification, color = temp, shape=rz)) +
  geom_point(aes(color = temp, shape=rz), size = 2.5)+
  ggtitle("PSTR T30")+
  geom_smooth(method = "lm", col = "black", aes(group=temp))+ #group=1 creates only one regression line
  scale_color_manual(values=c("28"="blue", "31"="red"))+
  scale_shape_manual(values=c("F"=15, "N"=17))+
  #facet_wrap(temp~.)+
  theme_cowplot()
plot
ggsave(file="~/Documents/BU/NEU_experiment/Regression_Analysis/by_treatment/PSTR/PSTR_calc_prot_bytemp+rz_T30_nofacet.pdf", plot, width=5, height=5, units=c("in"), useDingbats=FALSE)

  
#lmer to get significance of correlation
str(time90dip)
time90dip$gen_rz = factor(time90dip$gen_rz) #make genotype a factor

#full model
fits1 <- lmer(calcification ~ protein*temp + (1|gen_rz), data=time90dip)
#summarize model
summary(fits1)
AIC(fits1) #check AIC value of the model
Anova(fits1) #Type II Wald chisquare test
r.squaredGLMM(fits1) #MuMIn package to get the r-squared value
#R2m = marginal r squared, represents variance explained by fixed effects
#R2c = conditional r squared, variance explained by the entire model (both fixed and random effects)


fits2 <- lmer(calcification ~ protein+temp + (1|gen_rz), data=time60dip)
summary(fits2)

fits3 <- lmer(calcification ~ protein + (1|gen_rz), data=time60dip)
summary(fits3)

fits4 <- lmer(calcification ~ temp + (1|gen_rz), data=time60dip)
summary(fits4)

#Now we compare the increasingly parsimonious, nested models listed above with likelihood ratio tests 
#to estimate p values of the predictors and their interactions (with anova() below)
#this chi-squre pvalue is indicitive of interaction between protein and temperature, 
#because this is the only thing we changed between fits and fits2
anova(fits1,fits2)

#Get r-squared for each temperature
fits_28 <- lmer(calcification ~ protein + (1|gen_rz), data=time90dip_28)
r.squaredGLMM(fits_28)
fits_31 <- lmer(calcification ~ protein + (1|gen_rz), data=time90dip_31)
r.squaredGLMM(fits_31)

#post-hoc tests
lsmeans(fits1, pairwise~temp, adjust="tukey")


## Calcification ~ carbohydrate regression
plot <- ggplot(time30dip, aes(x=carbs,y=calcification, color = temp, shape=rz)) +
  geom_point(aes(color = temp, shape=rz), size = 2.5)+
  ggtitle("PSTR T30")+
  scale_color_manual(values=c("28"="blue", "31"="red"))+
  scale_shape_manual(values=c("F"=15, "N"=17))+
  #facet_wrap(temp~.)+
  stat_smooth(method = "lm", col = "black", aes(group=temp))+
  theme_cowplot()
plot
ggsave(file="~/Documents/BU/NEU_experiment/Regression_Analysis/by_treatment/PSTR/PSTR_calc_carb_T30_bytemp+rz_nofacet.pdf", plot, width=5, height=5, units=c("in"), useDingbats=FALSE)


#full model
fits1 <- lmer(calcification ~ carbs*temp + (1|gen_rz), data=time60dip)
summary(fits1)
AIC(fits1)
Anova(fits1) #Type II Wald chisquare test
r.squaredGLMM(fits1)

fits2 <- lmer(calcification ~ carbs+temp + (1|gen_rz), data=time60dip)
summary(fits2)

fits3 <- lmer(calcification ~ carbs + (1|gen_rz), data=time60dip)
summary(fits3)

fits4 <- lmer(calcification ~ temp + (1|gen_rz), data=time60dip)
summary(fits4)

#this pvalue is indicitive of interaction between protein and temperature, 
#because this is the only thing we changed between fits and fits2
anova(fits1,fits2)

#Get r-squared for each temperature
fits_28 <- lmer(calcification ~ carbs + (1|gen_rz), data=time90dip_28)
r.squaredGLMM(fits_28)
fits_31 <- lmer(calcification ~ carbs + (1|gen_rz), data=time90dip_31)
r.squaredGLMM(fits_31)


## Calcification ~ symbiont density regression
plot <- ggplot(time30dip, aes(x=syms,y=calcification, color = temp, shape=rz)) +
  geom_point(aes(color = temp, shape=rz), size = 2.5)+
  ggtitle("PSTR T30")+
  scale_color_manual(values=c("28"="blue", "31"="red"))+
  scale_shape_manual(values=c("F"=15, "N"=17))+
  #facet_wrap(temp~.)+
  stat_smooth(method = "lm", col = "black", aes(group=temp))+
  theme_cowplot()
plot
ggsave(file="~/Documents/BU/NEU_experiment/Regression_Analysis/by_treatment/PSTR/PSTR_calc_syms_bytemp+rz_T30_nofacet.pdf", plot, width=5, height=5, units=c("in"), useDingbats=FALSE)


#full model
fits1 <- lmer(calcification ~ syms*temp + (1|gen_rz), data=time30dip)
summary(fits1)
AIC(fits1)
Anova(fits1) #Type II Wald chisquare test
r.squaredGLMM(fits1)

fits2 <- lmer(calcification ~ syms+temp + (1|gen_rz), data=time30dip)
summary(fits2)

fits3 <- lmer(calcification ~ syms + (1|gen_rz), data=time30dip)
summary(fits3)

fits4 <- lmer(calcification ~ temp + (1|gen_rz), data=time30dip)
summary(fits4)

anova(fits1,fits2)

#Get r-squared for each temperature
fits_28 <- lmer(calcification ~ syms + (1|gen_rz), data=time90dip_28)
r.squaredGLMM(fits_28)
fits_31 <- lmer(calcification ~ syms + (1|gen_rz), data=time90dip_31)
r.squaredGLMM(fits_31)

## Calcification ~ Chlorophyll a regression

plot <- ggplot(time90dip, aes(x=Chla,y=calcification, color = temp, shape=rz)) +
  geom_point(aes(color = temp, shape=rz), size = 2.5)+
  ggtitle("PSTR T90")+
  scale_color_manual(values=c("28"="blue", "31"="red"))+
  scale_shape_manual(values=c("F"=15, "N"=17))+
  #facet_wrap(temp~.)+
  stat_smooth(method = "lm", col = "black", aes(group=temp))+
  theme_cowplot()
plot
ggsave(file="~/Documents/BU/NEU_experiment/Regression_Analysis/by_treatment/PSTR/PSTR_calc_chla_bytemp+rz_T90_nofacet.pdf", plot, width=5, height=5, units=c("in"), useDingbats=FALSE)

#full model
fits1 <- lmer(calcification ~ Chla*temp + (1|gen_rz), data=time90dip)
summary(fits1)
AIC(fits1)
Anova(fits1) #Type II Wald chisquare test
r.squaredGLMM(fits1)

fits2 <- lmer(calcification ~ Chla+temp + (1|gen_rz), data=time90dip)
summary(fits2)

fits3 <- lmer(calcification ~ Chla + (1|gen_rz), data=time90dip)
summary(fits3)

fits4 <- lmer(calcification ~ temp + (1|gen_rz), data=time90dip)
summary(fits4)

anova(fits1,fits2)

#Get r-squared for each temperature
fits_28 <- lmer(calcification ~ Chla + (1|gen_rz), data=time90dip_28)
r.squaredGLMM(fits_28)
fits_31 <- lmer(calcification ~ Chla + (1|gen_rz), data=time90dip_31)
r.squaredGLMM(fits_31)

