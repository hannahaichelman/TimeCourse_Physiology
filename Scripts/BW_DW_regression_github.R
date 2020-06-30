#Written by Hannah Aichelman
#Contact: hannahaichelman@gmail.com
#Last updated: 6/30/2020

# This script is used to check the correlation between buoyant weight (BW) and dry weight (DW) 
# for coral fragments associated with Aichelman et al. (2020)


#load necessary libraries
library(tidyverse)

#set working directory
setwd("~/Documents/BU/NEU_experiment/BW_DW_Regression")

#read in data
regression <- read.csv("BW_DW_calculation_alldata.csv")

head(regression)

#filter data by species
#including both Sarah's and Colleen's fragments in the regression

#SSID
SSIDreg <- regression %>%
	filter(species == "S")
head(SSIDreg)

#PSTR
PSTRreg <- regression %>%
	filter(species == "P")
head(PSTRreg)

#make SSID scatterplot and linear model
#note the difference between BW/DW and BW2/DW2 is just mg and g, respectively
#I use g for the regression, because that is what the original BW values were taken in
#Then, the measurement is converted to mg before correcting for surface area and days in experiment

scatter.smooth(x= SSIDreg$BW, y= SSIDreg$DW, main="SSID Regression")
SIDlinearmod <- lm(DW2~BW2, data=SSIDreg) #using g
SIDlinearmod
summary(SIDlinearmod)

# Call:
#   lm(formula = DW2 ~ BW2, data = SSIDreg)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.3988 -0.8358 -0.1513  0.4535  3.3896 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.60184    0.58944   6.111 5.69e-08 ***
#   BW2          1.94902    0.07785  25.037  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.128 on 67 degrees of freedom
# Multiple R-squared:  0.9034,	Adjusted R-squared:  0.902 
# F-statistic: 626.9 on 1 and 67 DF,  p-value: < 2.2e-16

# To interpret the output of lm() into an equation: 
# the intercept is the Estimate value in the (Intercept) row
# the slope is the Estimate value in the x, or BW2, row.
# http://www.sthda.com/english/articles/40-regression-analysis/167-simple-linear-regression-in-r/


#make PSTR regression and linear model
scatter.smooth(x= PSTRreg$DW, y= PSTRreg$BW, main="PSTR Regression")
DIPlinearmod <- lm(DW2~BW2, data=PSTRreg)
DIPlinearmod
summary(DIPlinearmod)

# Call:
#   lm(formula = DW2 ~ BW2, data = PSTRreg)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.4425 -1.1370 -0.5291  0.9183  6.9040 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  6.95614    0.56723   12.26   <2e-16 ***
#   BW2          1.63395    0.07939   20.58   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.748 on 95 degrees of freedom
# Multiple R-squared:  0.8168,	Adjusted R-squared:  0.8149 
# F-statistic: 423.6 on 1 and 95 DF,  p-value: < 2.2e-16
