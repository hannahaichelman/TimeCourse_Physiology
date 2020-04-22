#Written by Hannah Aichelman, based on script from Hannah G. Reich
#Contact: hannahaichelman@gmail.com
#Last updated: 4/22/2020

#This script is used to explore and plot PCAs that consider how holobiont physiology changes through time
#The dataset includes host physiology (calcification rate, total protein, and total carbohydrate) 
#as well as symbiont physiology (symbiont density and chlorophyll a content) from two species of 
#Caribbean reef-building corals (Siderastrea siderea and Pseudodiploria strigosa) across reef zone (forereef and nearshore).
#These corals were exposed to temperature (28 and 31C) and acidification (present day, next century, and extreme pco2 conditions)
#stressors for 95 days, and physiology measurements were taken every 30 days. 
#Throughout, SSID = Siderastrea siderea and PSTR = Pseudodiploria strigosa

#Some good resources used in creating this script, along with advice from Hannah Reich:
# http://www.sthda.com/english/wiki/print.php?id=202
# how to center title on grid: https://stackoverflow.com/questions/40675778/center-plot-title-in-ggplot2
# more facto mine r: http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
# for interpreting biplots: https://blog.bioturing.com/2018/06/18/how-to-read-pca-biplots-and-scree-plots/ https://stats.stackexchange.com/questions/61215/how-to-interpret-this-pca-biplot-coming-from-a-survey-of-what-areas-people-are-i
# explanaition on confidence ellipses: https://stats.stackexchange.com/questions/217374/real-meaning-of-confidence-ellipse

#Set working directory
setwd("~/Documents/BU/NEU_experiment/data_sheets")

#Load necessary libraries
library(readxl)
library(ggpubr)
library(ggfortify)
library(ggplot2)
library(cluster)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(dplyr)
library(ggbiplot)
library(PerformanceAnalytics) # this is for the scatter plot matrix
library(corrplot) # this is for the correlation matrix
library(cowplot)
library(vegan)

#Exported this csv file from the TimeCoursePhysiologyAnalysis_github.R script, just compiles all host and symbiont physiology
all_data <- read.csv("AllPhysData_combined.csv")

#### Reorganize and group the data ####

# Note that this is the raw data that is not log transformed
pca_nolog <- all_data[,c(1:3,21,22,10:12,8,15,16,19,23)]
head(pca_nolog)

# Rename the physiology metrics
colnames(pca_nolog)[colnames(pca_nolog)=="zoox_cm2"] <-"syms"
colnames(pca_nolog)[colnames(pca_nolog)=="prot_mg_cm2"] <-"protein"
colnames(pca_nolog)[colnames(pca_nolog)=="chlA_ug_cm2"] <-"Chla"
colnames(pca_nolog)[colnames(pca_nolog)=="carbs_mg_mL"] <-"carbs"
colnames(pca_nolog)[colnames(pca_nolog)=="Growthmg.cm2.day"] <-"calcification"

# Remove missing data and the two treatments we don't consider here
pca_na <- (na.omit(pca_nolog)) %>%
  filter(treatment!="325uatm, 31°C")%>%
  filter(treatment!="471uatm, 28°C")
head(pca_na)  

# Rename pco2 also so we only have 3 levels
pca_na$pco2 = factor(pca_na$pco2)
levels(pca_na$pco2) <- c("ambient","ambient","next century","extreme")

# Make temperature a factor
pca_na$temp = factor(pca_na$temp)

# Make treatment a factor and re-level based on treatment names
levels(pca_na$treatment)
pca_na$treatment = factor(pca_na$treatment,levels = c("358uatm, 28°C","424uatm, 31°C","674uatm, 28°C","606uatm, 31°C","2750uatm, 28°C","2917uatm, 31°C"))

# Add 2 to the calcification column so there aren't zeros when we log transform  
pca_na[,13] <- pca_na[,13] + 2

# log transform the data (selecting only the physiology data)
pca <- log(pca_na[,9:13])

# Merge descriptive stuff and log-transformed stuff back together
all <- cbind(pca_na[,1:8],pca)
head(all)

# Split dataframe by time point only - use these for combined species PCAs
pca_na_t30 <- all %>%
  filter(time == "T30")

pca_na_t60 <- all %>%
  filter(time == "T60")

pca_na_t90 <- all %>%
  filter(time == "T90")

# Split dataframe by species only - use these for separate species PCAs
pca_log_ssid <- all %>%
  filter(spp == "S")

pca_log_pstr <- all %>%
  filter(spp == "P")

# Split dataframe by time points - PSTR (use these for separate species PCAs)
pca_na_pstr_t30 <- pca_log_pstr %>%
  filter(time == "T30")

pca_na_pstr_t60 <- pca_log_pstr %>%
  filter(time == "T60")

pca_na_pstr_t90 <- pca_log_pstr %>%
  filter(time == "T90")

# Split dataframe by time points - SSID (use these for separate species PCAs)
pca_na_ssid_t30 <- pca_log_ssid %>%
  filter(time == "T30")

pca_na_ssid_t60 <- pca_log_ssid %>%
  filter(time == "T60")

pca_na_ssid_t90 <- pca_log_ssid %>%
  filter(time == "T90")


##### Adonis significance tests ####
## Use an Adonis test to get significance of factors on holobiont physiology
library(vegan)
library(MCMC.OTU)

# Dont want to use log transformed data here
head(pca_na)
str(pca_na)

#for combined species analysis, change time filter and then re-run to adonis code below
pca_na_spp <- pca_na %>%
  filter(time=="T90") #change time here for different stats

#for SSID only analysis
pca_na_ssid <- pca_na %>%
  filter(spp=="S") %>%
  filter(time=="T90") #change time here for different stats

#for PSTR only analysis
pca_na_pstr <- pca_na %>%
  filter(spp=="P") %>%
  filter(time=="T90") #change time here for different stats

# Change dataframe here based on the comparison you are interested in
nl=startedLog(data=pca_na_ssid,count.columns=9:13, logstart=1)

goods.dist=vegdist(nl, method="bray")
goods.pcoa=pcoa(goods.dist)

# PCA:
pcp=prcomp(nl, retx=TRUE, center=TRUE)
scores=goods.pcoa$vectors
summary(goods.pcoa)
conditions=pca_na_ssid[,1:8] #make sure to change dataframe here

# PERMANOVA  
head(scores)
head(conditions)

adonis(scores~pco2*temp*rz + gen_rz, data=conditions,method="euclidean", permutations = 10000)  


##### Build PCAs ####
# Make PCAs using facto package

# Make necessary variables factors
# Change the dataframe based on what time point and species being considered
species <- factor(pca_na_pstr_t30[,3]) 
reefzone <- factor(pca_na_pstr_t30[,2])
temp <- factor(pca_na_pstr_t30[,7])
pco2 <- factor(pca_na_pstr_t30[,6])
gen <- factor(pca_na_pstr_t30[,5])
treat <- factor(pca_na_pstr_t30[,4])
time <- factor(pca_na_pstr_t30[,8])

# Set different colors based on the factors you're interested in looking at
# For example, we plot PSTR colors based on temperature, so that would be the 'cols' uncommented here
cols <- c("31" = "red", "28" = "blue") #PSTR temperature comparison
#cols <- c("S"="#FF9500", "P"="#0C5DA5") #species comparison
#cols <- c("F"="#7846B4", "N"="#82B446") #reef zone comparison
#cols <- c("ambient"="#1B9E77", "next century"="#D95F02","extreme"="#7570B3") #SSID pCO2 comparison

# Create object that we will plot PCA with
# Doing PSTR T30 here, but just replace pca_na_pstr_t30 with whatever species/time/etc. comparison we are interested in.
facto3 <- PCA(pca_na_pstr_t30[,9:13], scale.unit = TRUE, ncp = 10, graph = TRUE)

# Facto package gives some cool data, explore that here
# correlations variables - dimensions
facto3$var$cor

# mean of the variables
facto3$call$centre

# standard error of variables
facto3$call$ecart.type

# get eigenvalues
eig.val1 <- get_eigenvalue(facto3)
eig.val1

# visualize percent of variance explained by each dimension
fviz_eig(facto3, addlabels = TRUE)

# bar graphs of contributions for all of the dimensions... for variables
c <- fviz_contrib(facto3, choice = "var", axes = 1, top = 10)
save_plot("pc1_contribution_of_var.png",c, base_aspect_ratio = 1.6)
d <- fviz_contrib(facto3, choice = "var", axes = 2, top = 10)
save_plot("pc2_contribution_of_var.png",d, base_aspect_ratio = 1.6)
# Here, the red line = expected average contribution... "In this case, the expected average contribution (cutoff) 
#is calculated as follows: As mentioned above, if the contributions of the 10 variables were uniform, 
#the expected average contribution on a given PC would be 1/10 = 10%. 
#The expected average contribution of a variable for PC1 and PC2 is : [(10* Eig1) + (10 * Eig2)]/(Eig1 + Eig2)"

# visualize graph of variables in a different way
fviz_pca_var(facto3, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
             )

### describe dimensions
dim_descript <- dimdesc(facto3, axes = c(1,2,3,4,5), proba = 0.05)

# dim1
dim_descript$Dim.1
# dim2
dim_descript$Dim.2
# dim3
dim_descript$Dim.3

#### Plot PCAs ####
# these figures fuse the facto package and ggplot features

# PCA colored by temperature 
# Used this for PSTR through time
#PSTR T30 is plotted here
temp_pca <- fviz_pca_biplot(facto3, 
                     label = "var",
                     col.var = "black", labelsize = 4,
                     alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(shape=pco2, colour=temp), size = 3, stroke = 1) +
  scale_color_manual(values = cols, breaks=c("31", "28"), name = "Temperature") +
  scale_shape_manual(values = c(24,22,21), breaks = c("ambient", "next century", "extreme"), labels = c("ambient", "next century", "extreme"), name = "pCO2")  +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.2, 
               aes(fill= temp), show.legend = FALSE) + scale_fill_manual(values=cols) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  labs(title = "Principle Component Analysis (PCA) of \n P. strigosa Physiology at T30",
       x = "PC1 (76.5% explained variance)", #get this from output of fviz_eig() above
       y = "PC2 (11.8% explained variance)") + #get this from output of fviz_eig() above
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold")) 
temp_pca


#Note that to make any of the PCAs below, would need to re-make facto3 object with the species/time/etc. info you are interested in
#Also, will need to change the colors, which is 'cols' 

# PCA by pCO2 
# Used this for SSID through time
pco2_pca <- fviz_pca_biplot(facto3, 
                     label = "var",
                     col.var = "black", labelsize = 4,
                     alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(shape=temp, colour=pco2), size = 3, stroke = 1) +
  scale_color_manual(values = cols, breaks=c("ambient", "next century", "extreme"), name = "pCO2") +
  scale_shape_manual(values = c(24,22), breaks=c("31", "28"), labels = c("31", "28"), name = "Temperature")  +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.2, 
               aes(fill= pco2), show.legend = FALSE) + scale_fill_manual(values=cols) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  labs(title = "Principle Component Analysis (PCA) of \n S. siderea Physiology at T90",
       x = "PC1 (65% explained variance)",
       y = "PC2 (19.9% explained variance)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold")) 
pco2_pca

# PCA by species 
spp_pca <- fviz_pca_biplot(facto3, 
                     label = "var",
                     col.var = "black", labelsize = 4,
                     alpha.ind = "0") + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(shape=reefzone, colour=species), size = 3, stroke = 1) +
  scale_color_manual(values = cols, breaks=c("P", "S"), name = "Species") +
  scale_shape_manual(values = c(24,22), breaks = c("F", "N"), labels = c("ForeReef", "NearShore"), name = "Reef Zone")  +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.2, 
               aes(fill= species), show.legend = FALSE) + scale_fill_manual(values=cols) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  labs(title = "Principle Component Analysis (PCA) \n Physiology at T90",
       x = "PC1 (71.3% explained variance)",
       y = "PC2 (15% explained variance)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold")) 
spp_pca

# Also explored PCA coloring by reef zone
rz_pca <- fviz_pca_biplot(facto3, 
                     label = "var",
                     col.var = "black", labelsize = 4,
                     alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(shape=temp, colour=reefzone), size = 3, stroke = 1) +
  scale_color_manual(values = cols, breaks=c("F", "N"), name = "Reef Zone") +
  scale_shape_manual(values = c(24,22), breaks = c("31", "28"), labels = c("31°C", "28°C"), name = "Temperature")  +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.125, 
               aes(fill= reefzone), show.legend = FALSE) + scale_fill_manual(values=cols) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  labs(title = "Principle Component Analysis (PCA) of \n S. siderea Physiology at T30",
       x = "PC1 (49.8% explained variance)",
       y = "PC2 (22.1% explained variance)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold")) 
rz_pca


#### Further exploring PCAs ####

#Note: I did not end up using any of the below analyses in the manuscript
#But, it is an interesting way to consider how individual coral physiology moves through time 

# Extract all results (coordinates, squared cosine, contributions) from individuals from PCA output
var1 <- get_pca_ind(facto3)
ssidt90coord <- var1$coord

#get coordinates for each individual by species and time point, then calculate distance moved between T30 and T90
PSTRt30coord <- cbind(pca_na_pstr_t30[,1:8],pstrt30coord)
PSTRt60coord <- cbind(pca_na_pstr_t60[,1:8],pstrt60coord)
PSTRt90coord <- cbind(pca_na_pstr_t90[,1:8],pstrt90coord)
#PSTRcoords <- rbind(PSTRt30coord, PSTRt60coord, PSTRt90coord)
#write.csv(PSTRcoords, "PCA_stuff/PSTR_PCA_coordinates.csv", row.names=FALSE)

#excluding individuals that don't have a gen_rz and treatment match between the two timepoints
PSTRt30coord_forcalc <- PSTRt30coord %>%
  group_by(rz, gen_rz, pco2,temp) %>%
  filter(id != "SFPX32") %>%
  filter(id != "SNPY7")
#write.csv(PSTRt30coord_forcalc, "PCA_stuff/PSTRt30coord_forcalc.csv")

PSTRt90coord_forcalc <- PSTRt90coord %>%
  group_by(rz, gen_rz, pco2,temp) %>%
  filter(id != "SFPX24") %>%
  filter(id != "SNPX10")
#write.csv(PSTRt90coord_forcalc, "PCA_stuff/PSTRt90coord_forcalc.csv")

#used excel to sort the t30 and t90coord_forcalc.csv outputs so they matched, then read back in here.
PSTR_t30coords <- read.csv("PCA_stuff/PSTRt30coord_forcalc.csv")
PSTR_t90coords <- read.csv("PCA_stuff/PSTRt90coord_forcalc.csv")

#to calculate distance moved by an individual from T30 to T90
PSTR_distances <- PSTR_t90coords[,2:7]
PSTR_distances$distances <- sqrt(((PSTR_t90coords$Dim.1-PSTR_t30coords$Dim.1)^2)+((PSTR_t90coords$Dim.2-PSTR_t30coords$Dim.2)^2))
head(PSTR_distances)

#to do centroid analysis
#calculate distance from center for each time point
PSTR_centroid <- PSTR_t30coords[,2:7]
PSTR_t30_x <- mean(PSTR_t30coords$Dim.1) #center x coordinate
PSTR_t30_y <- mean(PSTR_t30coords$Dim.2) #center y coordinate
PSTR_centroid$distcent_t30 <- sqrt(abs(((PSTR_t30_y-PSTR_t30coords$Dim.2)^2)-((PSTR_t30_x-PSTR_t30coords$Dim.1)^2)))

PSTR_t90_x <- mean(PSTR_t90coords$Dim.1) #center x coordinate T90
PSTR_t90_y <- mean(PSTR_t90coords$Dim.2) #center y coordinate T90
PSTR_centroid$distcent_t90 <- sqrt(abs(((PSTR_t90_y-PSTR_t90coords$Dim.2)^2)-((PSTR_t90_x-PSTR_t90coords$Dim.1)^2)))


#Now for SSID
SSIDt30coord <- cbind(pca_na_ssid_t30[,1:8],ssidt30coord)
SSIDt60coord <- cbind(pca_na_ssid_t60[,1:8],ssidt60coord)
SSIDt90coord <- cbind(pca_na_ssid_t90[,1:8],ssidt90coord)
SSIDcoords <- rbind(SSIDt30coord, SSIDt60coord, SSIDt90coord)
#write.csv(SSIDcoords, "PCA_stuff/SSID_PCA_coordinates.csv", row.names=FALSE)

SSIDt30coord_forcalc <- SSIDt30coord %>%
  group_by(rz,gen_rz,pco2,temp) %>%
  filter(id != "SFPZ18")
#write.csv(SSIDt30coord_forcalc, "PCA_stuff/SSIDt30coord_forcalc.csv")

SSIDt90coord_forcalc <- SSIDt90coord %>%
  group_by(rz, gen_rz, pco2,temp) %>%
  filter(id != "SFSY8") %>%
  filter(id != "SFSX10") %>%
  filter(id != "SFPZ24") %>%
  filter(id != "SNSY11") %>%
  filter(id != "SNSX22")
#write.csv(SSIDt90coord_forcalc, "PCA_stuff/SSIDt90coord_forcalc.csv")

#used excel to sort the t30 and t90coord_forcalc.csv outputs so they matched, then read back in here.
SSID_t30coords <- read.csv("PCA_stuff/SSIDt30coord_forcalc.csv")
SSID_t90coords <- read.csv("PCA_stuff/SSIDt90coord_forcalc.csv")

#calculate distance moved between T30 and T90
SSID_distances <- SSID_t90coords[,2:7]
SSID_distances$distances <- sqrt(((SSID_t90coords$Dim.1-SSID_t30coords$Dim.1)^2)+((SSID_t90coords$Dim.2-SSID_t30coords$Dim.2)^2))
head(SSID_distances)

#to do centroid analysis
#calculate distance from center for each time point
SSID_centroid <- SSID_t30coords[,2:7]
SSID_t30_x <- mean(SSID_t30coords$Dim.1) #center x coordinate
SSID_t30_y <- mean(SSID_t30coords$Dim.2) #center y coordinate
SSID_centroid$distcent_t30 <- sqrt(abs(((SSID_t30_y-SSID_t30coords$Dim.2)^2)-((SSID_t30_x-SSID_t30coords$Dim.1)^2)))

SSID_t90_x <- mean(SSID_t90coords$Dim.1) #center x coordinate T90
SSID_t90_y <- mean(SSID_t90coords$Dim.2) #center y coordinate T90
SSID_centroid$distcent_t90 <- sqrt(abs(((SSID_t90_y-SSID_t90coords$Dim.2)^2)-((SSID_t90_x-SSID_t90coords$Dim.1)^2)))

#Now plot
#First summarySE distance moved T30 to T90
distances <- rbind(PSTR_distances,SSID_distances)
head(distances)
levels(distances$treatment)
levels(distances$treatment) <- c("358uatm, 28°C","424uatm, 31°C","674uatm, 28°C","606uatm, 31°C","2750uatm, 28°C","2917uatm, 31°C")

distances.means <- summarySE(distances, measurevar="distances", groupvars=c("spp","rz","treatment"))

#Centroid summarySE
centroid <- rbind(PSTR_centroid,SSID_centroid)
head(centroid)
levels(centroid$treatment)
levels(centroid$treatment) <- c("358uatm, 28°C","424uatm, 31°C","674uatm, 28°C","606uatm, 31°C","2750uatm, 28°C","2917uatm, 31°C")

centroid_t30means <- summarySE(centroid, measurevar="distcent_t30", groupvars=c("spp", "rz"))
centroid_t90means <- summarySE(centroid, measurevar="distcent_t90", groupvars=c("spp", "rz"))


pos = position_dodge(width = 0.5)
p.distances <- ggplot(centroid, aes(x = spp, y = distcent_t90)) +
  geom_boxplot(aes(color=rz), outlier.shape=NA)+
  geom_jitter(aes(color=rz), alpha=0.5)+
  #geom_errorbar(aes(x = spp, ymax = distcent_t90+se, ymin = distcent_t90-se), position = pos, width = 0.2) +
  ggtitle("Distance from centroid T30")+
  ylab("Distance from Centroid")+
  xlab("Species")+
  #ylim(0,7)+  
  theme_cowplot()+
  scale_color_manual(values=c("F"="slategray4", "N"="black"), breaks = c("F","N", name = "reef zone"))+
  #scale_shape_manual(values=c("S"=17, "P"=19), breaks=c("S","P"), name = "species")+
  #facet_wrap(treatment~., nrow=3, ncol=2)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, 
        #legend.text.align = 0,
        legend.title = element_text(face = "bold")) 
p.distances
ggsave(file="PCA_stuff/centroid_T30.pdf", p.distances, width=5, height=6, units=c("in"), useDingbats=FALSE)

distanceaov <- aov(distances~treatment+rz, data = SSID_distances)
summary(distanceaov)

