############### Analysis of NAD-ciprofloxacin mixture growth plates  #####################


##################### Read in files and libraries#######################################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
library(growthcurver)
library(MetBrewer)
library(plotrix)
library(patchwork)
library(dunn.test)
library(MASS)
library(emmeans)
library(DHARMa)

######################## Make confidence interval function ##########################


conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}



#############################Read in the Data ################################

# list files
files <- list.files(path="data", full.names=TRUE)


# # create empty list to store results
output = list() #this has the raw dataframes

output2 = list() #this has the grouped dataframes

output3 = list()

# i=2

# run a for loop
for(i in 1:length(files)){
  cat(files[i])
  # read in file
  df <- read.table(files[i], header=TRUE, sep="\t")
  
  #Change names
  
  ### Rename all the columns
  names(df)[names(df)=="A1"] <- "LOEC_7.8_1"
  names(df)[names(df)=="A2"] <- "LOEC_7.8_2"
  names(df)[names(df)=="A3"] <- "LOEC_7.8_3"
  names(df)[names(df)=="A4"] <- "LOEC_7.8_4"
  names(df)[names(df)=="A5"] <- "NOEC_7.8_1"
  names(df)[names(df)=="A6"] <- "NOEC_7.8_2"
  names(df)[names(df)=="A7"] <- "NOEC_7.8_3"
  names(df)[names(df)=="A8"] <- "NOEC_7.8_4"
  names(df)[names(df)=="A9"] <- "NO_7.8_1"
  names(df)[names(df)=="A10"] <- "NO_7.8_2"
  names(df)[names(df)=="A11"] <- "NO_7.8_3"
  names(df)[names(df)=="A12"] <- "NO_7.8_4"
  
  #B
  names(df)[names(df)=="B1"] <- "LOEC_3.9_1"
  names(df)[names(df)=="B2"] <- "LOEC_3.9_2"
  names(df)[names(df)=="B3"] <- "LOEC_3.9_3"
  names(df)[names(df)=="B4"] <- "LOEC_3.9_4"
  names(df)[names(df)=="B5"] <- "NOEC_3.9_1"
  names(df)[names(df)=="B6"] <- "NOEC_3.9_2"
  names(df)[names(df)=="B7"] <- "NOEC_3.9_3"
  names(df)[names(df)=="B8"] <- "NOEC_3.9_4"
  names(df)[names(df)=="B9"] <- "NO_3.9_1"
  names(df)[names(df)=="B10"] <- "NO_3.9_2"
  names(df)[names(df)=="B11"] <- "NO_3.9_3"
  names(df)[names(df)=="B12"] <- "NO_3.9_4"
  
  names(df)[names(df)=="C1"] <- "LOEC_1.95_1"
  names(df)[names(df)=="C2"] <- "LOEC_1.95_2"
  names(df)[names(df)=="C3"] <- "LOEC_1.95_3"
  names(df)[names(df)=="C4"] <- "LOEC_1.95_4"
  names(df)[names(df)=="C5"] <- "NOEC_1.95_1"
  names(df)[names(df)=="C6"] <- "NOEC_1.95_2"
  names(df)[names(df)=="C7"] <- "NOEC_1.95_3"
  names(df)[names(df)=="C8"] <- "NOEC_1.95_4"
  names(df)[names(df)=="C9"] <- "NO_1.95_1"
  names(df)[names(df)=="C10"] <- "NO_1.95_2"
  names(df)[names(df)=="C11"] <- "NO_1.95_3"
  names(df)[names(df)=="C12"] <- "NO_1.95_4"
  
  names(df)[names(df)=="D1"] <- "LOEC_0.975_1"
  names(df)[names(df)=="D2"] <- "LOEC_0.975_2"
  names(df)[names(df)=="D3"] <- "LOEC_0.975_3"
  names(df)[names(df)=="D4"] <- "LOEC_0.975_4"
  names(df)[names(df)=="D5"] <- "NOEC_0.975_1"
  names(df)[names(df)=="D6"] <- "NOEC_0.975_2"
  names(df)[names(df)=="D7"] <- "NOEC_0.975_3"
  names(df)[names(df)=="D8"] <- "NOEC_0.975_4"
  names(df)[names(df)=="D9"] <- "NO_0.975_1"
  names(df)[names(df)=="D10"] <- "NO_0.975_2"
  names(df)[names(df)=="D11"] <- "NO_0.975_3"
  names(df)[names(df)=="D12"] <- "NO_0.975_4"
  
  names(df)[names(df)=="E1"] <- "LOEC_0.4875_1"
  names(df)[names(df)=="E2"] <- "LOEC_0.4875_2"
  names(df)[names(df)=="E3"] <- "LOEC_0.4875_3"
  names(df)[names(df)=="E4"] <- "LOEC_0.4875_4"
  names(df)[names(df)=="E5"] <- "NOEC_0.4875_1"
  names(df)[names(df)=="E6"] <- "NOEC_0.4875_2"
  names(df)[names(df)=="E7"] <- "NOEC_0.4875_3"
  names(df)[names(df)=="E8"] <- "NOEC_0.4875_4"
  names(df)[names(df)=="E9"] <- "NO_0.4875_1"
  names(df)[names(df)=="E10"] <- "NO_0.4875_2"
  names(df)[names(df)=="E11"] <- "NO_0.4875_3"
  names(df)[names(df)=="E12"] <- "NO_0.4875_4"
  
  names(df)[names(df)=="F1"] <- "LOEC_0.24375_1"
  names(df)[names(df)=="F2"] <- "LOEC_0.24375_2"
  names(df)[names(df)=="F3"] <- "LOEC_0.24375_3"
  names(df)[names(df)=="F4"] <- "LOEC_0.24375_4"
  names(df)[names(df)=="F5"] <- "NOEC_0.24375_1"
  names(df)[names(df)=="F6"] <- "NOEC_0.24375_2"
  names(df)[names(df)=="F7"] <- "NOEC_0.24375_3"
  names(df)[names(df)=="F8"] <- "NOEC_0.24375_4"
  names(df)[names(df)=="F9"] <- "NO_0.24375_1"
  names(df)[names(df)=="F10"] <- "NO_0.24375_2"
  names(df)[names(df)=="F11"] <- "NO_0.24375_3"
  names(df)[names(df)=="F12"] <- "NO_0.24375_4"
  
  names(df)[names(df)=="G1"] <- "LOEC_0.121875_1"
  names(df)[names(df)=="G2"] <- "LOEC_0.121875_2"
  names(df)[names(df)=="G3"] <- "LOEC_0.121875_3"
  names(df)[names(df)=="G4"] <- "LOEC_0.121875_4"
  names(df)[names(df)=="G5"] <- "NOEC_0.121875_1"
  names(df)[names(df)=="G6"] <- "NOEC_0.121875_2"
  names(df)[names(df)=="G7"] <- "NOEC_0.121875_3"
  names(df)[names(df)=="G8"] <- "NOEC_0.121875_4"
  names(df)[names(df)=="G9"] <- "NO_0.121875_1"
  names(df)[names(df)=="G10"] <- "NO_0.121875_2"
  names(df)[names(df)=="G11"] <- "NO_0.121875_3"
  names(df)[names(df)=="G12"] <- "NO_0.121875_4"
  
  names(df)[names(df)=="H1"] <- "LOEC_0_1"
  names(df)[names(df)=="H2"] <- "LOEC_0_2"
  names(df)[names(df)=="H3"] <- "LOEC_0_3"
  names(df)[names(df)=="H4"] <- "LOEC_0_4" 
  names(df)[names(df)=="H5"] <- "NOEC_0_1"
  names(df)[names(df)=="H6"] <- "NOEC_0_2"
  names(df)[names(df)=="H7"] <- "NOEC_0_3"
  names(df)[names(df)=="H8"] <- "NOEC_0_4"
  names(df)[names(df)=="H9"] <- "NO_0_1"
  names(df)[names(df)=="H10"] <- "NO_0_2"
  names(df)[names(df)=="H11"] <- "NO_0_3"
  names(df)[names(df)=="H12"] <- "NO_0_4"
  
  
  #Change the time into minutes/hours
  names(df)[names(df)=="Kinetic.read"] <- "time"
  
  #Change them to individual 10 min intervals 
  df$time <- seq(0,1440, by=10)
  
  #Melt dataframe
  df_melt <- melt(df, id=c("time"))
  
  #Separate
  df_melt <- separate(df_melt, variable, c("Treatment", "Cip_Concentration", "Rep"), sep="_")
  
  # need to add a column so you know which file and week this is
  df_melt <- mutate(df_melt, file = basename(tools::file_path_sans_ext(files[i])))
  
  
  #New for loop with these
  df_stats <- df_melt %>%
    group_by(Treatment, Cip_Concentration, time) %>%
    summarise(N=length(value),
              Average=mean(value),
              CI95=conf_int95(value),
              se=std.error(value))
  
  
  
  # probably need to add a column so you know which file this is
  df_stats <- mutate(df_stats, file = basename(tools::file_path_sans_ext(files[i])))
  
  
  
  # then add this iteration to our empty list
  output[[i]] <- df_melt
  
  # then add this iteration to our empty list
  output2[[i]] <- df_stats
  
  output3[[i]] <- df
  
  rm(df)
  rm(df_melt)
  rm(df_stats)
  
}

## Make one dataframe!


## Total Dfs in one!
df_melt_total <- plyr::ldply(output, data.frame)

df_stats_total <- plyr::ldply(output2, data.frame)

df_diclofenac <- plyr::ldply(output3[1], data.frame)

df_metformin <- plyr::ldply(output3[3], data.frame)

df_estradiol <- plyr::ldply(output3[2], data.frame)



######################## SELECT Working Diclofenac ##########################

df_SELECT_diclofenac <- df_diclofenac
rownames(df_SELECT_diclofenac) <- df_SELECT_diclofenac$time
df_SELECT_diclofenac <- as.data.frame(t(df_SELECT_diclofenac))
df_SELECT_diclofenac <- df_SELECT_diclofenac[-c(1), ]

#Add extra column of conc
df_SELECT_diclofenac$conc <- row.names(df_SELECT_diclofenac)
df_SELECT_diclofenac <- df_SELECT_diclofenac %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_diclofenac$Concentration <- as.numeric(df_SELECT_diclofenac$Concentration)


#Change column names of hour time points
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="0"] <- "t0"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="60"] <- "t60"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="120"] <- "t120"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="180"] <- "t180"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="240"] <- "t240"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="300"] <- "t300"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="360"] <- "t360"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="420"] <- "t420"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="480"] <- "t480"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="540"] <- "t540"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="600"] <- "t600"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="660"] <- "t660"
names(df_SELECT_diclofenac)[names(df_SELECT_diclofenac)=="720"] <- "t720"

#Select only the hour time points
df_SELECT_diclofenac <- df_SELECT_diclofenac %>% dplyr::select(starts_with("t"))

#Control 
df_SELECT_control_diclofenac <- subset(df_SELECT_diclofenac, Treatment=="NO")

#LOEC
df_SELECT_LOEC_diclofenac <- subset(df_SELECT_diclofenac, Treatment== "LOEC")

#This gives us concentration and treatment columns
df_SELECT_control_diclofenac$conc <- row.names(df_SELECT_control_diclofenac)
df_SELECT_control_diclofenac<- df_SELECT_control_diclofenac %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_control_diclofenac$Concentration <- as.numeric(df_SELECT_control_diclofenac$Concentration)

df_SELECT_control_diclofenac$conctreat <- paste(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$Treatment, sep = "_")

df_SELECT_nodrugcontrol_diclofenac <- subset(df_SELECT_control_diclofenac, Concentration=="0")

#This gives us concentration and treatment columns
df_SELECT_LOEC_diclofenac$conc <- row.names(df_SELECT_LOEC_diclofenac)
df_SELECT_LOEC_diclofenac<- df_SELECT_LOEC_diclofenac %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_LOEC_diclofenac$Concentration <- as.numeric(df_SELECT_LOEC_diclofenac$Concentration)

df_SELECT_LOEC_diclofenac$conctreat <- paste(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$Treatment, sep = "_")

#Add no Drug control to LOEC dataframe
df_SELECT_LOEC_diclofenac <- rbind(df_SELECT_LOEC_diclofenac, df_SELECT_control_diclofenac)


#Test which timepoint in hours (but technically in minutes) has the largest difference in the 
#If p>0,05 not normal
#Normal correlation = Pearsons
#Not normal correlation = Spearman

#Hour 4
shapiro.test(df_SELECT_LOEC_diclofenac$t240)
#p-value = 0.003593

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t240, method="spearman")
#p-value = 0.001532


#Hour 5
shapiro.test(df_SELECT_LOEC_diclofenac$t300)
#p-value = 0.008098

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t300, method="spearman")
#p-value = 0.0001508


#Hour 6
shapiro.test(df_SELECT_LOEC_diclofenac$t360)
#p-value = 0..02979


#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t360, method="spearman")
#p-value = 5.696e-07


#Hour 7
shapiro.test(df_SELECT_LOEC_diclofenac$t420)
#p-value = 0.1468

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t420, method="pearson")
#p-value = 3.696e-05


#Hour 8
shapiro.test(df_SELECT_LOEC_diclofenac$t480)
#p-value = 0.04747

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t480, method="spearman")
#p-value = 1.496e-07


#Hour 9
shapiro.test(df_SELECT_LOEC_diclofenac$t540)
#p-value = 0.04362

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t540, method="spearman")
#p-value = 9.095e-05


#Hour 10
shapiro.test(df_SELECT_LOEC_diclofenac$t600)
#p-value = 0.1892

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t600, method="pearson")
#p-value = 0.005399


#Hour 11
shapiro.test(df_SELECT_LOEC_diclofenac$t660)
#p-value = 0.09085

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t660, method="pearson")
#p-value = 0.008606


#Hour 12
shapiro.test(df_SELECT_LOEC_diclofenac$t720)
#p-value = 0.3857

#Not normal
cor.test(df_SELECT_LOEC_diclofenac$Concentration, df_SELECT_LOEC_diclofenac$t720, method="pearson")
#p-value = 0.02825


#Hour 8 - 480 has the greatest dose response
#TEst which concentrations are significantly different to 0
dunn.test(df_SELECT_LOEC_diclofenac$t480, df_SELECT_LOEC_diclofenac$conctreat, list = TRUE)

# The LOEC of cip in this mixture is 0.4875ug/L p=0.0229



################# Ciprofloxacin LOEC Working (No Diclofenac) 

#Test which timepoint in hours (but technically in minutes) has the largest difference in the
#If p>0,05 not normal
#Normal correlation = Pearsons
#Not normal correlation = Spearman

#Hour 4
shapiro.test(df_SELECT_control_diclofenac_diclofenac$t240)
#p-value = 3.56e-06

#Not normal
cor.test(df_SELECT_control_diclofenac_diclofenac$Concentration, df_SELECT_control_diclofenac_diclofenac$t240, method="spearman")
#p-value = 0.3617


#Hour 5
shapiro.test(df_SELECT_control_diclofenac$t300)
#p-value = 0.001825

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t300, method="spearman")
#p-value = 0.2898


#Hour 6
shapiro.test(df_SELECT_control_diclofenac$t360)
#p-value = 0.06814

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t360, method="spearman")
#p-value = 0.0006806


#Hour 7
shapiro.test(df_SELECT_control_diclofenac$t420)
#p-value = 0.0006806

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t420, method="spearman")
#p-value = 4.383e-06


#Hour 8
shapiro.test(df_SELECT_control_diclofenac$t480)
#p-value = 0.002248

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t480, method="spearman")
#p-value =2.034e-05


#Hour 9
shapiro.test(df_SELECT_control_diclofenac$t540)
#p-value = 1.576e-05

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t540, method="spearman")
#p-value = 0.06621


#Hour 10
shapiro.test(df_SELECT_control_diclofenac$t600)
#p-value = 8.777e-06

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t600, method="spearman")
#p-value = 0.02101

#Hour 11
shapiro.test(df_SELECT_control_diclofenac$t660)
#p-value = 0.0003342

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t660, method="spearman")
#p-value = 0.002652


#Hour 12
shapiro.test(df_SELECT_control_diclofenac$t720)
#p-value = 0.1022

#Not normal
cor.test(df_SELECT_control_diclofenac$Concentration, df_SELECT_control_diclofenac$t720, method="pearson")
#p-value = 0.05852


#Hour 7 has the greatest dose response
#TEst which concentrations are significantly different to 0
dunn.test(df_SELECT_control_diclofenac$t420, df_SELECT_control_diclofenac$Concentration, list = TRUE)


#LOEC of ciprofloxacin in the diclofenac mixture is 3.7ug/L 


################# SELECT Working Metformin ##############################

df_SELECT_metformin <- df_metformin
rownames(df_SELECT_metformin) <- df_SELECT_metformin$time
df_SELECT_metformin <- as.data.frame(t(df_SELECT_metformin))
df_SELECT_metformin <- df_SELECT_metformin[-c(1), ]

#Add extra column of conc
df_SELECT_metformin$conc <- row.names(df_SELECT_metformin)
df_SELECT_metformin <- df_SELECT_metformin %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_metformin$Concentration <- as.numeric(df_SELECT_metformin$Concentration)


#Change column names of hour time points
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="0"] <- "t0"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="60"] <- "t60"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="120"] <- "t120"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="180"] <- "t180"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="240"] <- "t240"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="300"] <- "t300"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="360"] <- "t360"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="420"] <- "t420"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="480"] <- "t480"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="540"] <- "t540"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="600"] <- "t600"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="660"] <- "t660"
names(df_SELECT_metformin)[names(df_SELECT_metformin)=="720"] <- "t720"

#Select only the hour time points
df_SELECT_metformin <- df_SELECT_metformin %>% dplyr::select(starts_with("t"))

#Control 
df_SELECT_control_metformin <- subset(df_SELECT_metformin, Treatment=="NO")

#LOEC
df_SELECT_LOEC_metformin <- subset(df_SELECT_metformin, Treatment== "LOEC")

#This gives us concentration and treatment columns
df_SELECT_control_metformin$conc <- row.names(df_SELECT_control_metformin)
df_SELECT_control_metformin<- df_SELECT_control_metformin %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_control_metformin$Concentration <- as.numeric(df_SELECT_control_metformin$Concentration)

df_SELECT_control_metformin$conctreat <- paste(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$Treatment, sep = "_")

df_SELECT_nodrugcontrol_metformin <- subset(df_SELECT_control_metformin, Concentration=="0")

#This gives us concentration and treatment columns
df_SELECT_LOEC_metformin$conc <- row.names(df_SELECT_LOEC_metformin)
df_SELECT_LOEC_metformin<- df_SELECT_LOEC_metformin %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_LOEC_metformin$Concentration <- as.numeric(df_SELECT_LOEC_metformin$Concentration)

df_SELECT_LOEC_metformin$conctreat <- paste(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$Treatment, sep = "_")

#Add no Drug control to LOEC dataframe
df_SELECT_LOEC_metformin <- rbind(df_SELECT_LOEC_metformin, df_SELECT_nodrugcontrol_metformin)


#Test which timepoint in hours (but technically in minutes) has the largest difference in the 
#If p>0,05 not normal
#Normal correlation = Pearsons
#Not normal correlation = Spearman

#Hour 4
shapiro.test(df_SELECT_LOEC_metformin$t240)
#p-value = 0.0005633

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t240, method="spearman")
#p-value = 0.002464


#Hour 5
shapiro.test(df_SELECT_LOEC_metformin_metformin$t300)
#p-value = 0.008668

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t300, method="spearman")
#p-value = 0.0005553


#Hour 6
shapiro.test(df_SELECT_LOEC_metformin$t360)
#p-value = 0.3593


#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t360, method="pearson")
#p-value = 0.0002765


#Hour 7
shapiro.test(df_SELECT_LOEC_metformin$t420)
#p-value = 0.7626

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t420, method="pearson")
#p-value = 6.515e-05


#Hour 8
shapiro.test(df_SELECT_LOEC_metformin$t480)
#p-value = 0.5847

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t480, method="pearson")
#p-value = 0.0007877


#Hour 9
shapiro.test(df_SELECT_LOEC_metformin$t540)
#p-value = 0.1774

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t540, method="pearson")
#p-value = 0.5117



#Not doing from here because the graph looks weird - the growth goes down lol
#Hour 10
shapiro.test(df_SELECT_LOEC_metformin$t600)
#p-value = 0.0428

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t600, method="spearman")
#p-value = 4.86e-05


#Hour 11
shapiro.test(df_SELECT_LOEC_metformin$t660)
#p-value = 0.09183

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t660, method="pearson")
#p-value = 0.003782


#Hour 12
shapiro.test(df_SELECT_LOEC_metformin$t720)
#p-value = 0.3527

#Not normal
cor.test(df_SELECT_LOEC_metformin$Concentration, df_SELECT_LOEC_metformin$t720, method="spearman")
#p-value = 3.442e-05


#Hour 7 - 420 has the greatest dose response
#TEst which concentrations are significantly different to 0
dunn.test(df_SELECT_LOEC_metformin$t420, df_SELECT_LOEC_metformin$conctreat, list = TRUE)

#With the mixture, (so compared to the control with NO cip, we have cip and metformin LOEC = 0.24375ug/L)



########### Ciprofloxacin LOEC Working (No Metformin)

#Test which timepoint in hours (but technically in minutes) has the largest difference in the
#If p>0,05 not normal
#Normal correlation = Pearsons
#Not normal correlation = Spearman

#Hour 4
shapiro.test(df_SELECT_control_metformin$t240)
#p-value = 0.1173

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t240, method="pearson")
#p-value = 0.0008689


#Hour 5
shapiro.test(df_SELECT_control_metformin_metformin$t300)
#p-value = 0.3954

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t300, method="spearman")
#p-value = 4.932e-10


#Hour 6
shapiro.test(df_SELECT_control_metformin$t360)
#p-value = 0.9251

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t360, method="pearson")
#p-value = 3.058e-09


#Hour 7
shapiro.test(df_SELECT_control_metformin$t420)
#p-value = 0.8895

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t420, method="pearson")
#p-value = 6.01e-07


#Hour 8
shapiro.test(df_SELECT_control_metformin$t480)
#p-value = 0.6185

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t480, method="pearson")
#p-value = 0.0001008


#Hour 9
shapiro.test(df_SELECT_control_metformin$t540)
#p-value = 0.07692

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t540, method="pearson")
#p-value = 0.6971

#Not from herebecause the growth is weird
#Hour 10
shapiro.test(df_SELECT_control_metformin$t600)
#p-value = 0.02455

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t600, method="spearman")
#p-value = 0.0009638

#Hour 11
shapiro.test(df_SELECT_control_metformin$t660)
#p-value = 2.435e-06

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t660, method="spearman")
#p-value = 0.0004947


#Hour 12
shapiro.test(df_SELECT_control_metformin$t720)
#p-value = 9.142e-08

#Not normal
cor.test(df_SELECT_control_metformin$Concentration, df_SELECT_control_metformin$t720, method="pearson")
#p-value = 0.493


#Hour 5 has the greatest dose response
#TEst which concentrations are significantly different to 0
dunn.test(df_SELECT_control_metformin$t300, df_SELECT_control_metformin$Concentration)


################# SELECT Working 17-beta-estradiol ##############################

df_SELECT_estradiol <- df_estradiol
rownames(df_SELECT_estradiol) <- df_SELECT_estradiol$time
df_SELECT_estradiol <- as.data.frame(t(df_SELECT_estradiol))
df_SELECT_estradiol <- df_SELECT_estradiol[-c(1), ]

#Add extra column of conc
df_SELECT_estradiol$conc <- row.names(df_SELECT_estradiol)
df_SELECT_estradiol <- df_SELECT_estradiol %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_estradiol$Concentration <- as.numeric(df_SELECT_estradiol$Concentration)


#Change column names of hour time points
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="0"] <- "t0"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="60"] <- "t60"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="120"] <- "t120"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="180"] <- "t180"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="240"] <- "t240"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="300"] <- "t300"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="360"] <- "t360"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="420"] <- "t420"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="480"] <- "t480"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="540"] <- "t540"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="600"] <- "t600"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="660"] <- "t660"
names(df_SELECT_estradiol)[names(df_SELECT_estradiol)=="720"] <- "t720"

#Select only the hour time points
df_SELECT_estradiol <- df_SELECT_estradiol %>% dplyr::select(starts_with("t"))

#Control 
df_SELECT_control_estradiol <- subset(df_SELECT_estradiol, Treatment=="NO")

#LOEC
df_SELECT_LOEC_estradiol <- subset(df_SELECT_estradiol, Treatment== "LOEC")

#This gives us concentration and treatment columns
df_SELECT_control_estradiol$conc <- row.names(df_SELECT_control_estradiol)
df_SELECT_control_estradiol<- df_SELECT_control_estradiol %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_control_estradiol$Concentration <- as.numeric(df_SELECT_control_estradiol$Concentration)

df_SELECT_control_estradiol$conctreat <- paste(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$Treatment, sep = "_")

df_SELECT_nodrugcontrol_estradiol <- subset(df_SELECT_control_estradiol, Concentration=="0")

#This gives us concentration and treatment columns
df_SELECT_LOEC_estradiol$conc <- row.names(df_SELECT_LOEC_estradiol)
df_SELECT_LOEC_estradiol<- df_SELECT_LOEC_estradiol %>% separate(conc, c("Treatment", "Concentration", "Rep"), sep="_")
df_SELECT_LOEC_estradiol$Concentration <- as.numeric(df_SELECT_LOEC_estradiol$Concentration)

df_SELECT_LOEC_estradiol$conctreat <- paste(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$Treatment, sep = "_")

#Add no Drug control to LOEC dataframe
df_SELECT_LOEC_estradiol <- rbind(df_SELECT_LOEC_estradiol, df_SELECT_nodrugcontrol_estradiol)

# df_SELECT_LOEC_working <- subset(df_SELECT_LOEC, conctreat != "0_LOEC")

#Test which timepoint in hours (but technically in minutes) has the largest difference in the 
#If p>0,05 not normal
#Normal correlation = Pearsons
#Not normal correlation = Spearman

#Hour 4
shapiro.test(df_SELECT_LOEC_estradiol$t240)
#p-value = 0.0002781

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t240, method="spearman")
#p-value = 0.001436


#Hour 5
shapiro.test(df_SELECT_LOEC_estradiol$t300)
#p-value = 0.0002252

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t300, method="spearman")
#p-value = 0.001753


#Hour 6
shapiro.test(df_SELECT_LOEC_estradiol$t360)
#p-value = 0.0002307


#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t360, method="spearman")
#p-value = 0.0009638


#Hour 7
shapiro.test(df_SELECT_LOEC_estradiol$t420)
#p-value = 0.00218

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t420, method="spearman")
#p-value = 0.0001029


#Hour 8
shapiro.test(df_SELECT_LOEC_estradiol$t480)
#p-value = 0.0523

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t480, method="pearson")
#p-value =0.002527


#Hour 9
shapiro.test(df_SELECT_LOEC_estradiol$t540)
#p-value = 0.1619

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t540, method="pearson")
#p-value = 0.005387

#Hour 10
shapiro.test(df_SELECT_LOEC_estradiol$t600)
#p-value = 0.04925

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t600, method="spearman")
#p-value = 0.00106


#Hour 11
shapiro.test(df_SELECT_LOEC_estradiol$t660)
#p-value = 0.1009

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t660, method="pearson")
#p-value = 0.06393


#Hour 12
shapiro.test(df_SELECT_LOEC_estradiol$t720)
#p-value = 0.1849

#Not normal
cor.test(df_SELECT_LOEC_estradiol$Concentration, df_SELECT_LOEC_estradiol$t720, method="pearson")
#p-value = 0.04112


#Hour 8 - 5480 has the greatest dose response
#TEst which concentrations are significantly different to 0
dunn.test(df_SELECT_LOEC_estradiol$t360, df_SELECT_LOEC_estradiol$conctreat, list = TRUE)


################### Ciprofloxacin LOEC Working (No 17-Beta) 

#Test which timepoint in hours (but technically in minutes) has the largest difference in the
#If p>0,05 not normal
#Normal correlation = Pearsons
#Not normal correlation = Spearman

#Hour 4
shapiro.test(df_SELECT_control_estradiol_estradiol$t240)
#p-value = 0.0001402

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t240, method="pearson")
#p-value = 0.0001402


#Hour 5
shapiro.test(df_SELECT_control_estradiol$t300)
#p-value = 0.107

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t300, method="pearson")
#p-value = 5.383e-05


#Hour 6
shapiro.test(df_SELECT_control_estradiol$t360)
#p-value = 0.20228

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t360, method="pearson")
#p-value = 7.392e-06


#Hour 7
shapiro.test(df_SELECT_control_estradiol$t420)
#p-value = 0.06542

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t420, method="pearson")
#p-value = 4.682e-09


#Hour 8
shapiro.test(df_SELECT_control_estradiol$t480)
#p-value = 0.008131

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t480, method="spearman")
#p-value = 1.577e-05


#Hour 9
shapiro.test(df_SELECT_control_estradiol$t540)
#p-value = 0.0005904

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t540, method="spearman")
#p-value = 5.921e-06


#Hour 10
shapiro.test(df_SELECT_control_estradiol$t600)
#p-value = 1.398e-05

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t600, method="spearman")
#p-value = 0.001134

#Hour 11
shapiro.test(df_SELECT_control_estradiol$t660)
#p-value = 0.1011

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t660, method="pearson")
#p-value = 0.0002964


#Hour 12
shapiro.test(df_SELECT_control_estradiol$t720)
#p-value = 0.7937

#Not normal
cor.test(df_SELECT_control_estradiol$Concentration, df_SELECT_control_estradiol$t720, method="pearson")
#p-value = 0.0002643

#Hour 7has the greatest dose response
#Test which concentrations are significantly different to 0
dunn.test(df_SELECT_control_estradiol$t600, df_SELECT_control_estradiol$Concentration, list = TRUE)

#LOEC is 3.9ug/L


######################### AUC Diclofenac ###################################

diclofenac_expo <- subset(df_diclofenac, time>0)
diclofenac_expo <- subset(diclofenac_expo, time<720)


diclofenac_expo <- SummarizeGrowthByPlate(diclofenac_expo, plot_fit = TRUE, plot_file = "GrowthCurvesPlateDiclofenac_xpo.pdf")

# Here we are using area under the exponential curve to determine overall changes in growth capacity

diclofenac_expo <- diclofenac_expo %>% separate(sample, c("Treatment", "Concentration", "Rep"), sep="_")

diclofenac_auc <- diclofenac_expo[c("Treatment", "Concentration", "Rep", "auc_e")]

# make an averaged dataframe 
diclofenac_xpo_avg <- diclofenac_auc %>%
  group_by(Concentration, Treatment) %>%
  summarise(
    mean=mean(auc_e),
    sd=sd(auc_e),
    se=std.error(auc_e)) %>%
  mutate(sd=sd)


diclofenac_expo$Concentration <- as.numeric(diclofenac_expo$Concentration)

#Make a model for the auc for xpo phase 
diclo_model1 <- lm(auc_e ~ Concentration * Treatment, data=diclofenac_expo)

summary(diclo_model1)

drop1(diclo_model1, test="Chisq")  
#interaction term is not significant 

diclo_model2 <- lm(auc_e ~ Concentration + Treatment, data=diclofenac_expo)
summary(diclo_model2)

drop1(diclo_model2, test="Chisq")
#Both are significant here

anova(diclo_model2, test="F")

#DHARMa stuff 
simulationOutput <- simulateResiduals(fittedModel = diclo_model2, plot = T)

ggplot(data = diclofenac_expo, aes(x = diclo_model2$residuals)) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

#Pairwise comparison 
emmeans(model2, specs = pairwise ~ Treatment|Concentration, type="response")

# $contrasts
# Concentration = 1.93:
#   contrast    estimate   SE df t.ratio p.value
# LOEC - NO     -54.16 6.25 92  -8.671  <.0001
# LOEC - NOEC   -52.08 6.25 92  -8.338  <.0001
# NO - NOEC       2.08 6.25 92   0.333  0.9406


### Make the Plot ###
# make an average for the mean points on the plot

#Make sure that gc_auc$Concentration is a numeric
diclofenac_xpo_avg$Concentration <- as.numeric(diclofenac_xpo_avg$Concentration)
diclofenac_xpo_avg$Treatment <- as.factor(diclofenac_xpo_avg$Treatment)


#Change the name of the mean column so that it matches the auc_e for the graph in gc_auc_avg
colnames(diclofenac_xpo_avg)[which(names(diclofenac_xpo_avg)=="mean")] <- "auc_e"


diclofenac_expo$TreatRep <- paste(diclofenac_expo$Treatment, diclofenac_expo$Rep)

diclofenac_expo$Treatment <- diclofenac_expo$Treatment %>% factor(levels =c("NO", "LOEC", "NOEC"))

# plot the graph 
diclofenac_line <- ggplot() +
  geom_point(data=diclofenac_expo, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, fill=Treatment), size=4, alpha=0.4)+
  geom_line(data=diclofenac_expo, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, group= TreatRep),
            linewidth=1.22, alpha = 0.1)+
  geom_point(data=diclofenac_xpo_avg, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, fill=Treatment), size=6)+
  geom_line(data=diclofenac_xpo_avg, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, group= Treatment), 
            linewidth=1.22)+
  scale_color_manual(values = c("skyblue4", "orange2", "rosybrown3" ), name = "Treatment",
                     labels = c("Ciprofloxacin Only", "High NAD mixture",
                                "Low NAD mixture"))+
  scale_fill_manual(values = c("skyblue4", "orange2", "rosybrown3"), name = "Treatment",
                    labels = c("Ciprofloxacin Only", "High NAD mixture",
                               "Low NAD mixture"))+
  labs(y="Area under the curve", x="Ciprofloxacin Concentration (µg/L)", 
       title = "A) Diclofenac Mixtures")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = 'none')  +
  theme(text=element_text(size=20))

diclofenac_line

######################### AUC Metformin #####################################

metformin_expo <- subset(df_metformin, time>0)
metformin_expo <- subset(metformin_expo, time<720)


metformin_expo <- SummarizeGrowthByPlate(metformin_expo, plot_fit = TRUE, plot_file = "GrowthCurvesPlateMetformin_xpo.pdf")

# Here we are using area under the exponential curve to determine overall changes in growth capacity

metformin_expo <- metformin_expo %>% separate(sample, c("Treatment", "Concentration", "Rep"), sep="_")

metformin_auc <- metformin_expo[c("Treatment", "Concentration", "Rep", "auc_e")]

metformin_xpo_avg <- metformin_auc %>%
  group_by(Concentration, Treatment) %>%
  summarise(
    mean=mean(auc_e),
    sd=sd(auc_e),
    se=std.error(auc_e)) %>%
  mutate(sd=sd)


metformin_auc$Concentration <- as.numeric(metformin_auc$Concentration)

#MAke a model for the auc for xpo phase 
metformin_model1 <- lm(auc_e ~ Concentration * Treatment, data=metformin_auc)

summary(metformin_model1)

drop1(metformin_model1, test="Chisq")  
#interaction term is not significant

metformin_model2 <- lm(auc_e ~ Concentration + Treatment, data=metformin_auc)
summary(metformin_model2)

drop1(metformin_model2, test="Chisq")
#Both are significant here

anova(metformin_model2, test="F")

#DHARMa stuff 
simulationOutput <- simulateResiduals(fittedModel = metformin_model2, plot = T)

ggplot(data = metformin_auc, aes(x = metformin_model2$residuals)) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


#Pairwise comparison thing
emmeans(metformin_model2, specs = pairwise ~ Treatment|Concentration, type="response")

# $contrasts
# Concentration = 1.93:
#   contrast    estimate  SE df t.ratio p.value
# LOEC - NO     -62.67 6.1 92 -10.278  <.0001
# LOEC - NOEC   -56.83 6.1 92  -9.319  <.0001
# NO - NOEC       5.84 6.1 92   0.958  0.6049

### Make the Plot ###

metformin_auc$Treatment <- metformin_auc$Treatment %>% factor(levels =c("NO", "LOEC", "NOEC"))


#Make sure that gc_auc$Concentration is a numeric
metformin_xpo_avg$Concentration <- as.numeric(metformin_xpo_avg$Concentration)
metformin_xpo_avg$Treatment <- as.factor(metformin_xpo_avg$Treatment)


#Change the name of the mean column so that it matches the auc_e for the graph in gc_auc_avg
colnames(metformin_xpo_avg)[which(names(metformin_xpo_avg)=="mean")] <- "auc_e"


metformin_auc$TreatRep <- paste(metformin_auc$Treatment, metformin_auc$Rep)


metformin_line <- ggplot() +
  geom_point(data=metformin_auc, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, fill=Treatment), size=4, alpha=0.4)+
  geom_line(data=metformin_auc, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, group= TreatRep),
            linewidth=1.22, alpha = 0.1)+
  geom_point(data=metformin_xpo_avg, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, fill=Treatment), size=6)+
  geom_line(data=metformin_xpo_avg, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, group= Treatment), 
            linewidth=1.22)+
  scale_color_manual(values = c("skyblue4", "orange2", "rosybrown3" ), name = "Treatment",
                     labels = c("Ciprofloxacin Only", "High NAD mixture",
                                "Low NAD mixture"))+
  scale_fill_manual(values = c("skyblue4", "orange2", "rosybrown3"), name = "Treatment",
                    labels = c("Ciprofloxacin Only", "High NAD mixture",
                               "Low NAD mixture"))+
  labs(y="Area under the curve", x="Ciprofloxacin Concentration (µg/L)", 
       title = "B) Metformin Mixtures")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = 'bottom')  +
  theme(text=element_text(size=20))


metformin_line


##################### AUC 17-beta-estradiol #####################################

estradiol_expo <- subset(df_estradiol, time>0)
estradiol_expo <- subset(estradiol_expo, time<720)


estradiol_expo <- SummarizeGrowthByPlate(estradiol_expo, plot_fit = TRUE, plot_file = "GrowthCurvesPlateEstradiol_xpo.pdf")

# Here we are using area under the exponential curve to determine overall changes in growth capacity

estradiol_expo <- estradiol_expo %>% separate(sample, c("Treatment", "Concentration", "Rep"), sep="_")

estradiol_auc <- estradiol_expo[c("Treatment", "Concentration", "Rep", "auc_e")]

estradiol_xpo_avg <- estradiol_auc %>%
  group_by(Concentration, Treatment) %>%
  summarise(
    mean=mean(auc_e),
    sd=sd(auc_e),
    se=std.error(auc_e)) %>%
  mutate(sd=sd)


estradiol_auc$Concentration <- as.numeric(estradiol_auc$Concentration)

#MAke a model for the auc for xpo phase 
estradiol_model1 <- lm(auc_e ~ Concentration * Treatment, data=estradiol_auc)

summary(estradiol_model1)

drop1(estradiol_model1, test="Chisq")  
#interaction term is not significant

estradiol_model2 <- lm(auc_e ~ Concentration + Treatment, data=estradiol_auc)
summary(estradiol_model2)

drop1(estradiol_model2, test="Chisq")
#Both are significant here

anova(estradiol_model2, test="F")

#DHARMa stuff 
simulationOutput <- simulateResiduals(fittedModel = estradiol_model2, plot = T)

ggplot(data = estradiol_auc, aes(x = estradiol_model2$residuals)) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


#Pairwise comparison thing
emmeans(estradiol_model2, specs = pairwise ~ Treatment|Concentration, type="response")

# $contrasts
# Concentration = 1.93:
#   contrast    estimate   SE df t.ratio p.value
# LOEC - NO     -71.38 8.09 92  -8.828  <.0001
# LOEC - NOEC   -66.01 8.09 92  -8.165  <.0001
# NO - NOEC       5.36 8.09 92   0.663  0.7852



### Make the Plot ###

estradiol_auc$Treatment <- estradiol_auc$Treatment %>% factor(levels =c("NO", "LOEC", "NOEC"))


#Make sure that gc_auc$Concentration is a numeric
estradiol_xpo_avg$Concentration <- as.numeric(estradiol_xpo_avg$Concentration)
estradiol_xpo_avg$Treatment <- as.factor(estradiol_xpo_avg$Treatment)
estradiol_xpo_avg$Treatment <- estradiol_xpo_avg$Treatment %>% factor(levels =c("NO", "LOEC", "NOEC"))


#Change the name of the mean column so that it matches the auc_e for the graph in gc_auc_avg
colnames(estradiol_xpo_avg)[which(names(estradiol_xpo_avg)=="mean")] <- "auc_e"


estradiol_auc$TreatRep <- paste(estradiol_auc$Treatment, estradiol_auc$Rep)

# plot it 
estradiol_line <- ggplot() +
  geom_point(data=estradiol_auc, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, fill=Treatment), size=4, alpha=0.4)+
  geom_line(data=estradiol_auc, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, group= TreatRep),
            linewidth=1.22, alpha = 0.1)+
  geom_point(data=estradiol_xpo_avg, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, fill=Treatment), size=6)+
  geom_line(data=estradiol_xpo_avg, mapping = aes(x=Concentration, y=auc_e, colour=Treatment, group= Treatment), 
            linewidth=1.22)+
  scale_color_manual(values = c("skyblue4", "orange2", "rosybrown3" ), name = "Treatment",
                     labels = c("Ciprofloxacin Only", "High NAD mixture",
                                "Low NAD mixture"))+
  scale_fill_manual(values = c("skyblue4", "orange2", "rosybrown3"), name = "Treatment",
                    labels = c("Ciprofloxacin Only", "High NAD mixture",
                               "Low NAD mixture"))+
  labs(y="Area under the curve", x="Ciprofloxacin Concentration (µg/L)", 
       title = "C) 17-β-estradiol Mixtures")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = 'none')  +
  theme(text=element_text(size=20))


estradiol_line


########################## AUC Graph Together ######################

diclofenac_line + metformin_line + estradiol_line +
  # plot_annotation(tag_levels = 'A') + 
  plot_layout(axes = "collect")

ggsave("figures/Figure 2.tiff", dpi = 300, height = 8, width = 14)




