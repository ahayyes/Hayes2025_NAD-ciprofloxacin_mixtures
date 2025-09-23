################################ QPCR Analysis #############################

##################### Read in files and libraries#######################################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(lme4)
library(emmeans) 
library(DHARMa) 
library(rstatix)
library(ggpubr)
library(patchwork)

## colours 
# diclofenac "#79ad41"
# metformimn #34b6c6
# estradiol #F95639
# cip #B476FA

############################# Read in the Data and sort ################################

# list files
files <- list.files(path="data", full.names=TRUE)

#order = 
# ciprofloxacin, diclofenac, estradiol, metformin


# # create empty list to store results
output = list() #df

output2 = list() #df_melt

# output3 = list()


# run a for loop
for(i in 1:length(files)){
  cat(files[i])
  # read in file
  df <- read.csv(files[i], header=TRUE)
  
  #Give each specific sample tube a unique ID
  df$Rep <- paste(df$Sample, df$Rep) 
  
  #add in the prevalence data
  df <- df %>%
    mutate(d0_previntI1 = df$d0_intI1 / df$d0_16S) %>%
    mutate(d7_previntI1 = df$d7_intI1 / df$d7_16S)
  
  #Average the technical replicates for each dataset
  df_avg <- aggregate(df[c(7:8)], list(df$Rep), mean)
  
  #Melt for use
  df_melt <- melt(df_avg, id=c("Group.1"))
  
  #Change column name
  names(df_melt)[colnames(df_melt) == "Group.1"] <- "Sample"
  
  #Make Sample a factor so that it reads properly
  df_melt <- separate(df_melt, Sample, c("Sample", "Rep"), sep = " ")
  
  df_melt$Sample <- as.factor(df_melt$Sample)
  
  # need to add a column so you know which file and week this is
  df_melt <- mutate(df_melt, file = basename(tools::file_path_sans_ext(files[i])))
  
  # need to add a column so you know which file and week this is
  df_avg <- mutate(df_avg, file = basename(tools::file_path_sans_ext(files[i])))
  
  # then add this iteration to our empty list
  output[[i]] <- df
  
  # then add this iteration to our empty list
  output2[[i]] <- df_melt
  
  
  rm(df)
  rm(df_avg)
  rm(df_melt)
  
}



## Total Dfs in one!
# df_melt_total <- plyr::ldply(output2, data.frame)

df_ciprofloxacin <- plyr::ldply(output2[1], data.frame)

df_diclofenac <- plyr::ldply(output2[2], data.frame)

df_estradiol <- plyr::ldply(output2[3], data.frame)

df_metformin <- plyr::ldply(output2[4], data.frame)

# Set the levels of all dataframes
df_ciprofloxacin$Sample <- df_ciprofloxacin$Sample %>% factor(levels = c("0", "5", "10", "20", "40"))

df_diclofenac$Sample <- df_diclofenac$Sample %>% factor(levels = c("0","2.5", "5", "10", "20", "40"))

df_metformin$Sample <- df_metformin$Sample %>% factor(levels = c("0", "5", "10", "20", "40"))

df_estradiol$Sample <- df_estradiol$Sample %>% factor(levels = c("0", "5", "10", "20", "40"))


################################ Ciprofloxacin Analysis #######################
df_ciprofloxacin$concrep <- paste(df_ciprofloxacin$Sample, df_ciprofloxacin$Rep)

m1_cip <- lmer(value ~ variable * Sample + (1|concrep), data=df_ciprofloxacin)

summary(m1_cip)

plot(m1_cip)
qqnorm(resid(m1_cip))
qqline(resid(m1_cip))

#DHARMa stuff 
simulationOutput <- simulateResiduals(fittedModel = m1_cip, plot = T)

ggplot(data = df_ciprofloxacin, aes(x = resid(m1_cip))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

m2_cip <- lmer(log(value) ~ variable * Sample + (1|concrep), data=df_ciprofloxacin)

summary(m2_cip)

plot(m2_cip)
qqnorm(resid(m2_cip))
qqline(resid(m2_cip))

simulationOutput <- simulateResiduals(fittedModel = m2_cip, plot = T)

ggplot(data = df_ciprofloxacin, aes(x = resid(m2_cip))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

drop1(m2_cip, test="Chisq")
#interaction term not significant when data is better 


m3_cip <- lmer(log(value) ~ variable + Sample + (1|concrep), data=df_ciprofloxacin)

summary(m3_cip)

plot(m3_cip)
qqnorm(resid(m3_cip))
qqline(resid(m3_cip))

ggplot(data = df_ciprofloxacin, aes(x = resid(m3_cip))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

#Histogram looks somewhat whack but dharma says okay

simulationOutput <- simulateResiduals(fittedModel = m3_cip, plot = T)

Anova(m3_cip, test="Chisq", type="III") 


m4_cip <- lmer(log(value) ~ Sample + (1|concrep), data=df_ciprofloxacin)

plot(m4_cip)
qqnorm(resid(m4_cip))
qqline(resid(m4_cip))

ggplot(data = df_ciprofloxacin, aes(x = resid(m4_cip))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


simulationOutput <- simulateResiduals(fittedModel = m4_cip, plot = T)

summary(m4_cip)

drop1(m4_cip, test="Chisq")

Anova(m4_cip, test="Chisq", type="III") 
#Sample p=0.001054 

emmeans(m4_cip, list(pairwise ~ Sample), adjust ="fdr")
#40 is significant p=0.0115


#p<0.01 **
#p<0.05 *

stat.test1_cip <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj,
  "0", "40", "d7_previntI1", '*',
  "0", "5", "d7_previntI1", "NS",
  "0", "10", "d7_previntI1", "NS",
  "0", "20", "d7_previntI1", "NS")

#### Ciprofloxacin Boxplot

cip_plot <- ggplot(subset(df_ciprofloxacin, variable=="d7_previntI1"), aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "A) Ciprofloxacin Only")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#B476FA")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  guides(fill = "none") +
  stat_pvalue_manual(stat.test1_cip, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)


cip_plot


supp_cip_plot <- ggplot(df_ciprofloxacin, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "Ciprofloxacin Only")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#d2b4f3", "#B476FA"), name = "Day", labels = c("D0", "D7")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  stat_pvalue_manual(stat.test1_cip, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)



################################ Diclofenac Mixture Analysis #######################
df_diclofenac$concrep <- paste(df_diclofenac$Sample, df_diclofenac$Rep)

m1_diclo <- lmer(value ~ variable * Sample + (1|concrep), data=df_diclofenac)

summary(m1_diclo)

plot(m1_diclo)
qqnorm(resid(m1_diclo))
qqline(resid(m1_diclo))

ggplot(data = df_diclofenac, aes(x = resid(m1_diclo))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m1_diclo, plot = T)

m2_diclo <- lmer(log(value) ~ variable * Sample + (1|concrep), data=df_diclofenac)

summary(m2_diclo)

plot(m2_diclo)
qqnorm(resid(m2_diclo))
qqline(resid(m2_diclo))

simulationOutput <- simulateResiduals(fittedModel = m2_diclo, plot = T)

ggplot(data = df_diclofenac, aes(x = resid(m2_diclo))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

drop1(m2_diclo, test="Chisq")

Anova(m2_diclo, test="Chisq", type="III") 

emmeans(m2_diclo, list(pairwise ~ variable * Sample), adjust ="none")

# read in unadjusted and adjust only those comparisons we are interested in

pvalues_diclo <- read.csv("unadjusted_emmeans/Diclofenac_emmeansoutput_unadjusted.csv", header=TRUE)

pvalues_diclo$p.adj <- p.adjust(pvalues_diclo$p.value, method = "fdr")

# 10 = 0.0085 ***
# 20 = 0.028 *
# 40 = 0.0085 ***

#p<0.001 ***
#p<0.01 **
#p<0.05 *

stat.test1_diclo <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj,
  "0", "40", "d7_previntI1", '***',
  "0", "5", "d7_previntI1", "NS",
  "0", "10", "d7_previntI1", "***",
  "0", "20", "d7_previntI1", "*")

####  Boxplot
#remove 2.5 for plotting so they are the same 
df_diclofenac <- subset(df_diclofenac, Sample !="2.5")


diclo_plot <- ggplot(subset(df_diclofenac, variable=="d7_previntI1" & 
                              Sample !="2.5"), aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "B) Diclofenac Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#79ad41")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  guides(fill = "none") +
  stat_pvalue_manual(stat.test1_diclo, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)

diclo_plot

#supplementary figure

supp_diclo_plot <- ggplot(df_diclofenac, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "Diclofenac Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#b5ea7c", "#79ad41"), name="Day", labels=c("D0", "D7")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  stat_pvalue_manual(stat.test1_diclo, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)



################################ Metformin Mixture Analysis #######################

df_metformin$concrep <- paste(df_metformin$Sample, df_metformin$Rep)

m1_met <- lmer(value ~ variable * Sample + (1|concrep), data=df_metformin)

summary(m1_met)

plot(m1_met)
qqnorm(resid(m1_met))
qqline(resid(m1_met))

ggplot(data = df_metformin, aes(x = resid(m1_met))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m1_met, plot = T)


m2_met <- lmer(log(value) ~ variable * Sample + (1|concrep), data=df_metformin)

summary(m2_met)

plot(m2_met)
qqnorm(resid(m2_met))
qqline(resid(m2_met))

simulationOutput <- simulateResiduals(fittedModel = m2_diclo, plot = T)

ggplot(data = df_metformin, aes(x = resid(m2_met))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

drop1(m2_met, test="Chisq")

Anova(m2_met, test="Chisq", type="III") 

emmeans(m2_met, list(pairwise ~ variable * Sample), adjust ="none")

# read in unadjusted and adjust only those comparisons we are interested in

pvalues_met <- read.csv("unadjusted_emmeans/Metformin_emmeansoutput_unadjusted.csv", header=TRUE)

pvalues_met$p.adj <- p.adjust(pvalues_met$p.value, method = "fdr")

# 10 = 0.0016 ***
# 20 = 0.0238 *
# 40 = 0.063 NS

#p<0.001 ***
#p<0.01 **
#p<0.05 *

stat.test1_met <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj,
  "0", "40", "d7_previntI1", "NS",
  "0", "5", "d7_previntI1", "NS",
  "0", "10", "d7_previntI1", "***",
  "0", "20", "d7_previntI1", "*")


####  Boxplot

met_plot <- ggplot(subset(df_metformin, variable=="d7_previntI1"), aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "C) Metformin Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#34b6c6")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  guides(fill = "none") +
  stat_pvalue_manual(stat.test1_met, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)

met_plot

supp_met_plot <- ggplot(df_metformin, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "Metformin Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#bcedf3", "#34b6c6"), name="Day", labels=c("D0", "D7")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  stat_pvalue_manual(stat.test1_met, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)



######################### 17-beta-estradiol Mixture Analysis #######################

df_estradiol$concrep <- paste(df_estradiol$Sample, df_estradiol$Rep)


m1_estradiol <- lmer(value ~ variable * Sample + (1|concrep), data=df_estradiol)

summary(m1_estradiol)

plot(m1_estradiol)
qqnorm(resid(m1_estradiol))
qqline(resid(m1_estradiol))

ggplot(data = df_estradiol, aes(x = resid(m1_estradiol))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m1_estradiol, plot = T)


m2_estradiol <- lmer(log(value) ~ variable * Sample + (1|concrep), data=df_estradiol)

summary(m2_estradiol)

plot(m2_estradiol)
qqnorm(resid(m2_estradiol))
qqline(resid(m2_estradiol))

simulationOutput <- simulateResiduals(fittedModel = m2_estradiol, plot = T)

ggplot(data = df_estradiol, aes(x = resid(m2_estradiol))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

drop1(m2_estradiol, test="Chisq")

Anova(m2_estradiol, test="Chisq", type="III") 

emmeans(m2_estradiol, list(pairwise ~ variable * Sample), adjust ="fdr")

# 10 = <.0001 ***
# 20 = <.0001 ***
# 40 = <.0001 ***

#p<0.001 ***
#p<0.01 **
#p<0.05 *

stat.test1_estradiol <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj,
  "0", "40", "d7_previntI1", "***",
  "0", "5", "d7_previntI1", "NS",
  "0", "10", "d7_previntI1", "***",
  "0", "20", "d7_previntI1", "***")


####  Boxplot

estradiol_plot <- ggplot(subset(df_estradiol, variable=="d7_previntI1"), aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "D) 17-β-estradiol Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#F95639")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  guides(fill = "none") +
  stat_pvalue_manual(stat.test1_estradiol, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)

estradiol_plot

supp_estradiol_plot <- ggplot(df_estradiol, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",
       title = "17-β-estradiol Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#f59c8c", "#F95639"), name = "Day", labels= c("D0", "D7")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  stat_pvalue_manual(stat.test1_estradiol, y.position = 0.5, label = "p.adj", xmin="group2", xmax=NULL, size = 6)



############################## Patchworking the Plots ##############################

cip_plot + diclo_plot +  met_plot + estradiol_plot +
  plot_layout(axes = "collect")

# ggsave("figures/Figure 3.tiff", dpi = 300, height = 14, width = 14)


## Supplementary Plot with D0 data

supp_cip_plot + supp_diclo_plot +  supp_met_plot + supp_estradiol_plot +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(axes = "collect")

ggsave("figures/Supplementary Figure 1.tiff", dpi = 300, height = 14, width = 14)


################### Comparison Graphs between ciprofloxacin and mixtures ###################

## add extra column on all dfs to include file, and concrep
df_ciprofloxacin$file_rep <- paste(df_ciprofloxacin$file, df_ciprofloxacin$concrep)

df_diclofenac$file_rep <- paste(df_diclofenac$file, df_diclofenac$concrep)

#remove 2.5ug/L from diclofenac otherwise won't plot properly
df_diclofenac <- subset(df_diclofenac, Sample !="2.5")

df_metformin$file_rep <- paste(df_metformin$file, df_metformin$concrep)

df_estradiol$file_rep <- paste(df_estradiol$file, df_estradiol$concrep)

## make group df 

df_all <- rbind(df_ciprofloxacin, df_diclofenac, df_metformin, df_estradiol)

df_all$file <- df_all$file %>% factor(levels = c("ciprofloxacin" ,"diclofenac", "metformin", "estradiol"))

df_ciprofloxacin$Sample <- df_ciprofloxacin$Sample %>% factor(levels = c("0", "5", "10", "20", "40"))


total_comparison <- ggplot(subset(df_all, variable == "d7_previntI1"),aes(x=Sample, y=value, fill=file))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, alpha = 0.9)+
  geom_point(position=position_jitterdodge(), alpha=0.6, size = 4, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Ciprofloxacin concentration (µg/L)", y="intI1 prevalence",  title = "E) Comparison between Ciprofloxacin and Mixtures")+
  ylab(expression(paste(italic("intI1"), " prevalence")))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"),
                    name = "Treatment",
                    labels= c("Ciprofloxacin Alone", "Diclofenac Mixtures",
                              "Metformin Mixtures", "17-β-estradiol Mixtures")) +
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))

# ggsave("figures/Figure 3E.tiff", dpi = 300, height = 8, width = 10)


# incorporating day 0 and day 7 data 

m3_all <- lmer(log(value) ~ file * Sample + variable + (1|file_rep), data=df_all)


plot(m3_all)
qqnorm(resid(m3_all))
qqline(resid(m3_all))

simulationOutput <- simulateResiduals(fittedModel = m3_all, plot = T)

ggplot(data = df_all, aes(x = resid(m2_all))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

drop1(m3_all, test="Chisq")
#interaction is important

Anova(m3_all, test="Chisq", type="III") 
# file         20.259  3    0.00015 ***
#   Sample       28.361  4  1.054e-05 ***
#   variable     20.967  1  4.673e-06 ***
#   file:Sample  54.637 12  2.102e-07 ***

summary(m3_all)
# file_rep (Intercept) 0.0000   0.0000  

emmeans(m3_all, list(pairwise ~ file | Sample | variable), adjust ="fdr", regrid = "response")



################################# Figure 3 #############################

#make one grouped plot
fig3D <- cip_plot + diclo_plot +  met_plot + estradiol_plot +
  plot_layout(axes = "collect")


free(fig3D) / total_comparison +  plot_layout(axes = "collect") +
  plot_layout(heights = c(2, 1))

ggsave("figures/Figure 3.tiff", dpi = 300, height = 18, width = 14)

