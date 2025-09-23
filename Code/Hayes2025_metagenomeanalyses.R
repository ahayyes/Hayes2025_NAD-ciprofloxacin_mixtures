########################## Analysing AMR++ Output ##############################

#Notes on colours for graphs
#Diclofenac = green 
#Metformin = blue
#Beta estradiol = red
#Ciprofloxacin = purple


######################### Read in Libraries ###################################
library(tidyverse)
library(data.table)
library(ggplot2)
library(vegan)
# library(cowplot)
library(MetBrewer)
library(emmeans) 
library(DHARMa) 
library(ggpubr)
library(rstatix)
library(dunn.test)
library(broom)
library(multcomp)
library(reshape2)
library(patchwork)
library(apeglm) # might not work in all versions of R

#to install several packages for the first time - 

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")


# BiocManager::install("phyloseq", force=TRUE)
library(phyloseq)

# BiocManager::install("metagenomeSeq")
library(metagenomeSeq)

# BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

# BiocManager::install("DESeq2")
library(DESeq2)



###################### Read in Data Resistome ##################################

list_files_amr <- list.files(path="data/resistome", full.names=TRUE)

for (file in list_files_amr){
  #if the merged dataset doesn't exist, create it
  if (!exists("amr_data")){
    amr_data <- read.csv(file, header=T,sep=',', quote = "")
  }
  #if the merged dataset does exist, append to it 
  if (exists("amr_data")){
    temp_dataset <- read.csv(file, header=T, sep=',', quote = "")
    amr_data <- merge(amr_data, temp_dataset)
    rm(temp_dataset)
  }
}

#Make gene accession the row names
amr_data <- amr_data %>% remove_rownames %>% column_to_rownames(var="gene_accession") 

# convert this into an 'otu_table' format required for phyloseq
amr_data <- otu_table(amr_data, taxa_are_rows = T)

#read in gene annotations
annotations <- read.table('data/megares_annotations_v3.00.csv', header=T, row.names=1, sep=",", quote = "")

#convert this into a 'taxonomy table' format required for phyloseq
annotations <- tax_table(as.matrix(annotations))

# read in  metadata
metadata <- read.table('data/metadata/working_metadata.csv', header=T, sep=',', row.names = 1, quote = "")

# convert to 'sample_data' format for phyloseq
metadata <- sample_data(metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
amr.ps <- merge_phyloseq(amr_data, annotations, metadata)



################################ Taxonomy Work ########################################

#Read in Metadata
# metadata <- read.table('data/metadata/working_metadata.csv', header=T, sep=',',  row.names = 1, quote = "")

#Read in multiple files - taxonomy
list_files_tax <- list.files(path="data/taxonomy", full.names=TRUE)


for (file in list_files_tax){
  #if the merged dataset doesn't exist, create it
  if (!exists("microbiome_data")){
    microbiome_data <- read.csv(file, header=T,  sep=',', quote = "")
  }
  #if the merged dataset does exist, append to it 
  if (exists("microbiome_data")){
    temp_dataset <- read.csv(file, header=T, sep=',', quote = "")
    microbiome_data <- merge(microbiome_data, temp_dataset)
    rm(temp_dataset)
  }
}

#Make gene accession the row names
microbiome_data <- microbiome_data %>% remove_rownames %>% column_to_rownames(var="X") 

kraken_microbiome <- otu_table(microbiome_data, taxa_are_rows = TRUE)

# Repeat similar steps to what we did with the qiime2 taxonomy
kraken_taxonomy <- data.table(id=rownames(kraken_microbiome))
kraken_taxonomy[, c('domain',
                    'kingdom',
                    'phylum',
                    'class',
                    'order',
                    'family',
                    'genus',
                    'species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]

kraken_taxonomy[is.na(kraken_taxonomy)] <- "" # replacing the NAs with an empty string

# Conver to data.frame
kraken_taxonomy <- as.data.frame(kraken_taxonomy)

# Use the id variable to rename the row.names
row.names(kraken_taxonomy) <- kraken_taxonomy$id

# Remove the "id" and "kingdom" columns
kraken_taxonomy <- within(kraken_taxonomy, rm(id))
kraken_taxonomy <- within(kraken_taxonomy, rm(kingdom))

kraken_microbiome.ps <- merge_phyloseq(kraken_microbiome, tax_table(as.matrix(kraken_taxonomy)),sample_data(metadata))



################### Phyloseq Objects, with and without the double control ###############

#With the double control and inoculum  - AKA total data
kraken_microbiome.ps
amr.ps

#Without the double control and inoculum
kraken_microbiome_small.ps <- subset_samples(kraken_microbiome.ps, Treatment!="beta control")
kraken_microbiome_small.ps <- subset_samples(kraken_microbiome_small.ps, Treatment!="diclo control")
kraken_microbiome_small.ps <- subset_samples(kraken_microbiome_small.ps, Treatment!="met control")
kraken_microbiome_small.ps <- subset_samples(kraken_microbiome_small.ps, Treatment!="inoculum")

amr_small.ps <- subset_samples(amr.ps, Treatment!="beta control")
amr_small.ps <- subset_samples(amr_small.ps, Treatment!="diclo control")
amr_small.ps <- subset_samples(amr_small.ps, Treatment!="met control")
amr_small.ps <- subset_samples(amr_small.ps, Treatment!="inoculum")

# we want to have datasets where we are able to compare the mixtures themselves to the 
# ciprofloxacin treatments, hence making so many!


################################## Alpha Diversity ################################

# calculating richness and Shannon
microbiome_alpha <- estimate_richness(kraken_microbiome_small.ps, measures = c("Observed","Shannon"))
microbiome_alpha

amr_alpha <- estimate_richness(amr_small.ps, measures = c("Observed","Shannon"))
amr_alpha

# to make meaningful comparisons we need to include the metadata

# convert metadata into a data.frame
microbiome_metadata <- as(sample_data(kraken_microbiome_small.ps),"data.frame")

amr_metadata <- as(sample_data(amr_small.ps),"data.frame")

# combine the metadata and alpha diversity values
microbiome_alpha_meta <- cbind(microbiome_metadata, microbiome_alpha)

amr_alpha_meta <- cbind(amr_metadata, amr_alpha)


#### Testing effect of treatment on richness

#Resistome

m_rich_amr <- lm(Observed ~ Sample_type, data=subset(amr_alpha_meta, Treatment!="inoculum"))

summary(m_rich_amr)

emmeans(m_rich_amr, pairwise ~ Sample_type, adjust = "fdr")

# contrast                         estimate   SE df t.ratio p.value
# (beta-estradiol) - ciprofloxacin   -3.081 1.52 55  -2.031  0.0707
# (beta-estradiol) - diclofenac      -4.014 1.52 55  -2.646  0.0212
# (beta-estradiol) - metformin        0.986 1.52 55   0.650  0.5339
# ciprofloxacin - diclofenac         -0.933 1.49 55  -0.626  0.5339
# ciprofloxacin - metformin           4.067 1.49 55   2.728  0.0212
# diclofenac - metformin              5.000 1.49 55   3.354  0.0087

#p<0.05 = *
#p<0.01 = **
#p<0.001 = ***


simulationOutput <- simulateResiduals(fittedModel = m_rich_amr, plot = T)
#The model fits well within the parameters

#Mixtures has some effect on richness of the resistome in the treatments

stat.test_amr_richness<-  tibble::tribble(
  ~group1, ~group2, ~p.adj, ~y, ~Treatment,
  "beta-estradiol", "diclofenac", "*", "Observed", "beta-estradiol",
  "ciprofloxacin", "metformin", "*", "Observed", "ciprofloxacin",
  "diclofenac", "metformin", "**", "Observed", "diclofenac")

stat.test_amr_richness<-  tibble::tribble(
  ~group1, ~group2, ~p.adj, ~y, ~Treatment,
  "ciprofloxacin", "metformin", "*", "Observed", "ciprofloxacin")

#Microbiome
m_richmicrobiome <- lm(Observed ~ Sample_type + Cip_concentration, data=subset(microbiome_alpha_meta, Treatment!="inoculum"))

summary(m_richmicrobiome)

drop1(m_richmicrobiome, test="Chisq")

emmeans(m_richmicrobiome, pairwise ~ Sample_type, adjust = "fdr")

simulationOutput <- simulateResiduals(fittedModel = m_richmicrobiome, plot = T)


############################## Plotting richness ####################################

amr_alpha_meta$Sample_type <- factor(amr_alpha_meta$Sample_type, levels = c("ciprofloxacin", "diclofenac", "metformin", "beta-estradiol"))
amr_alpha_meta$Cip_concentration <- as.numeric(amr_alpha_meta$Cip_concentration)

ggplot(amr_alpha_meta) +
  geom_boxplot(aes(x= Sample_type, y= Observed), alpha = 0.9, outlier.shape = NA, linewidth = 1) +
  geom_point(aes(x=Sample_type, y=Observed, colour=Cip_concentration), size = 4, alpha = 0.8, position=position_dodge2(width = 0.5)) +
  scale_color_gradientn(colours = met.brewer("Hokusai3"))+
  theme_classic()+
  stat_pvalue_manual(stat.test_amr_richness, 
                     label = "p.adj", tip.length = 0.04, y.position = 375, size = 8,
                     bracket.size = 0.9, step.increase = 0.2)+
  theme(text=element_text(size=14))+
  labs(y= "Observed Genes", x="Treatment", color="Ciprofloxacin \nConcentration")+
  scale_x_discrete(labels = c("Ciprofloxacin \nAlone", "Diclofenac \nMixture", "Metformin \nMixture", "17-β-estradiol \nMixture"))

ggsave("figures/Figure 4.tiff", dpi=300, width = 8, height = 6)


########## Testing evenness between treatments 

#Microbiome
m_even_microbiome <- lm(Shannon ~ Sample_type * Cip_concentration, data=subset(microbiome_alpha_meta))

summary(m_even_microbiome)

drop1(m_even_microbiome, test="Chisq")
#interaction is important in the model p=0.0023

simulationOutput <- simulateResiduals(fittedModel = m_even_microbiome, plot = T)

anova(m_even_microbiome)
#p=0.0054
emmeans(m_even_microbiome, pairwise ~ Sample_type | Cip_concentration, , adjust = "fdr")

#No difference

#Resistome
m_even_amr <- lm(Shannon ~ Sample_type, data=subset(amr_alpha_meta))
summary(m_even_amr)

simulationOutput <- simulateResiduals(fittedModel = m_even_amr, plot = T)



############### Normalisation and Cumulative Sum Scaling (CSS) ###################

# First, convert the phyloseq object to metagenomeSeq
microbiome.metaseq <- phyloseq_to_metagenomeSeq(kraken_microbiome.ps) 
amr.metaseq <- phyloseq_to_metagenomeSeq(amr.ps)

# Perform normalization
microbiome_css.metaseq <- cumNorm(microbiome.metaseq)
amr_css.metaseq <- cumNorm(amr.metaseq)

# Here,  need to use MRcounts() and re-make the phyloseq object with the normalized counts
CSS_microbiome_counts <- MRcounts(microbiome.metaseq, norm = TRUE)
CSS_amr_counts <- MRcounts(amr.metaseq, norm = TRUE)

# Use the new counts and merge with components from the original phyloseq object.
microbiome_css.ps <- merge_phyloseq(otu_table(CSS_microbiome_counts, taxa_are_rows = TRUE),sample_data(kraken_microbiome.ps),
                                    tax_table(kraken_microbiome.ps))

amr_css.ps <- merge_phyloseq(otu_table(CSS_amr_counts, taxa_are_rows = TRUE),sample_data(amr.ps),tax_table(amr.ps))

######################### Phyloseq Objects #######################################

# Here are the phyloseq objects we have:
kraken_microbiome.ps # not normalized; 16S ASV counts, tree, rep-seqs, taxonomy, tree, metadata
amr.ps # not normalized; ARG counts, taxonomy, metadata

kraken_microbiome_small.ps #Small non normalised version
amr_small.ps # small non normalised version

microbiome_css.ps # CSS normalized; 16S ASV counts, tree, rep-seqs, taxonomy, tree, metadata
amr_css.ps # CCS normalized; ARG counts, taxonomy, metadata


## No Inoc ##

kraken_microbiome_noinoc.ps <- subset_samples(kraken_microbiome.ps, Treatment!="inoculum") # not normalized; 16S ASV counts, tree, rep-seqs, taxonomy, tree, metadata
amr_noinoc.ps <- subset_samples(amr.ps, Treatment!="inoculum")   # not normalized; ARG counts, taxonomy, metadata

microbiome_css_noinoc.ps <- subset_samples(microbiome_css.ps, Treatment!="inoculum") # CSS normalized; 16S ASV counts, tree, rep-seqs, taxonomy, tree, metadata
amr_css_noinoc.ps <- subset_samples(amr_css.ps, Treatment!="inoculum") # CCS normalized; ARG counts, taxonomy, metadata


#Make data sets of these 

amr_css_df <- psmelt(amr_css_noinoc.ps)
microbiome_css_df <- psmelt(microbiome_css_noinoc.ps)


################################# Beta Diversity #############################

## Ordination ##

microbiome_bray.dist <- vegdist(t(otu_table(microbiome_css_noinoc.ps)), method = "bray")
amr_bray.dist <- vegdist(t(otu_table(amr_css_noinoc.ps)), method = "bray")


#ordinate the distance matrices
microbiome_bray.ord <- ordinate(microbiome_css_noinoc.ps, method = "NMDS", distance = microbiome_bray.dist)
amr_bray.ord <- ordinate(amr_css_noinoc.ps, method = "NMDS", distance = amr_bray.dist)


#plot

#### 16S
plot_ordination(microbiome_css_noinoc.ps, microbiome_bray.ord, type = "samples",
                color = "Sample_type") +
  theme_bw() +
  stat_ellipse(lty = 2) +
  geom_point(size = 4, shape = 18) +
  scale_color_manual(values = c("#F95639", "#B476FA", "#79ad41", "#34b6c6"), name="Treatment", 
                     labels = c("17-β-estradiol mixture",
                                "Ciprofloxacin alone", 
                                "Diclofenac mixture", "Metformin mixture"))+
  theme(text=element_text(size=14))

ggsave("figures/Supplementary Figure 2.tiff", dpi = 300, height = 6, width = 8)

# resistome
plot_ordination(amr_css_noinoc.ps, amr_bray.ord, type = "samples", color = "Sample_type") +
  theme_bw() +
  stat_ellipse(lty = 2) +
  geom_point(size = 4, shape = 18) +
  scale_color_manual(values = c("#F95639", "#B476FA", "#79ad41", "#34b6c6"), name="Treatment", 
                     labels = c("17-β-estradiol mixture",
                                "Ciprofloxacin alone", 
                                "Diclofenac mixture", "Metformin mixture"))+  theme(text=element_text(size=14))+
  labs(colour = "Treatment")

ggsave("figures/Supplementary Figure 3.tiff", dpi = 300, height = 6, width = 8)



### Analysis of similarities (ANOSIM) ###

# The anosim function performs a non-parametric test of the significance of the 
# sample-grouping variable you provide against a permutation-based null distribution.

## Resistome

# First, make the R object with the sample grouping variable.
# ANOSIM only allows for the comparison of samples from 2 different groups.
group_variable = get_variable(amr_css_noinoc.ps,"Sample_type")

anosim_amr_by_group = anosim(phyloseq::distance(amr_css_noinoc.ps, "bray"), group_variable)
anosim_amr_by_group

plot(anosim_amr_by_group)



#Taxonomy
group_variable = get_variable(microbiome_css_noinoc.ps,"Sample_type")

anosim_phylum_by_group = anosim(phyloseq::distance(microbiome_css_noinoc.ps, "bray"), group_variable)
anosim_phylum_by_group

plot(anosim_phylum_by_group)


########################## Filter out low abundance species ######################

#make relative abundance dataset

# Convert OTU abundances to relative abundances (I prefer making it a % opposed to proportion, but up to you)
microbiome_ra.ps <- transform_sample_counts(microbiome_css_noinoc.ps, function(x) {x/sum(x)}*100)

#have full dataset too
microbiome_ra.melt <- psmelt(microbiome_ra.ps)

microbiome_high_abun_genus<- microbiome_ra.melt %>%
  group_by(Treatment, genus) %>%
  mutate(mean_abun = mean(Abundance))

microbiome_high_abun_genus$genus <- as.character(microbiome_high_abun_genus$genus)

microbiome_high_abun_genus$mean_abun<- as.numeric(microbiome_high_abun_genus$mean_abun)

microbiome_high_abun_genus$genus[microbiome_high_abun_genus$mean_abun < 0.1] <- "Low abundance genera"

unique(microbiome_high_abun_genus$genus)


microbiome_high_abun_genus <- subset(microbiome_high_abun_genus, Treatment!="diclo control")
microbiome_high_abun_genus <- subset(microbiome_high_abun_genus, Treatment!="met control")
microbiome_high_abun_genus <- subset(microbiome_high_abun_genus, Treatment!="beta control")

microbiome_high_abun_genus$Treatment <- factor(microbiome_high_abun_genus$Treatment, 
                                               levels=c("cip 0","cip 5", "cip 10", "cip 20", "cip 40",
                                                        "diclo 0", "diclo 5", "diclo 10", "diclo 20", "diclo 40",
                                                        "met 0", "met 5","met 10", "met 20", "met 40",
                                                        "beta 0", "beta 5", "beta 10", "beta 20", "beta 40"))



genera_highabun <- ggplot(subset(microbiome_high_abun_genus, genus!="" & domain =="Bacteria"),
                          aes(x=Treatment, y=mean_abun, fill=genus)) +
  geom_bar(position="fill", stat="identity")+ 
  theme_classic()  +
  theme(text=element_text(size=20)) +
  scale_fill_manual(values=met.brewer("Hiroshige", 17))+
  geom_bar(position="fill", stat="identity")+
  labs(y= "Relative Abundance", x= "Treatment", fill = "Genus")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(legend.text=element_text(face="italic"))


phylum_plot <- ggplot(subset(microbiome_high_abun_genus, phylum!="" & domain=="Bacteria"), 
                      aes(x=Treatment, y= Abundance, fill = phylum)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values=met.brewer("Hiroshige", 23))+
  theme_classic()  +
  theme(text=element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Relative Abundance", x= "Treatment", fill = "Phylum")

genera_highabun + phylum_plot +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(axes = "collect")


ggsave("figures/Supplementary Figure 4.tiff", height = 8, width = 20, dpi = 300)





amr_css_df$Treatment <- as.factor(amr_css_df$Treatment)
#Sort the levels
amr_css_df$Treatment <- factor(amr_css_df$Treatment, 
                               levels=c("cip 0","cip 5", "cip 10", "cip 20", "cip 40",
                                        "diclo control", "diclo 0", "diclo 5", "diclo 10", "diclo 20", "diclo 40",
                                        "met control", "met 0", "met 5","met 10", "met 20", "met 40",
                                        "beta control", "beta 0", "beta 5", "beta 10", "beta 20", "beta 40"))




############# use these datasets for stats #############################

amr_css_df
microbiome_css_df 

# these are normalised data not including the inoculum - we want to know the difference in the evolved treatments

######################### Testing to see if any genes differ between treatments ################ 

# run a Kruskal wallis test on each gene ###

# set up empty results dataframe
models_kruskal1 <- dplyr::select(amr_css_df, group) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))


# kruskal.test(Abundance ~ Treatment, data=subset(amr_ra.melt, mechanism=="Drug and biocide and metal RND efflux pumps"))


for(i in 1:nrow(models_kruskal1)){
  
  # grab amr class
  temp_amr_class <- models_kruskal1$group[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(amr_css_df, group==temp_amr_class) %>%
    dplyr::select(Treatment, Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Abundance ~ Treatment, data=temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal1$data[[i]] <- temp_data
  models_kruskal1$model[[i]] <- temp_model
  models_kruskal1$summary[[i]] <- temp_output
  
}

# look at one of the models
models_kruskal1$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal1, data) %>% dplyr::select(-model, -summary)

# 2. get the output of the model
model_output1 <- unnest(models_kruskal1, summary) %>% dplyr::select(-model, -data) 

# add in a column for p adjustments
# can look at this column to determine which are significant across treatments
model_output1 <- mutate(model_output1, padj = p.adjust(p.value, method = 'fdr'))


model_output1_sig <- subset(model_output1, padj<0.05)

# write.csv(model_output1_sig, "Genes that significantly alter with treatment.csv")


#Make vector of the genes of interest
all_genes_of_interest <- model_output1_sig$group

genes_of_interest_df <- subset(amr_css_df, group %in% all_genes_of_interest)


#Sort the levels
genes_of_interest_df$Treatment <- factor(genes_of_interest_df$Treatment, 
                                         levels=c("cip 0","cip 5", "cip 10", "cip 20", "cip 40",
                                                  "diclo control", "diclo 0", "diclo 5", "diclo 10", "diclo 20", "diclo 40",
                                                  "met control", "met 0", "met 5","met 10", "met 20", "met 40",
                                                  "beta control", "beta 0", "beta 5", "beta 10", "beta 20", "beta 40"))


#Separate mx and conc
genes_of_interest_df <- genes_of_interest_df %>% separate(Treatment, c("mix", "conc"), sep = " ")

genes_of_interest_df$mix <- as.factor(genes_of_interest_df$mix)
# genes_of_interest_df$conc <- as.factor(genes_of_interest_df$conc)

#set levels
genes_of_interest_df$mix <- factor(genes_of_interest_df$mix,
                                   levels=c("cip", "diclo", "met", "beta"))



genes_of_interest_df <- subset(genes_of_interest_df, conc!="control")
genes_of_interest_df$conc <- as.numeric(as.character(genes_of_interest_df$conc))

#set levels
genes_of_interest_df$mix <- factor(genes_of_interest_df$mix,
                                   levels=c("cip", "diclo", "met", "beta"))


genes_of_interest_summ <- genes_of_interest_df %>%
  group_by(group, mix, conc) %>%
  summarise(avg = mean(Abundance),
            sd = sd(Abundance))

# write.csv(genes_of_interest_summ, "summarised genes of interest.csv")


############################## CTX  #####################################

m_ctx <- lm(Abundance ~ conc * mix,
            data=subset(genes_of_interest_df, group=="CTX"))

summary(m_ctx)

m_ctx_unimodal <- lm(Abundance ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="CTX"))

summary(m_ctx_unimodal)

anova(m_ctx, m_ctx_unimodal)
# unimodal is better than linear

drop1(m_ctx_unimodal, test="Chisq")

m_ctx_unimodal2 <- update(m_ctx_unimodal, ~. -mix:I(conc^2))

anova(m_ctx_unimodal, m_ctx_unimodal2)

drop1(m_ctx_unimodal2, test="Chisq")

m_ctx_unimodal3 <- update(m_ctx_unimodal2, ~. -conc:mix)

anova(m_ctx_unimodal2, m_ctx_unimodal3)
# no difference

drop1(m_ctx_unimodal3, test="Chisq")
#mix is not significant

## mixture type does not alter ctx abundance. no further analysis

summary(m_ctx_unimodal3)
Anova(m_ctx_unimodal3, test="F")
# mix is not significant, no further analysis 

ggplot(data = subset(genes_of_interest_df, group=="CTX"), aes(x = resid(m_ctx_unimodal3))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

qqnorm(resid(m_ctx_unimodal3))
qqline(resid(m_ctx_unimodal3))

simulationOutput <- simulateResiduals(fittedModel = m_ctx_unimodal3, plot = T)


m_ctx_unimodal4 <- lm(Abundance ~ conc + I(conc^2), data = subset(genes_of_interest_df, 
                                                                  group == "CTX"))

summary(m_ctx_unimodal4)

drop1(m_ctx_unimodal4, test="Chisq")


ggplot(data = subset(genes_of_interest_df, group=="CTX"), aes(x = resid(m_ctx_unimodal4))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_ctx_unimodal4, plot = T)

# this data fits ok. mix is not sig, no further anlaysis 




#######################################  QNRB  ###################################

m_qnrB <- lm(Abundance ~ conc* mix, data=subset(genes_of_interest_df, group=="QNRB"))

summary(m_qnrB)

m_qnrB_unimodal <- lm(Abundance ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="QNRB"))

summary(m_qnrB_unimodal)

anova(m_qnrB, m_qnrB_unimodal)

# unimodal is not better fit, stick with non-unimodal

drop1(m_qnrB, test="Chisq")
#Interaction is significant for the model
#p=0.0318, the interaction explains some of the variance

#check residuals

# par(mfrow = c(1, 2)) # combine plots

# 1. Homogeneity of variances
plot(m_qnrB, which = 3)

# 2. Normality
plot(m_qnrB, which = 2)


ggplot(data = subset(genes_of_interest_df, group=="QNRB"), aes(x = resid(m_qnrB))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_qnrB, plot = T)

#Overall we will use this model as pretty okay fit and transforming it makes it worse

Anova(m_qnrB, test = "F")
#mix is significant p=0.047
# conc is sig p<0.001
#interaction is p=0.053

#make predictions from model to plot more accurately

preds_qnB <- subset(genes_of_interest_df, group=="QNRB") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_qnB <- cbind(preds_qnB, predict(m_qnrB, newdata = preds_qnB, interval = 'confidence'))


#Plot this

plot_qnrB <- ggplot() +
  geom_ribbon(data = preds_qnB, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_line(data = preds_qnB, aes(x=conc, y=fit, colour = mix), linewidth=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="QNRB"), 
             mapping = aes(x=conc, y= Abundance, colour = mix), 
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  theme_classic()  +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20))+
  facet_wrap(~group) +
  theme(legend.position="none")

plot_qnrB

################################### APH3-DPRIME #################################

#Stats
m_aph3 <- lm(log(Abundance+1) ~ conc * mix,
             data=subset(genes_of_interest_df, group=="APH3-DPRIME"))

m_aph3 <- lm(Abundance ~ conc * mix,
             data=subset(genes_of_interest_df, group=="APH3-DPRIME"))

summary(m_aph3)


m_aph3_unimodal <- lm(Abundance ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="APH3-DPRIME"))

summary(m_aph3_unimodal)

anova(m_aph3, m_aph3_unimodal)
# unimodal is no better than linear, stay linear

drop1(m_aph3, test="Chisq")
#Interaction is  significant for the model

#check residuals

# par(mfrow = c(1, 2)) # combine plots

# 1. Homogeneity of variances
plot(m_aph3, which = 3)

# 2. Normality
plot(m_aph3, which = 2)

ggplot(data=subset(genes_of_interest_df, group=="APH3-DPRIME"), aes(x = resid(m_aph3))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_aph3, plot = T)

#do anova
Anova(m_aph3, test = "F")
# conc p<0.001
# mix p<0.001
# interaction p=0.49


#make predictions from model to plot more accurately

preds_aph3 <- subset(genes_of_interest_df, group=="APH3-DPRIME") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_aph3 <- cbind(preds_aph3, predict(m_aph3, newdata = preds_aph3, interval = 'confidence'))


#Plot this


# very small effect, can only see when 'zoomed in' on y axis

plot_aph3 <- ggplot() +
  geom_ribbon(data = preds_aph3, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_line(data = preds_aph3, aes(x=conc, y=fit, colour = mix), size=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="APH3-DPRIME"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 4), size = 4, alpha = 0.4)+
  theme_classic()  +
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Log Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~group) +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20))+
  theme(legend.position="none") +
  scale_y_log10()



########################################## APH6 ##########################################
m_aph6 <- lm(log(Abundance+1) ~ mix * conc, data=subset(genes_of_interest_df, group=="APH6"))

m_aph6_unimodal <- lm(log(Abundance+1) ~ mix * conc + mix * I(conc^2),
                      data=subset(genes_of_interest_df, group=="APH6"))

summary(m_aph6)

anova(m_aph6, m_aph6_unimodal)

# unimodal is better than not

drop1(m_aph6_unimodal, test="Chisq")

# all important

plot(m_aph6_unimodal, which = 3)

# 2. Normality
plot(m_aph6_unimodal, which = 2)

ggplot(data = subset(genes_of_interest_df, group=="APH6"), aes(x = resid(m_aph6_unimodal))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_aph6_unimodal, plot = T)

#do anova
Anova(m_aph6_unimodal, test = "F")

#unimodal term of conc : mix p=0.0060
# mix:con p=0.05
#mix is significant p<0.0001


#make predictions from model to plot more accurately

preds_aph6 <- subset(genes_of_interest_df, group=="APH6") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_aph6 <- cbind(preds_aph6, predict(m_aph6_unimodal, newdata = preds_aph6, interval = 'confidence'))


#Plot

plot_aph6 <- ggplot() +
  geom_ribbon(data = preds_aph6, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_line(data = preds_aph6, aes(x=conc, y=fit, colour = mix), size=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="APH6"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  theme_classic()  +
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Log Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~group) +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20))+
  theme(legend.position="bottom")+
  scale_y_log10()


###################################### FECE ################################
m_fece <- lm(Abundance ~ mix * conc,
             data=subset(genes_of_interest_df, group=="FECE"))

summary(m_fece)

drop1(m_fece, test="Chisq")
#Interaction is  significant for the model

plot(m_fece, which = 3)

plot(m_fece, which = 2)

ggplot(data = subset(genes_of_interest_df, group=="FECE"), aes(x = resid(m_fece))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_fece, plot = T)


#do anova
Anova(m_fece, test = "F")
# mixture is significant in explaiing the data as both main effect p=0.017
# and as an interaction with cip concentration p<0.001


#make predictions to plot more accurately
preds_fece <- subset(genes_of_interest_df, group=="FECE") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_fece <- cbind(preds_fece, predict(m_fece, newdata = preds_fece, interval = 'confidence'))


#Plot this

plot_fece <- ggplot() +
  geom_ribbon(data = preds_fece, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_line(data = preds_fece, aes(x=conc, y=fit, colour = mix), size=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="FECE"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  theme_classic()  +
  theme(text=element_text(size=20))+
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~group) +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20))+
  theme(legend.position="none")

## very strong effect


######################################### SILC ###########################
m_silc <- lm(Abundance ~ mix * conc,
             data=subset(genes_of_interest_df, group=="SILC"))

summary(m_silc)

drop1(m_silc, test="Chisq")
#Interaction is not significant for the model


m_silc_unimodal <- lm(Abundance ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="SILC"))

summary(m_silc_unimodal)

anova(m_silc, m_silc_unimodal)
# unimodal is not better stay linear

drop1(m_silc, test="Chisq")
#mix is not significant

m_silc <- lm(Abundance ~ mix + conc,
             data=subset(genes_of_interest_df, group=="SILC"))

drop1(m_silc, test="Chisq")
# mix not sig p=0.29
# conc is 
# no further analysis 

plot(m_silc, which = 3)

plot(m_silc, which = 2)


ggplot(data = subset(genes_of_interest_df, group=="SILC"), aes(x = resid(m_silc))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_silc, plot = T)



################################## SULII ######################################
m_sulII <- lm(Abundance ~ mix * conc,
              data=subset(genes_of_interest_df, group=="SULII"))

summary(m_sulII)


m_sulII_unimodal <- lm(Abundance ~ mix * conc + mix * I(conc^2),
                       data=subset(genes_of_interest_df, group=="SULII"))

summary(m_sulII_unimodal)

anova(m_sulII, m_sulII_unimodal)
# unimodal is better than not

drop1(m_sulII_unimodal, test="Chisq")
#Interactions are not significant for the model

m_sulII_unimodal2 <- update(m_sulII_unimodal, ~. -mix:conc)

drop1(m_sulII_unimodal2, test="Chisq")

m_sulII_unimodal3 <- update(m_sulII_unimodal2, ~. -mix:I(conc^2))

drop1(m_sulII_unimodal3, test="Chisq")


plot(m_sulII_unimodal3, which = 3)

# 2. Normality
plot(m_sulII_unimodal3, which = 2)

ggplot(data = subset(genes_of_interest_df, group=="SULII"), aes(x = resid(m_sulII_unimodal3))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_sulII_unimodal3, plot = T)

Anova(m_sulII_unimodal3, test = "F")

# mix is not significant p=0.67, no further analysis 

############################## TETA #########################################


m_teta <- lm(Abundance ~ mix * conc,
             data=subset(genes_of_interest_df, group=="TETA"))


m_teta_unimodal <- lm(Abundance ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="TETA"))

summary(m_teta_unimodal)

anova(m_teta, m_teta_unimodal)
# unimodal is better than not

drop1(m_teta_unimodal, test="Chisq")

# all sig 

simulationOutput <- simulateResiduals(fittedModel = m_teta_unimodal, plot = T)
# not amaze

m_teta_unimodal <- lm(log(Abundance+1) ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="TETA"))

drop1(m_teta_unimodal, test="Chisq")

simulationOutput <- simulateResiduals(fittedModel = m_teta_unimodal, plot = T)

# best fit we get

plot(m_teta_unimodal, which = 3)

# 2. Normality
plot(m_teta_unimodal, which = 2)

ggplot(data = subset(genes_of_interest_df, group=="TETA"), aes(x = resid(m_teta_unimodal))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

Anova(m_teta_unimodal, test = "F")
# mix is sig as a main effect p=0.0009
# and as a interaction iwith linear conc p=0.0063
# and as an interaction with unimodal conc p=0.00026

#make predictions to plot more accurately
preds_tetA <- subset(genes_of_interest_df, group=="TETA") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_tetA <- cbind(preds_tetA, predict(m_teta_unimodal, newdata = preds_tetA, interval = 'confidence'))



#Plot this

plot_tetA <- ggplot() +
  geom_ribbon(data = preds_tetA, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_line(data = preds_tetA, aes(x=conc, y=fit, colour = mix), size=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="TETA"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  theme_classic()  +
  theme(text=element_text(size=20))+
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Log Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~group) +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20))+
  theme(legend.position="none")+
  scale_y_log10()


plot_tetA

############################## TETQ #########################################

m_tetq <- lm(Abundance ~ mix * conc,
             data=subset(genes_of_interest_df, group=="TETQ"))

summary(m_tetq)

m_tetq_unimodal <- lm(Abundance ~ mix * conc + mix * I(conc^2),
                      data=subset(genes_of_interest_df, group=="TETQ"))

anova(m_tetq, m_tetq_unimodal)

# linear no worse than unimodal, stay linear
drop1(m_tetq, test="Chisq")

m_tetq <- lm(Abundance ~ mix + conc,
             data=subset(genes_of_interest_df, group=="TETQ"))

drop1(m_tetq, test="Chisq")

m_tetq <- lm(Abundance ~ mix,
             data=subset(genes_of_interest_df, group=="TETQ"))

drop1(m_tetq, test="Chisq")

# mix not conc important 

plot(m_tetq, which = 3)

# 2. Normality
plot(m_tetq, which = 2)

ggplot(data = subset(genes_of_interest_df, group=="TETQ"), aes(x = resid(m_tetq))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_tetq, plot = T)

m_tetq <- lm(log(Abundance+1) ~ mix,
             data=subset(genes_of_interest_df, group=="TETQ"))

ggplot(data = subset(genes_of_interest_df, group=="TETQ"), aes(x = resid(m_tetq))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_tetq, plot = T)


Anova(m_tetq, test="F")
#mix = 0.0023

#make predictions to plot more accurately
preds_tetQ <- subset(genes_of_interest_df, group=="TETQ") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_tetQ <- cbind(preds_tetQ, predict(m_tetq, newdata = preds_tetQ, interval = 'confidence'))



#Plot this

plot_tetQ <- ggplot() +
  geom_ribbon(data = preds_tetQ, aes(ymin=lwr, ymax=upr, x=conc, fill = mix),alpha=0.1) +
  geom_line(data = preds_tetQ, aes(x=conc, y=fit, colour = mix), size=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="TETQ"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  theme_classic()  +
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Log Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~group) +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20))+
  theme(legend.position="none")+
  scale_y_log10()

########################### TOLC ###############################################

m_tolc <- lm(Abundance ~ mix * conc,
             data=subset(genes_of_interest_df, group=="TOLC"))

m_tolc_unimodal <- lm(Abundance ~  conc * mix + I(conc^2) * mix, data=subset(genes_of_interest_df, group=="TOLC"))

summary(m_tolc_unimodal)

anova(m_tolc, m_tolc_unimodal)
# unimodal not better

summary(m_tolc)

drop1(m_tolc, test = "Chisq")

plot(m_tolc, which = 3)

plot(m_tolc, which = 2)

ggplot(data = subset(genes_of_interest_df, group=="TOLC"), aes(x = resid(m_tolc))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_tolc, plot = T)

qqnorm(m_tolc$resid)
qqline(m_tolc$resid)

Anova(m_tolc, test = "F")
# mixx main effect = p=0.0091
# mix interaction p<0.001
# conc main effect p=0.0012


#make predictions to plot more accurately
preds_tolC <- subset(genes_of_interest_df, group=="TOLC") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_tolC <- cbind(preds_tolC, predict(m_tolc, newdata = preds_tolC, interval = 'confidence'))



#Plot this

plot_tolC <-  ggplot() +
  geom_ribbon(data = preds_tolC, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_line(data = preds_tolC, aes(x=conc, y=fit, colour = mix), size=1.22)+
  geom_point(data = subset(genes_of_interest_df, group=="TOLC"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  theme_classic()  +
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~group) +
  theme(strip.text = element_text(face = "italic"),
        text=element_text(size=20)) +
  theme(legend.position = "none")


######################### Plot Genes together ###################################

plot_fece + plot_qnrB + plot_tolC + plot_spacer() +  plot_aph3 + plot_aph6 + plot_tetA + plot_tetQ + 
  plot_layout(axis_titles = "collect", ncol = 4) 

ggsave("figures/Figure 5.tiff", height = 8, width = 14)


################ Testing to see if any species differs between treatments ###############


# run a Kruskal wallis test on each gene ###

# set up empty results dataframe
models_kruskal2 <- dplyr::select(microbiome_css_df, order) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))


for(i in 1:nrow(models_kruskal2)){
  
  # grab amr class
  temp_species_class <- models_kruskal2$order[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(microbiome_css_df, order==temp_species_class) %>%
    dplyr::select(Treatment, Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Abundance ~ Treatment, data=temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal2$data[[i]] <- temp_data
  models_kruskal2$model[[i]] <- temp_model
  models_kruskal2$summary[[i]] <- temp_output
  
}

# models_kruskal2

# look at one of the models
models_kruskal2$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal2, data) %>% dplyr::select(-model, -summary)

# 2. get the output of the model
model_output2 <- unnest(models_kruskal2, summary) %>% dplyr::select(-model, -data) 

# add in a column for p adjustments
# can look at this column to determine which are significant across treatments
model_output2 <- mutate(model_output2, padj = p.adjust(p.value, method = 'fdr'))


model_output2_sig <- subset(model_output2, padj<0.05)

# write.csv(model_output2_sig, "all order, only significant genus.csv")

#Make vector of the genes of interest
all_order_of_interest <- model_output2_sig$order

order_of_interest_df <- subset(microbiome_css_df, order %in% all_order_of_interest)


#Sort the levels
order_of_interest_df$Treatment <- factor(order_of_interest_df$Treatment, 
                                         levels=c("cip 0","cip 5", "cip 10", "cip 20", "cip 40",
                                                  "diclo control", "diclo 0", "diclo 5", "diclo 10", "diclo 20", "diclo 40",
                                                  "met control", "met 0", "met 5","met 10", "met 20", "met 40",
                                                  "beta control", "beta 0", "beta 5", "beta 10", "beta 20", "beta 40"))


#Separate mx and conc
order_of_interest_df <- order_of_interest_df %>% separate(Treatment, c("mix", "conc"), sep = " ")

order_of_interest_df$mix <- as.factor(order_of_interest_df$mix)



order_of_interest_df <- subset(order_of_interest_df, conc!="control")

order_of_interest_df$conc <- as.numeric(order_of_interest_df$conc)
#set levels
order_of_interest_df$mix <- factor(order_of_interest_df$mix,
                                   levels=c("cip", "diclo", "met", "beta"))


order_of_interest_summ <- order_of_interest_df %>%
  group_by(genus, mix, conc) %>%
  summarise(avg = mean(Abundance),
            sd = sd(Abundance))

# write.csv(order_of_interest_summ, "summarised order of interest.csv")


############## Analysis on Genera of interest #######################


models_kruskal3 <- dplyr::select(microbiome_css_df, genus) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))


for(i in 1:nrow(models_kruskal3)){
  
  # grab amr class
  temp_species_class <- models_kruskal3$genus[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(microbiome_css_df, genus==temp_species_class) %>%
    dplyr::select(Treatment, Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Abundance ~ Treatment, data=temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal3$data[[i]] <- temp_data
  models_kruskal3$model[[i]] <- temp_model
  models_kruskal3$summary[[i]] <- temp_output
  
}


# look at one of the models
models_kruskal3$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal3, data) %>% dplyr::select(-model, -summary)

# 2. get the output of the model
model_output3 <- unnest(models_kruskal3, summary) %>% dplyr::select(-model, -data) 

# add in a column for p adjustments
# can look at this column to determine which are significant across treatments
model_output3 <- mutate(model_output3, padj = p.adjust(p.value, method = 'fdr'))


model_output3_sig <- subset(model_output3, padj<0.05)

#Make vector of the genes of interest
all_genera_of_interest <- model_output3_sig$genus

genera_of_interest_df <- subset(microbiome_css_df, genus %in% all_genera_of_interest)



#Sort the levels
genera_of_interest_df$Treatment <- factor(genera_of_interest_df$Treatment, 
                                          levels=c("cip 0","cip 5", "cip 10", "cip 20", "cip 40",
                                                   "diclo control", "diclo 0", "diclo 5", "diclo 10", "diclo 20", "diclo 40",
                                                   "met control", "met 0", "met 5","met 10", "met 20", "met 40",
                                                   "beta control", "beta 0", "beta 5", "beta 10", "beta 20", "beta 40"))


#Separate mx and conc
genera_of_interest_df <- genera_of_interest_df %>% separate(Treatment, c("mix", "conc"), sep = " ")

genera_of_interest_df$mix <- as.factor(genera_of_interest_df$mix)


genera_of_interest_df <- subset(genera_of_interest_df, conc!="control")

genera_of_interest_df$conc <- as.numeric(genera_of_interest_df$conc)
#set levels
genera_of_interest_df$mix <- factor(genera_of_interest_df$mix,
                                    levels=c("cip", "diclo", "met", "beta"))

### of these ##caudovirales and thiotrichales and nitrosphaerales seem to change 
## same order as of interest, keep to order level analysis


############################ Nitrososphaerales ##########################

m_nitro <- lm(Abundance ~ mix * conc, data=subset(order_of_interest_df, order=="Nitrososphaerales"))

summary(m_nitro)

drop1(m_nitro, test="Chisq")
#ineraction not important

m_nitro <- lm(Abundance ~ mix + conc, data=subset(order_of_interest_df, order=="Nitrososphaerales"))

summary(m_nitro)

drop1(m_nitro, test="Chisq")

#check residuals

# par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(m_nitro, which = 3)

# 2. Normality
plot(m_nitro, which = 2)

ggplot(data = m_nitro, aes(x = resid(m_nitro))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_nitro, plot = T)


#do anova
nitro_aov <- aov(Abundance ~ mix + conc, data=subset(order_of_interest_df, order=="Nitrososphaerales"))

summary(nitro_aov)


#make predictions from model to plot more accurately

preds_nitro <- subset(order_of_interest_df, order =="Nitrososphaerales") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_nitro <- cbind(preds_nitro, predict(m_nitro, newdata = preds_nitro, interval = 'confidence'))


#Plot this


plot_nitro <- ggplot() +
  geom_ribbon(data = preds_nitro, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix), alpha=0.1) +
  geom_point(data = subset(order_of_interest_df, order=="Nitrososphaerales"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+
  geom_line(data = preds_nitro, aes(x=conc, y=fit, colour = mix), linewidth=1.22)+
  theme_classic()  +
  theme(text=element_text(size=20))+
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~order)+
  theme(legend.position="bottom")


########################### Caulobacterales ##################

m_caulo <- lm(Abundance ~ mix * conc, data=subset(order_of_interest_df, order=="Caulobacterales"))

summary(m_caulo)

drop1(m_caulo, test="Chisq")

#check residuals
# 1. Homogeneity of variances
plot(m_caulo, which = 3)

# 2. Normality
plot(m_caulo, which = 2)

ggplot(data = m_caulo, aes(x = resid(m_caulo))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m_caulo, plot = T)


#do anova
caulo_aov <- aov(Abundance ~ mix * conc, data=subset(order_of_interest_df, order=="Caulobacterales"))

summary(caulo_aov)


#make predictions from model to plot more accurately

preds_caulo <- subset(order_of_interest_df, order =="Caulobacterales") %>%
  group_by(conc, mix) %>%
  do(data.frame(expand.grid(conc = seq(0,40, length.out = 100),
                            mix = c("cip", "diclo", "met", "beta")))) %>%
  ungroup()


preds_caulo <- cbind(preds_caulo, predict(m_caulo, newdata = preds_caulo, interval = 'confidence'))



plot_caulo <- ggplot() +
  geom_ribbon(data = preds_caulo, aes(ymin=lwr, ymax=upr, x=conc,  fill = mix),alpha=0.1) +
  geom_point(data = subset(order_of_interest_df, order=="Caulobacterales"), 
             mapping = aes(x=conc, y= Abundance, colour = mix),
             position = position_jitterdodge(dodge.width = 3), size = 4, alpha = 0.6)+  theme_classic()  +
  geom_line(data = preds_caulo, aes(x=conc, y=fit, colour = mix), size=1.22)+
  theme(text=element_text(size=20))+
  scale_colour_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                      labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#B476FA", "#79ad41", "#34b6c6", "#F95639"), name="Treatment",
                    labels=c('Ciprofloxacin', "Diclofenac", "Metformin", "17-β-estradiol"))+
  labs(y= "Normalised Abundance", x= "Ciprofloxacin Concentration")+
  guides(size = "none")+
  facet_wrap(~order)+
  theme(legend.position="bottom")


## Plot the two together
(plot_caulo + plot_nitro) + plot_layout(axis_titles = "collect") +
  plot_layout(guides = "collect")  & theme(legend.position = 'bottom') 


ggsave("figures/Figure 6.tiff", dpi= 300, height = 6, width = 10)



############################ Fluoroquinolones ########################################

fluoro_df <- amr_css_df[which(amr_css_df$class=="Fluoroquinolones"),]


# run a Kruskal wallis test on each gene ###

# set up empty results dataframe
models_kruskal_quino <- dplyr::select(fluoro_df, group) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))


for(i in 1:nrow(models_kruskal_quino)){
  
  # grab amr class
  temp_amr_class <- models_kruskal_quino$group[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(fluoro_df, group==temp_amr_class) %>%
    dplyr::select(Treatment, Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Abundance ~ Treatment, data=temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal_quino$data[[i]] <- temp_data
  models_kruskal_quino$model[[i]] <- temp_model
  models_kruskal_quino$summary[[i]] <- temp_output
  
}

# models_kruskal_quino

# look at one of the models
models_kruskal_quino$model[[1]]

# 1. get the data frame again
unnest(models_kruskal_quino, data) %>% dplyr::select(-model, -summary)

# 2. get the output of the models_kruskal_quino
model_output_quino <- unnest(models_kruskal_quino, summary) %>% dplyr::select(-model, -data) 

# add in a column for p adjustments
# can look at this column to determine which are significant across treatments
model_output_quino <- mutate(model_output_quino, padj = p.adjust(p.value, method = 'fdr'))

model_output_quino_sig <- subset(model_output_quino, padj<0.05)


#None of the others are significant. we already test QNRB earlier


################################## Separate mix and conc ######################

# amr_df_split <- amr_css_df %>% separate(Treatment, c("mix", "conc"), sep = " ")
# 
# #Set levels 
# amr_df_split$mix <- as.factor(amr_df_split$mix)
# 
# amr_df_split$mix <- factor(amr_df_split$mix,
#                            levels=c("cip", "diclo", "met", "beta"))
# 
# 
# 
# amr_df_split <- subset(amr_df_split, conc!="control")
# 
# 
# amr_df_split$conc <- as.numeric(amr_df_split$conc)


############################ Comparing Means - Wilcoxon ############################

# 16S
pairwise.wilcox.test(microbiome_css_df$Abundance, microbiome_css_df$Sample_type, p.adjust.method = "fdr")

pairwise.wilcox.test(microbiome_css_df$Abundance, microbiome_css_df$Treatment, p.adjust.method = "fdr")

# AMR
pairwise.wilcox.test(amr_css_df$Abundance, amr_css_df$Sample_type, p.adjust.method = "fdr")


################################## DESeq 2 AMR #########################################

#This will tell us any log fold changes in AMR genes  

#We will use mixture as the variable of change 
#Do samples change with the mixture or not

#Use these datasets
kraken_microbiome_small.ps
amr_small.ps


#Make separate DFs
diclo_DA_df <- subset_samples(amr_small.ps, Sample_type=="diclofenac" | Sample_type =="ciprofloxacin")

#check
head(sample_data(diclo_DA_df)$Sample_type, 250)

#Create a DESeq2 object
desq_amr_diclo = phyloseq_to_deseq2(diclo_DA_df, ~ Sample_type)

#Test this 
desq_amr_diclo = DESeq(desq_amr_diclo, test="Wald", fitType="local")

#Look at the results
desq_amr_results_diclo = results(desq_amr_diclo, cooksCutoff = FALSE)

desq_amr_results_df_diclo <-  as.data.frame(desq_amr_results_diclo)
#None of these are significant



#### Metformin ###

#Make separate DFs
met_DA_df <- subset_samples(amr_small.ps, Sample_type=="metformin" | Sample_type =="ciprofloxacin")

#check
head(sample_data(met_DA_df)$Sample_type, 250)

#Create a DESeq2 object
desq_amr_met = phyloseq_to_deseq2(met_DA_df, ~ Sample_type)

#Test this 
desq_amr_met = DESeq(desq_amr_met, test="Wald", fitType="local")

#Look at the results
desq_amr_results_met = results(desq_amr_met, cooksCutoff = FALSE)

desq_amr_results_df_met <-  as.data.frame(desq_amr_results_met)
#None of these are significant


#### Beta-Estradiol ###

#Make separate DFs
beta_DA_df <- subset_samples(amr_small.ps, Sample_type=="beta-estradiol" | Sample_type =="ciprofloxacin")

#check
head(sample_data(beta_DA_df)$Sample_type, 250)

#Create a DESeq2 object
desq_amr_beta = phyloseq_to_deseq2(beta_DA_df, ~ Sample_type)

#Test this 
desq_amr_beta = DESeq(desq_amr_beta, test="Wald", fitType="local")

#Look at the results
desq_amr_results_beta = results(desq_amr_beta, cooksCutoff = FALSE)

desq_amr_results_df_beta <-  as.data.frame(desq_amr_results_beta)
#None of these are significant


# alpha = 0.05
# 
# sigtab_amr = desq_amr_results[which(desq_amr_results$padj < alpha), ]
# 
# sigtab_amr_df <- as.data.frame(sigtab_amr)
# 
# sigtab_amr = cbind(sigtab_amr_df, as(tax_table(amr_noinoc.ps)[rownames(sigtab_amr), ], "matrix"))
# 
# head(sigtab_amr)

plotDispEsts(desq_amr_beta)



# test_log <- lfcShrink(desq_amr, coef=2, type="apeglm") 

deseq_amr_beta <- data.frame(desq_amr_results_beta) %>%
  rownames_to_column(var = "Genes") %>%
  dplyr::arrange(padj)                                 

fdr_deseq_amr_beta <- deseq_amr_beta %>%
  dplyr::filter(padj < 0.05)

pvalue_deseq_amr_beta <- deseq_amr_beta %>%
  dplyr::filter(pvalue< 0.05)

dim(pvalue_deseq_amr_beta)


##Metformin ##

deseq_amr_met <- data.frame(desq_amr_results_met) %>%
  rownames_to_column(var = "Genes") %>%
  dplyr::arrange(padj)                                 

fdr_deseq_amr_met <- deseq_amr_met %>%
  dplyr::filter(padj < 0.05)

pvalue_deseq_amr_met <- deseq_amr_met %>%
  dplyr::filter(pvalue< 0.05)

dim(pvalue_deseq_amr_met)


## Diclofenac  ##


deseq_amr_diclo <- data.frame(desq_amr_results_diclo) %>%
  rownames_to_column(var = "Genes") %>%
  dplyr::arrange(padj)                                 

fdr_deseq_amr_diclo <- deseq_amr_diclo %>%
  dplyr::filter(padj < 0.05)

pvalue_deseq_amr_diclo <- deseq_amr_diclo %>%
  dplyr::filter(pvalue< 0.05)

dim(pvalue_deseq_amr_diclo)

######################## DeSeq2 Taxonomy ##################################


diclo_kraken.ps <- subset_samples(kraken_microbiome_small.ps, Sample_type=="diclofenac" | Sample_type =="ciprofloxacin")

#check
head(sample_data(diclo_kraken.ps)$Sample_type, 250)

#Create a DESeq2 object
desq_microbiome_diclo = phyloseq_to_deseq2(diclo_kraken.ps, ~ Sample_type)

#Test this 
desq_microbiome_diclo = DESeq(desq_microbiome_diclo, test="Wald", fitType="local")

#Look at the results
desq_microbiome_results_diclo = results(desq_microbiome_diclo, cooksCutoff = FALSE)

desq_microbiome_results_diclo <-  as.data.frame(desq_microbiome_results_diclo)

plotDispEsts(desq_microbiome_diclo)

deseq_microbiome_diclo_df <- data.frame(desq_microbiome_results_diclo) %>%
  rownames_to_column(var = "species") %>%
  dplyr::arrange(padj)                                 

fdr_microbiome_diclo <- deseq_microbiome_diclo_df %>%
  dplyr::filter(padj < 0.05)

dim(fdr_microbiome_diclo)
#151

# create a function to split each string by an identifier and grab the second index
quick_strsplit <- function(x, index, split){
  # split the string
  temp <- strsplit(x, split = split)
  # iterate through the list and grab each indexed element
  temp <- sapply(temp, "[", index)
  return(temp)
}

#Split up the species column

fdr_microbiome_diclo <- fdr_microbiome_diclo %>%
  mutate(species_only = quick_strsplit(fdr_microbiome_diclo$species, index = 8, split = "\\|"))


fdr_microbiome_diclo <- fdr_microbiome_diclo %>%
  mutate(genus = quick_strsplit(fdr_microbiome_diclo$species, index = 7, split = "\\|"))

fdr_microbiome_diclo <- fdr_microbiome_diclo %>%
  mutate(family = quick_strsplit(fdr_microbiome_diclo$species, index = 6, split = "\\|"))

fdr_microbiome_diclo <- fdr_microbiome_diclo %>%
  mutate(order = quick_strsplit(fdr_microbiome_diclo$species, index = 5, split = "\\|"))


#Add column for plotting

fdr_microbiome_diclo <- fdr_microbiome_diclo %>%
  mutate(gene_type = case_when(log2FoldChange > 0 ~ "up",
                               log2FoldChange < 0 ~ "down" )) 

logfold_diclo <- ggplot(subset(fdr_microbiome_diclo, species_only!="NA"), 
                        aes(x = reorder(species_only, log2FoldChange), y = log2FoldChange, colour = gene_type)) +
  geom_point(size = 4) +
  labs(y = "\nLog2 Fold-Change", x = "") +
  theme(legend.position = "none") +
  coord_flip() +
  theme_classic()  +
  theme(axis.text.y = element_text(face = 'italic'))+
  scale_x_discrete(limits=rev)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "none") +  
  scale_color_manual(values = c("#00CFC1", "#FD5200"))


logfold_diclo


#### Metformin ###
met_kraken.ps <- subset_samples(kraken_microbiome_small.ps, Sample_type=="metformin" | Sample_type =="ciprofloxacin")

#check
head(sample_data(met_kraken.ps)$Sample_type, 250)

#Create a DESeq2 object
desq_microbiome_met = phyloseq_to_deseq2(met_kraken.ps, ~ Sample_type)

#Test this 
desq_microbiome_met = DESeq(desq_microbiome_met, test="Wald", fitType="local")

#Look at the results
desq_microbiome_results_met = results(desq_microbiome_met, cooksCutoff = FALSE)

desq_microbiome_results_met <-  as.data.frame(desq_microbiome_results_met)
#Lots are signiciant lol


plotDispEsts(desq_microbiome_met)

deseq_microbiome_met_df <- data.frame(desq_microbiome_results_met) %>%
  rownames_to_column(var = "taxonomy") %>%
  dplyr::arrange(padj)                                 

fdr_microbiome_met <- deseq_microbiome_met_df %>%
  dplyr::filter(padj < 0.05)

dim(fdr_microbiome_met)
#157

#Split up the species column

fdr_microbiome_met <- fdr_microbiome_met %>%
  mutate(species_only = quick_strsplit(fdr_microbiome_met$taxonomy, index = 8, split = "\\|"))


fdr_microbiome_met <- fdr_microbiome_met %>%
  mutate(genus = quick_strsplit(fdr_microbiome_met$taxonomy, index = 7, split = "\\|"))

fdr_microbiome_met <- fdr_microbiome_met %>%
  mutate(family = quick_strsplit(fdr_microbiome_met$taxonomy, index = 6, split = "\\|"))

fdr_microbiome_met <- fdr_microbiome_met %>%
  mutate(order = quick_strsplit(fdr_microbiome_met$taxonomy, index = 5, split = "\\|"))


#Add column for plotting

fdr_microbiome_met <- fdr_microbiome_met %>%
  mutate(gene_type = case_when(log2FoldChange > 0 ~ "up",
                               log2FoldChange < 0 ~ "down" )) 



logfold_met <- ggplot(subset(fdr_microbiome_met, species_only!="NA"), 
                      aes(x = reorder(species_only, log2FoldChange), y = log2FoldChange, colour = gene_type)) +
  geom_point(size = 4) +
  labs(y = "\nLog2 Fold-Change", x = "") +
  theme(legend.position = "none") +
  coord_flip() +
  theme_classic()  +
  theme(axis.text.y = element_text(face = 'italic'))+
  scale_x_discrete(limits=rev)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "none") +  
  scale_color_manual(values = c("#00CFC1", "#FD5200"))




#### 17-beta-estradiol ###
beta_kraken.ps <- subset_samples(kraken_microbiome_small.ps, Sample_type=="beta-estradiol" | Sample_type =="ciprofloxacin")

#check
head(sample_data(beta_kraken.ps)$Sample_type, 250)

#Create a DESeq2 object
desq_microbiome_beta = phyloseq_to_deseq2(beta_kraken.ps, ~ Sample_type)

#Test this 
desq_microbiome_beta = DESeq(desq_microbiome_beta, test="Wald", fitType="local")

#Look at the results
desq_microbiome_results_beta = results(desq_microbiome_beta, cooksCutoff = FALSE)

desq_microbiome_results_beta <-  as.data.frame(desq_microbiome_results_beta)

plotDispEsts(desq_microbiome_beta)

deseq_microbiome_beta_df <- data.frame(desq_microbiome_results_beta) %>%
  rownames_to_column(var = "taxonomy") %>%
  dplyr::arrange(padj)                                 

fdr_microbiome_beta <- deseq_microbiome_beta_df %>%
  dplyr::filter(padj < 0.05)

dim(fdr_microbiome_beta)
#110

#Split up the species column
fdr_microbiome_beta <- fdr_microbiome_beta %>%
  mutate(species_only = quick_strsplit(fdr_microbiome_beta$taxonomy, index = 8, split = "\\|"))

fdr_microbiome_beta <- fdr_microbiome_beta %>%
  mutate(genus = quick_strsplit(fdr_microbiome_beta$taxonomy, index = 7, split = "\\|"))

fdr_microbiome_beta <- fdr_microbiome_beta %>%
  mutate(family = quick_strsplit(fdr_microbiome_beta$taxonomy, index = 6, split = "\\|"))

fdr_microbiome_beta <- fdr_microbiome_beta %>%
  mutate(order = quick_strsplit(fdr_microbiome_beta$taxonomy, index = 5, split = "\\|"))


#Add column for plotting

fdr_microbiome_beta <- fdr_microbiome_beta %>%
  mutate(gene_type = case_when(log2FoldChange > 0 ~ "up",
                               log2FoldChange < 0 ~ "down" )) 



logfold_beta <- ggplot(subset(fdr_microbiome_beta, species_only!="NA"), 
                       aes(x = reorder(species_only, log2FoldChange), y = log2FoldChange, colour = gene_type)) +
  geom_point(size = 4) +
  labs(y = "\nLog2 Fold-Change", x = "") +
  theme(legend.position = "none") +
  coord_flip() +
  theme_classic()  +
  theme(axis.text.y = element_text(face = 'italic'))+
  scale_x_discrete(limits=rev)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "none") +  
  scale_color_manual(values = c("#00CFC1", "#FD5200"))



################### Put the graphs together log fold microbione ###############


logfold_diclo + logfold_met + logfold_beta +   plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = 'A')

ggsave("figures/Supplementary Figure 5", height = 28, width = 22, dpi = 300)




######################### Comparing the three ##################################

#We have 
fdr_microbiome_diclo
fdr_microbiome_met
fdr_microbiome_beta

#We want to add these together. But to do that we'll need to alter the column names slightly and make the dfs smaller

deseq2_microbiome_diclo <- fdr_microbiome_diclo %>%
  dplyr::select(species_only, log2FoldChange, padj, gene_type, order) %>%
  na.omit(deseq2_microbiome_diclo) 

#137 species 

deseq2_microbiome_diclo$nad <-  c("diclofenac")


deseq2_microbiome_met <- fdr_microbiome_met %>%
  dplyr::select(species_only, log2FoldChange, padj, gene_type, order) %>%
  na.omit(deseq2_microbiome_met)
#140 species

deseq2_microbiome_met$nad <-  c("metformin")



deseq2_microbiome_beta <- fdr_microbiome_beta %>%
  dplyr::select(species_only, log2FoldChange, padj, gene_type, order) %>%
  na.omit(deseq2_microbiome_beta)
#97

deseq2_microbiome_beta$nad <-  c("beta-estradiol")


### Join these dataframes together?

deseq2_all <- rbind(deseq2_microbiome_diclo, deseq2_microbiome_met, deseq2_microbiome_beta)

#Make a new dataframe of species that occur in multiple mixtures
deseq2_all_multiple <-  deseq2_all %>% group_by(species_only) %>% filter(n()>2) 

unique(deseq2_all_multiple$species_only)
#106 things in at least two mixture types
#43 things in all the mixture types


deseq2_all_multiple$nad <- factor(deseq2_all_multiple$nad, levels=c("diclofenac", "metformin", "beta-estradiol"))

#Make a facet label
nad.labs <- c("Diclofenac", "Metformin", "17-β-estradiol")
names(nad.labs) <- c("diclofenac", "metformin", "beta-estradiol")

ggplot(deseq2_all_multiple, 
       aes(x = reorder(species_only, log2FoldChange), y = log2FoldChange, colour = gene_type)) +
  geom_point(size = 4) +
  labs(y = "\nLog2 Fold-Change", x = "") +
  theme(legend.position = "none") +
  coord_flip() +
  theme_classic()  +
  theme(axis.text.y = element_text(face = 'italic'))+
  scale_x_discrete(limits=rev)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.position = "none") +  
  facet_wrap(~nad, labeller = labeller(nad = nad.labs))+
  scale_color_manual(values = c("#00CFC1", "#FD5200"), 
                     labels=c("Diclofenac", "Metformin", "17-β-estradiol"))


ggsave("figures/Figure 7.tiff", height = 14, width = 14, dpi = 300)
