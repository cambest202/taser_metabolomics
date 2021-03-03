### Investigating the outlier metabolites which show 
###libraries------
library(ggplot2)
library(tidyverse) 
library (readr)
library(PMCMRplus)
library (PMCMR)
library(ggpubr)
library (mosaic)
library (dplyr)
library (data.table)
library(reshape2)
library (gtools)
library(plyr)
library(limma)
library(ggrepel)
library(amap)
library(rstatix)
library(broom)
library(ggprism)
library(parallel)
library(doParallel)
library(caret)
library(rsample)
library(sandwich)
library(rpart)
library(rpart.plot)
library(randomForest)
library(RColorBrewer)
library(plotly)
library(purrr)
library(devtools)
library(e1071)
library(ggraph)
library(igraph)


# high levels in association with greatest reductions in DAS44 in samples
resp_ints <- read.csv('20201105_resp_ints.csv', header=TRUE)
peak_metadata <- read.csv(file='20210217_neg_peak_metadata.csv', header=TRUE)
sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)

patient_metadata$Sample_Name <- rownames(patient_metadata)
names(sample_sheet)[2] <- 'Sample_Name'
resp_ints_A <- subset(resp_ints,resp_ints$Sample == 'A')
resp_A_melt <- resp_ints_A[,c(2,9:1466)]
resp_A_melt <- melt(resp_A_melt)
names(resp_A_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_A_melt$DAS44 <- resp_ints_A$DAS44
resp_A_melt$Peak_ID <- gsub('X','', resp_A_melt$Peak_ID)
resp_A_melt$Peak_ID <- as.numeric(resp_A_melt$Peak_ID)
resp_A_melt <- inner_join(resp_A_melt,peak_metadata, by='Peak_ID')
resp_A_melt <- subset(resp_A_melt, resp_A_melt$Putative_Metabolite != 'NA')

# Linear models for Peak intensities and ΔDAS44
resp_melt_top <- resp_A_melt
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~DAS44, data = .x)))
model_coef_nested <- ints_lm %>%
  mutate(coef=map(model, ~tidy(.x)))
model_coef <- model_coef_nested %>%
  unnest(coef)
model_perf_nested <- ints_lm %>%
  mutate(fit=map(model, ~glance(.x)))
model_perf <- model_perf_nested%>%
  unnest(fit)
best_fit <- model_perf %>%
  top_n(n=4, wt=r.squared)
bestest_fit <- with(model_perf,model_perf[order(-r.squared),])

best_augmented <- bestest_fit %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)

best_augmented <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented$Peak_ID <- as.numeric(best_augmented$Peak_ID)
adj_p <- p.adjust(best_augmented$p.value, method='BH')
best_augmented$adj_p <- adj_p
best_augmented$Peak_ID <- as.numeric(best_augmented$Peak_ID)
best_augmented_ID <- inner_join(best_augmented, peak_metadata, by='Peak_ID')
best_augmented_ID <- subset(best_augmented_ID, best_augmented_ID$Putative_Metabolite != 'NA')

best_augmented_ID%>%
  subset(Putative_Metabolite != 'Traumatic acid')%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(size=1, alpha=0.7) + 
  stat_cor(method = "spearman", 
           vjust=1, hjust=0,
           size=2.5)+
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=8))+
  labs(x='ΔDAS44',
       y='Peak Intensity')+
  theme_minimal()

# Select the metabolites of interest for outliers
resp_melt_outlier <- inner_join(resp_melt_top, patient_metadata, by= 'Sample_Name') 
resp_melt_outlier$Peak_ID <- as.numeric(resp_melt_outlier$Peak_ID)

resp_outlier_sig <- subset(resp_melt_outlier, resp_melt_outlier$Peak_ID %in% best_augmented_ID$Peak_ID)
resp_outlier_sig <- resp_outlier_sig%>%
  subset(Putative_Metabolite != "Fumarate?")%>%
  subset(Putative_Metabolite != 'Traumatic acid')
  
# initial plot
resp_outlier_sig%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(size=1, alpha=0.7) + 
  stat_cor(method = "spearman", 
           vjust=8, hjust=0,
           size=2.5)+
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=8))+
  labs(x='ΔDAS44',
       y='Peak Intensity')+
  theme_minimal()

resp_outlier_sig$Sample_ID <- 0
resp_outlier_sig$Sample_ID <- ifelse(resp_outlier_sig$Peak_Intensity >18, resp_outlier_sig$Sample_Name,'')

# 
resp_outlier_sig%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(aes(colour=Smoking),
             size=1, alpha=0.7) + 
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=8))+
  geom_vline(xintercept= -2, linetype="dotted")+
  geom_hline(yintercept= 17.5, linetype="dotted")+
  labs(x='ΔDAS44',
       y='Peak Intensity')+
  geom_text_repel(aes(x = DAS44, y=Peak_Intensity, label = Sample_ID),
                  box.padding =1,
                  size=3,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.0015, "npc"))) +    
  theme_minimal()

resp_outlier_dist <- distinct(resp_outlier_sig,Peak_ID, .keep_all = TRUE)

resp_outlier_sig%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(size=1, alpha=0.7) + 
  stat_cor(method = "spearman", 
           vjust=5, hjust=-0.7,
           size=3)+
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=8))+   
  geom_hline(yintercept= 17.5, linetype="dotted")+
  labs(x='ΔDAS44',
       y='Peak Intensity')+
    theme_minimal()

resp_outlier_sig%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(aes(colour=Smoking),
             size=2, alpha=0.5) + 
  
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=8))+
  geom_vline(xintercept= -2, linetype="dotted")+
  geom_hline(yintercept= 17.5, linetype="dotted")+
  labs(x='ΔDAS44',
       y='Peak Intensity',
       fill='Smoker status')+
  scale_colour_manual(values=c("#009900", "#CC0099", "#0066FF"),
                      name="Smoker status",
                    breaks=c("C", "F", "N"),
                    labels=c("Current", "Former", "Non-Smoker"))+
  theme_minimal()

resp_outlier_sig$Sym_Dur <- 0
resp_outlier_sig$Sym_Dur <- ifelse(resp_outlier_sig$SymptomDuration >9, '> 9', '< 9')

resp_outlier_sig%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(aes(colour= Group),
             size=2, alpha=0.5) + 
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=8))+
  geom_vline(xintercept= -2, linetype="dotted")+
  geom_hline(yintercept= 17.5, linetype="dotted")+
  labs(x='ΔDAS44',
       y='Peak Intensity',
       fill='Smoker status')+
  theme_minimal()
