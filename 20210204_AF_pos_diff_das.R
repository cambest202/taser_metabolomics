# ## AF baseline analysis and association of metabolites with DAS44. Positive Ion Mode

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
library(HDDesign)
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
library(qvalue)
library(pscl)
library(parallel)
library(doParallel)

####-------

setwd('/Users/cameronbest/University/MRC_DTP_in_Precision_Medicine/Project/RA/Taser_data/Analysis/taser_metabolomics')

### files------
resp_ints <- read.csv('20210204_AF_resp_int_pos.csv', header=TRUE)
resp_diff <- read.csv('20210204_AF_resp_diff_pos.csv', header=TRUE)
#peak_IDs <- read.table (file="20200117_Taser_PeakIDs.csv", header=TRUE, row.names=1, sep= "\t")
#peak_ID_HMDB <- read.csv(file='peak_IDs_HMDB.csv', header=TRUE, row.names=1)
peak_metadata <- read.csv(file='20200521_Taser_POS_Peakdata.csv', header=TRUE, row.names=1)

sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)

## Add uncertainty around peak identification
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'DG', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'LysoPC', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'LysoPE', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'PI', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'PC', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'GPC', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)
peak_metadata$Putative_Metabolite <- ifelse(peak_metadata$Putative_Metabolite == 'PE', paste0(peak_metadata$Putative_Metabolite,'?'),peak_metadata$Putative_Metabolite)

### Differential Analysis: Metabolic changes between positive and negative responders ----------
resp_PN_ints <- subset(resp_ints, resp_ints$time =='A')
resp_PN_ints <- resp_PN_ints[12:1596]
resp_PN_ints$DAS44_Response[resp_PN_ints$DAS44_Response %like% 'Good'|resp_PN_ints$DAS44_Response %like% 'Remission'] <- 'Positive'
resp_PN_ints$DAS44_Response[resp_PN_ints$DAS44_Response %like% 'Poor'|resp_PN_ints$DAS44_Response %like% 'Mild'] <- 'Negative'
resp_PN_ints$DAS44_Response <- as.factor(resp_PN_ints$DAS44_Response)

resp_PN_limma <- resp_PN_ints
rownames(resp_PN_limma) <- paste0(rownames(resp_PN_limma), resp_PN_limma$DAS44_Response)
resp_PN_limma <- resp_PN_limma[,-1]
resp_PN_limma_t <- as.data.frame(t(resp_PN_limma))
names(resp_PN_limma_t)
names(resp_PN_limma_t)[names(resp_PN_limma_t) %like% 'Positive'] <- 'Positive'
names(resp_PN_limma_t)[names(resp_PN_limma_t) %like% 'Negative'] <- 'Negative'
resp_PN_limma_t <- as.data.frame(resp_PN_limma_t)

#Limma----
Group <- factor(colnames(resp_PN_limma_t), levels = c('Positive', 'Negative'))
design <- model.matrix (~Group)
colnames(design) <- c('Positive', 'PositivevsNegative')
eset <- resp_PN_limma_t
fit <- lmFit(eset, design)
fit <- eBayes(fit)
toptable <- topTable(fit, coef = 'PositivevsNegative', adjust = 'BH', number = 1584)
toptable <- as.data.frame(toptable)
toptable$Peak_ID <- rownames(toptable)
toptable$Peak_ID <- gsub('X', '', toptable$Peak_ID)
toptable <- toptable[,c(ncol(toptable),1:(ncol(toptable)-1))]
toptable <- join(toptable, peak_metadata, by = 'Peak_ID')
##-----
toptable$Sig <- 0
toptable$Sig <- ifelse(toptable$adj.P.Val <0.05, 'Significant', 'Not significant') 
#toptable$Sig_Peaks <- ifelse(toptable$Putative_Metabolite %in% LR_mets, toptable$Putative_Metabolite, '')
toptable$Sig_Peaks <- ifelse(toptable$adj.P.Val <0.05 & toptable$identification != '', toptable$Putative_Metabolite, '')

toptable_dist <- distinct(toptable, Putative_Metabolite, .keep_all = TRUE)

ggplot(data=toptable_dist, aes(x=logFC, y=-log10(P.Value), 
                                    colour=Sig, 
                                    group=Sig)) +
  geom_point (alpha=0.7) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Signficance')+
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                  box.padding =1,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.0015, "npc"))) +  
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  scale_color_brewer(palette = "Set1")


### Differential Analysis: Metabolic changes between baseline and 18 months ----------
resp_AF_ints <- resp_ints[,c(3,13:1596)]
rownames(resp_AF_ints) <- paste0(rownames(resp_AF_ints), resp_AF_ints$time)
resp_AF_limma <- resp_AF_ints[,-1]
resp_AF_limma_t <- as.data.frame(t(resp_AF_limma))
names(resp_AF_limma_t)
names(resp_AF_limma_t)[names(resp_AF_limma_t) %like% 'A'] <- 'A'
names(resp_AF_limma_t)[names(resp_AF_limma_t) %like% 'F'] <- 'F'
resp_AF_limma_t <- as.data.frame(resp_AF_limma_t)

#Limma----
Group_AF <- factor(colnames(resp_AF_limma_t), levels = c('A', 'F'))
design_AF <- model.matrix (~Group_AF)
colnames(design_AF) <- c('A', 'AvsF')
eset_AF <- resp_AF_limma_t
fit_AF <- lmFit(eset_AF, design_AF)
fit_AF <- eBayes(fit_AF)
toptable_AF <- topTable(fit_AF, coef ='AvsF', adjust = 'BH', number = 1584)
toptable_AF <- as.data.frame(toptable_AF)
toptable_AF$Peak_ID <- rownames(toptable_AF)
toptable_AF$Peak_ID <- gsub('X', '', toptable_AF$Peak_ID)
toptable_AF <- toptable_AF[,c(ncol(toptable_AF),1:(ncol(toptable_AF)-1))]
toptable_AF$Peak_ID <- as.numeric(toptable_AF$Peak_ID)
toptable_AF <- join(toptable_AF, peak_metadata, by = 'Peak_ID')
##-----
select_mets_syno <-c('Inosine', 'L-carnitine', 'L-Glutamine', 'Pyroglutamate', 'Thymine', 'L-valine')
select_mets_syno_ID <- c(58, 23, 1232, 1069, 57, 1015)

select_mets_das_ID <- c(602,6,168,558,303)
select_mets_das_ID_all <- c(1038,243,400,25,423,1076)

toptable_AF$Sig <- 0
toptable_AF$Sig <- ifelse(toptable_AF$adj.P.Val <0.05, '< 0.05', '> 0.05') 
#toptable_AF$Sig_Peaks <- ifelse(toptable_AF$Peak_ID %in% select_mets_syno_ID, toptable_AF$Putative_Metabolite, '')
toptable_AF$Sig_Peaks <- ifelse(toptable_AF$P.Value<0.01 & toptable_AF$identification != '', toptable_AF$Putative_Metabolite, '')

toptable_AF_ID <- subset(toptable_AF,toptable_AF$Putative_Metabolite != '')
toptable_AF_sig <- subset(toptable_AF, toptable_AF$adj.P.Val <0.05)
toptable_AF_sig_dist <- distinct(toptable_AF_sig, Putative_Metabolite, .keep_all = TRUE)

#write.csv(toptable_AF_ID, '20210218_AF_pos_diff.csv')

ggplot(data=toptable_AF_ID, aes(x=logFC, y=-log10(P.Value), 
                          colour=Sig, 
                          group=Sig)) +
  geom_point (alpha=0.7) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Adjusted \np-value',
        title= 'Positive Ion Mode')+
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                  box.padding =1.2,
                  size=2.5,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.005, "npc"))) +  
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  scale_color_brewer(palette = "Set1",direction=-1)+
  ylim(0,4)+
  xlim(-1,1)

 ## P-value histogram---
# Adding local FDR, qvalue and π0 values to p-value histogram
pi0 <- 2*mean(toptable_AF$P.Value > 0.05)
lfdrvals <- lfdr(toptable_AF$P.Value, pi0)
qobj <- qvalue(toptable_AF$P.Value)
hist(qobj)

# PCA of samples
scaled_intensities <- scale(t(resp_AF_limma_t))
scaled_intensities[do.call(cbind, lapply(scaled_intensities, is.nan))] <- 0

pca_data <- prcomp(scaled_intensities)
pca_coord <- data.frame(pca_data$x)
var_explained <- pca_data$sdev^2/sum(pca_data$sdev^2)
var_explained[1:5]

pca_coord$group <-as.factor(row.names(pca_coord))
head(pca_coord$group )
pca_coord$group[pca_coord$group %like% 'A'] <- 'A'
pca_coord$group[pca_coord$group %like% 'F'] <- 'F'

ggplot(pca_coord) + 
  geom_point(size=3, 
             aes(x=PC1,y=PC2, colour= group, fill= group))+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
  theme_classic2()


### Correlations of the differential metabolites with the DAS44
resp_diff_melt <- resp_diff[,c(2,13:1597)]
resp_diff_melt <- melt(resp_diff_melt)
names(resp_diff_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_diff_melt$DAS44 <- resp_diff$DAS44
resp_diff_melt$Peak_ID <- gsub('X','', resp_diff_melt$Peak_ID)
resp_diff_melt$CRP <- resp_diff$CRP
resp_diff_melt$ESR <- resp_diff$ESR
resp_diff_melt$HAQ <- resp_diff$HAQ
resp_diff_melt$GHVAS <- resp_diff$GHVAS
resp_diff_melt$PVAS <- resp_diff$PVAS
resp_diff_melt$RAI <- resp_diff$RAI
resp_diff_melt$SJC <- resp_diff$SJC

# DAS44----
resp_melt_top <- subset(resp_diff_melt, resp_diff_melt$Peak_ID %in% toptable_AF_sig$Peak_ID)
# make sure all disease measures are included in the melted df

ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

top_20_peaks <- head(bestest_fit, 20)
top_50_peaks <- head(bestest_fit, 50)
top_100_peaks <- head(bestest_fit, 100)
top_200_peaks <- head(bestest_fit, 200)

best_augmented <- bestest_fit %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
sig_peaks <- inner_join(best_adj, peak_metadata, by='Peak_ID')
sig_peaks<- subset(sig_peaks,sig_peaks$Putative_Metabolite != '')
sig_peaks <- distinct(sig_peaks, Peak_ID, .keep_all = TRUE)

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
best_adj_hmdb$Peak_ID2 <- best_adj_hmdb$Peak_ID
best_adj_hmdb<- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID != '168')
das_metabolites <-best_adj_hmdb 
colnames(das_metabolites)[17] <-'Disease_Measure'
das_metabolites$Disease_Measure_Type <- 'DAS'

ggplot(best_adj_hmdb,aes(x = DAS44, y=Peak_Intensity)) +
  geom_point() + 
  stat_cor(method = "pearson", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  
  labs(x='ΔDAS44',
       y='ΔPeak Intensity')+
  theme_minimal()


## ESR-----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~ESR, data = .x)))
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
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.5) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
esr_metabolites <-best_adj_hmdb 
colnames(esr_metabolites)[17] <-'Disease_Measure'
esr_metabolites$Disease_Measure_Type <- 'ESR'
ggplot(best_adj_hmdb,aes(x = ESR, y=Peak_Intensity)) +
  geom_point() + 
  stat_cor(method = "pearson", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  
  labs(x='ΔESR',
       y='ΔPeak Intensity')+
  theme_minimal()

## HAQ-----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~HAQ, data = .x)))
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
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.5) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')

haq_metabolites <-best_adj_hmdb 
colnames(haq_metabolites)[17] <-'Disease_Measure'
haq_metabolites$Disease_Measure_Type <- 'HAQ'

ggplot(good_,aes(x = HAQ, y=Peak_Intensity)) +
  geom_point() + 
  stat_cor(method = "pearson", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  
  labs(x='ΔHAQ',
       y='ΔPeak Intensity')+
  theme_minimal()

## PVAS-----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~PVAS, data = .x)))
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
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.5) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')

pvas_metabolites <-best_adj_hmdb 
colnames(pvas_metabolites)[17] <-'Disease_Measure'
pvas_metabolites$Disease_Measure_Type <- 'PVAS'

ggplot(best_adj_hmdb,aes(x = PVAS, y=Peak_Intensity)) +
  geom_point() + 
  stat_cor(method = "pearson", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  
  labs(x='ΔPVAS',
       y='ΔPeak Intensity')+
  theme_minimal()
## GHVAS-----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~GHVAS, data = .x)))
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
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.5) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)
sig_peaks$ID <- sig_peaks$Peak_ID
sig_peaks <- subset(sig_peaks,sig_peaks$Putative_Metabolite != '')

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')

ghvas_metabolites <-best_adj_hmdb 
colnames(ghvas_metabolites)[17] <-'Disease_Measure'
ghvas_metabolites$Disease_Measure_Type <- 'GHVAS'

ggplot(best_adj_hmdb,aes(x = GHVAS, y=Peak_Intensity)) +
  geom_point() + 
  stat_cor(method = "pearson", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  
  labs(x='ΔGHVAS',
       y='ΔPeak Intensity')+
  theme_minimal()

# Combining the correlation results
mets_reduce <- function(disease_df){
  disease_met <- disease_df[c(31:32)]
  disease_met <- distinct(disease_met, Putative_Metabolite, .keep_all = TRUE)
}

DAS_met <- mets_reduce(das_metabolites)
ESR_met <- mets_reduce(esr_metabolites)
HAQ_met <- mets_reduce(haq_metabolites)
PVAS_met <- mets_reduce(pvas_metabolites)
GHVAS_met <- mets_reduce(ghvas_metabolites)

shared_mets <- rbind.data.frame(das_metabolites,
                                esr_metabolites,
                                haq_metabolites,
                                pvas_metabolites,
                                ghvas_metabolites)

shared_mets_reduced <- rbind.data.frame(DAS_met,
                                        ESR_met,
                                        HAQ_met,
                                        PVAS_met,
                                        GHVAS_met)

metabolite_counts <- as.data.frame(as.matrix(table(shared_mets_reduced$Putative_Metabolite, shared_mets_reduced$Disease_Measure_Type)))
names(metabolite_counts)<- c('Putative_Metabolite', 'Disease_Measure','Frequency')

write.csv(shared_mets_reduced, '20210205_AF_pos_diff_das_mets_count.csv')

metabolite_counts%>%
  ggplot()+
  geom_col(aes(reorder(Putative_Metabolite, Frequency),
               Frequency,
               fill=Disease_Measure))+
  theme_pubclean()+
  coord_flip()+
  labs(x='Putative Metabolite',
       y='Frequency')


## Determine the most important features to predict clinical outcomes
# Reduce number of features through those correlating peaks
resp_split <- resp_diff[c(11,13:1597)]
resp_split <- resp_split[,-1]
resp_split_t <- as.data.frame(t(resp_split))
resp_split_t$Peak_ID <- rownames(resp_split_t)
resp_split_t$Peak_ID <- gsub('X','', resp_split_t$Peak_ID)
resp_split_done <- subset(resp_split_t, resp_split_t$Peak_ID %in% toptable_AF_sig$Peak_ID)
#resp_split_done <- subset(resp_split_done, resp_split_done$Peak_ID %in% toptable_AF_sig$Peak_ID)

resp_red <- resp_split_done[,-64]
resp_red_t <- as.data.frame(t(resp_red))
resp_red_t$Response <- resp_PN_ints$DAS44_Response
resp_red_t$Response <- as.character(resp_red_t$Response)
#resp_red_t$Response[resp_red_t$Response %like% 'Positive'] <- 1
#resp_red_t$Response[resp_red_t$Response %like% 'Negative'] <- 0
#resp_red_t$Response <- as.numeric(resp_red_t$Response)

set.seed(42)
index <- createDataPartition(resp_red_t$Response, p = 0.8, list = FALSE)
train_data <- resp_red_t[index, ]
test_data  <- resp_red_t[-index, ]
test_data$Response <- as.factor(test_data$Response)
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
set.seed(42)
model_rf <- train(Response~.,
                         data = train_data,
                         method = "sv",
                         metric = "Accuracy",
                         tuneGrid=tuneGrid,
                         trControl = trainControl(method = "repeatedcv",
                                                  number =5,
                                                  repeats = 5, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE, 
                                                  allowParallel = TRUE),
                         importance = TRUE,
                         ntree = 500)

plot(model_rf)
print(model_rf)

test_results <- predict(model_rf, newdata = test_data)
summary(test_results)
summary(test_data$Response)
confusionMatrix(test_results, test_data$Response)
# 85% split, number=10, repeats=5 --> 77.8% acc

final <- data.frame(actual = test_data$Response,
                    predict(model_rf, newdata = test_data, type = "prob"))
final$predicted <- 0
final$predicted[final$Positive > final$Negative] <- 'Positive'
final$predicted[final$Positive < final$Negative] <- 'Negative'
final$correct <- 0
final$correct <- ifelse(final$predict == final$actual, 'Correct', 'Incorrect')
summary(final$predicted == final$actual)

imp <- model_rf$finalModel$imprtance
imp <- as.data.frame(imp)
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
imp_peaks <- varImp(model_rf, scale = TRUE)
plot(imp_peaks, top=20)
imps <- as.data.frame(imp_peaks$importance)
imps$Peak_ID <- rownames(imps)
imps <- imps[,-1]
colnames(imps)[1] <- 'Importance'
#imps <- subset(imps,imps$Importance >10)
imps$Peak_ID <- gsub('X','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_metadata, by='Peak_ID')
imps_hmdb <- distinct(imps_hmdb, Putative_Metabolite, .keep_all = TRUE)
imps_hmdb <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='')
# Plot the annotated peaks from feature importance
ggplot(imps_hmdb)+
  geom_col(aes(reorder(Putative_Metabolite, Importance), 
               Importance),
           fill=0x3a5e84,
           colour='black')+
  coord_flip()+
  theme_minimal()+
  labs(y='Relative Importance',
       x='Putative Metabolite')+
  theme(axis.text.y = element_text(size = 8))

### Logistic regression model
resp_split <- resp_PN_ints
resp_split <- resp_split[,-1]
resp_split_t <- as.data.frame(t(resp_split))
resp_split_t$Peak_ID <- rownames(resp_split_t)
resp_split_t$Peak_ID <- gsub('X','', resp_split_t$Peak_ID)
toptable_AF_sig_FC <- subset(toptable_AF_sig, toptable_AF_sig$adj.P.Val < 0.05 & toptable_AF_sig$Sig_Peaks != '')
toptable_AF_sig_FC$logFCsqrd <- toptable_AF_sig_FC$logFC *toptable_AF_sig_FC$logFC
toptable_AF_sig_FC <- distinct(toptable_AF_sig_FC, Putative_Metabolite, .keep_all = TRUE)
toptable_AF_sig_FC <- subset(toptable_AF_sig_FC, toptable_AF_sig_FC$logFCsqrd > 0.5)
resp_split_done <- subset(resp_split_t, resp_split_t$Peak_ID %in% toptable_AF_sig_FC$Peak_ID)
#resp_split_done <- subset(resp_split_t, resp_split_t$Peak_ID %in% imps_hmdb$Peak_ID)

resp_red <- resp_split_done[,-64]
resp_red_t <- as.data.frame(t(resp_red))
resp_red_t$Response <- resp_PN_ints$DAS44_Response
resp_red_t$Response <- as.character(resp_red_t$Response)
resp_red_t$Response[resp_red_t$Response=='Positive'] <- 1
resp_red_t$Response[resp_red_t$Response=='Negative'] <- 0
resp_red_t$Response <- as.integer(resp_red_t$Response)

set.seed(42)
index <- createDataPartition(resp_red_t$Response, p = 0.8, list = FALSE)
train_data <- resp_red_t[index, ]
test_data  <- resp_red_t[-index, ]

model <- glm(Response ~    X157 + X168 +  X322 + X414 ,
             family=binomial(link='logit'),
             data=train_data)

?glm
summary(model)
anova(model, test="Chisq")
pR2(model)
fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_data$Response)
print(paste('Accuracy',1-misClasificError))

1-pchisq(58.086-46.034, 53-43)

library(ROCR) 
p <- predict(model,newdata=test_data,type='response')
pr <- prediction(p, test_data$Response)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

## Assess selected metabolites for predicting response
resp_diff_melt2 <- resp_diff_melt
resp_diff_melt2$Response <- resp_PN_ints$DAS44_Response
select <- c('932', '168', '1137', '157','37')
resp_diff_melt2 %>%
  subset(Peak_ID %in% select) %>%
  ggplot(aes(Response, Peak_Intensity,
             fill=Response))+
  theme_light()+
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  stat_compare_means(method= 'wilcox.test',
                     label = "p.format",
                     vjust=1, 
                     hjust=-1)+
  theme(legend.position = 'none')+
  labs(y='Peak Intensity') 
