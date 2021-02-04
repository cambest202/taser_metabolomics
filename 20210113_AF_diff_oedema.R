## AF differential analysis and association of metabolites with MRI Oedema Score Outcomes

### ggplot theme------
theme_fish <- function () { 
  theme_bw() %+replace% 
    theme(
      panel.grid.major = element_line(colour = "black"),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.background = element_blank(), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      plot.title = element_text(size=16, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial"),
      title = element_text(size = 16, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial"),
      axis.text.y = element_text(size = 14, margin = margin(r = 5),hjust=1,vjust=0.5, family="Arial",colour="black"),
      axis.text.x = element_text(size = 14, margin = margin(t = 5),hjust=0.5,vjust=1, family="Arial",colour="black"), 
      axis.title.y = element_text(size = 16, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, family="Arial", face="bold"),
      axis.title.x = element_text(size = 16, margin = margin(t = 10),hjust=0.5,vjust=1, family="Arial", face="bold"),
      legend.text=element_text(size=14, family="Arial", face="bold"),
      legend.title=element_blank(), 
      legend.key.size=unit(2.5,"line"),
      plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm")
    )
}

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

setwd('/Users/cameronbest/University/MRC_DTP_in_Precision_Medicine/Project/RA/Taser_data/Analysis/taser_metabolomics')

### files------
resp_ints <- read.csv('20201105_resp_ints.csv', header=TRUE)
resp_diff <- read.csv('20201105_resp_diff.csv', header=TRUE)
peak_IDs <- read.table (file="20200117_Taser_PeakIDs.csv", header=TRUE, row.names=1, sep= "\t")
peak_ID_HMDB <- read.csv(file='peak_IDs_HMDB.csv', header=TRUE, row.names=1)
peak_metadata <- read.csv(file='20200427_Taser_PeakMetadata.csv', header=TRUE, row.names=1)
sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
mri_erosion <- read.csv('20201202_MRI_mean.csv', header=TRUE, na.strings=c("","NA"))
xray_erosion <- read.csv('20201202_xray_mean.csv', header=TRUE, na.strings=c("","NA"))

peak_ID_HMDB[13,5] <- '4-Hydroxybutyric acid'
names(sample_sheet)[2] <- 'Sample_Name'
peak_IDs$Peak_ID <- rownames(peak_IDs)
peak_IDs$Peak_ID  <- as.numeric(peak_IDs$Peak_ID )
peak_IDs$moleculeName  <- as.character(peak_IDs$moleculeName )
colnames(peak_metadata)[6] <- 'Peak_ID'
peak_metadata <- left_join(peak_metadata, peak_IDs, by='Peak_ID')
names(sample_sheet)[2] <- 'Sample_Name'

peak_metadata_IDs <- subset(peak_metadata, peak_metadata$identification != '')
peak_metadata$Peak_ID_test <- peak_metadata$Peak_ID

# Which metabolites over the 18 months are associated with the erosion imaging score outcomes after 18 months of treatment?
# Xray and MRI should be done independently
# MRI Oedema
names(mri_erosion)[1] <- 'Sample_Name'
mri_resp <- inner_join(mri_erosion, resp_diff, by='Sample_Name')
mri_melt <- mri_resp[,c(1,17:1474)]
mri_melt <- melt(mri_melt)
mri_melt$DAS44 <- mri_resp$DAS44
mri_melt$MRI_Oedema <- mri_resp$ΔMRI_Oedema
mri_melt$CRP <- mri_resp$CRP
mri_melt$ESR <- mri_resp$ESR
mri_melt$HAQ <- mri_resp$HAQ
mri_melt$GHVAS <- mri_resp$GHVAS
mri_melt$PVAS <- mri_resp$PVAS
mri_melt$RAI <- mri_resp$RAI
mri_melt$SJC <- mri_resp$SJC
names(mri_melt)[1:3] <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')

mri_resp$Oedema_Response <- 0
mri_pre_limma <- mri_resp[, c(4,1475,17:1474)]
mean(mri_pre_limma$ΔMRI_Oedema)
median(mri_pre_limma$ΔMRI_Oedema)
histogram(mri_pre_limma$ΔMRI_Oedema)

mri_pre_limma$Oedema_Response[mri_pre_limma$ΔMRI_Oedema >=-2] <- 'Negative'
mri_pre_limma$Oedema_Response[mri_pre_limma$ΔMRI_Oedema < -2] <- 'Positive'
mri_pre_limma_2 <- mri_pre_limma

mri_pre_limma %>%
  ggplot(aes(length(Oedema_Response),
             fill=Oedema_Response))+
  geom_histogram()+
  facet_wrap(~Oedema_Response)+
  theme(legend.position = 'none')

mri_pre_limma <- mri_pre_limma[,-1]
rownames(mri_pre_limma) <- paste0(rownames(mri_pre_limma), mri_pre_limma$Oedema_Response)
mri_pre_limma <- mri_pre_limma[,-1]
mri_pre_limma_t <- t(mri_pre_limma)

colnames(mri_pre_limma_t)[colnames(mri_pre_limma_t) %like% 'Negative'] <- 'Negative'
colnames(mri_pre_limma_t)[colnames(mri_pre_limma_t) %like% 'Positive'] <- 'Positive'
rownames(mri_pre_limma_t) <- gsub('X', '', rownames(mri_pre_limma_t))
### Differential analysis of baseline metabolites from above across response groups. Which metabolites are differentially produced?------
## Limma functions----------------------------------------------------------------------------------------
limma_fun <- function(matrix_AB, no., var1, var2){
  Group <- factor(colnames(matrix_AB), levels = c(`var1`, `var2`))
  design <- model.matrix (~Group)
  colnames(design) <- c('var1', 'var1vsvar2')
  eset <- matrix_AB
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  toptable <- topTable(fit, coef = 'var1vsvar2', adjust = 'BH', number = no.)
  toptable <- as.data.frame(toptable)
  toptable$Peak_ID <- rownames(toptable)
  toptable$Peak_ID <- as.numeric(as.character(toptable$Peak_ID))
  toptable <- toptable[,c(ncol(toptable),1:(ncol(toptable)-1))]
  toptable <- join(toptable, peak_metadata, by = 'Peak_ID')
  toptable <- toptable[, -c(3,4,7,10:12)]
}

mri_limma <- limma_fun(mri_pre_limma_t, 1500, 'Negative', 'Positive')
mri_limma$Sig <- 0
mri_limma$Sig <- ifelse(mri_limma$adj.P.Val <0.05, 1, 0) 
mri_limma$Sig_Peaks <- ifelse(mri_limma$P.Value<0.05 & mri_limma$identification != '', mri_limma$Peak_ID, '')

mri_limma_hmdb <- inner_join(mri_limma, peak_ID_HMDB, by='Peak_ID')
mri_limma_hmdb$Sig_Peaks <- ifelse(mri_limma_hmdb$P.Value<0.05, mri_limma_hmdb$Putative_Metabolite, '')

ggplot(data=mri_limma_hmdb, aes(x=logFC, y=-log10(P.Value), 
                                colour=Sig, 
                                group=Sig)) +
  geom_point (alpha=0.7) +
  theme_light() +
  labs (x='LogFC',
        y='-Log p-value',
        title='MRI Oedema: Differential Analysis of Metabolites \nBetween Oedema Scoring Outcomes') +
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                  box.padding =1,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.0015, "npc"))) +  
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Histogram of p-values 
mri_limma%>%
  mutate(Colour = ifelse(P.Value < 0.05, "P.Value < 0.05", "P.Value > 0.05")) %>%
  ggplot(aes(P.Value))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  labs(x='p-value',
       y='Frequency')

# PCA of samples
scaled_intensities <- scale(t(mri_pre_limma_t))
scaled_intensities[do.call(cbind, lapply(scaled_intensities, is.nan))] <- 0

pca_data <- prcomp(scaled_intensities)
pca_coord <- data.frame(pca_data$x)
var_explained <- pca_data$sdev^2/sum(pca_data$sdev^2)
var_explained[1:5]

pca_coord$group <-as.factor(row.names(pca_coord))
head(pca_coord$group )
pca_coord$group[pca_coord$group %like% 'Positive'] <- 'Positive'
pca_coord$group[pca_coord$group %like% 'Negative'] <- 'Negative'

ggplot(pca_coord) + 
  geom_point(size=3, 
             aes(x=PC1,y=PC2, colour= group, fill= group))+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
  theme_classic2()

# Generate the model
# Set training and testing data
mri_pre_limma_2 <- mri_pre_limma_2[,-1]
set.seed(42)
index <- createDataPartition(mri_pre_limma_2$Oedema_Response, p = 0.85, list = FALSE)
train_data <- mri_pre_limma_2[index, ]
test_data  <- mri_pre_limma_2[-index, ]
tunegrid <- expand.grid(.mtry=c(1:sqrt(1458)))
set.seed(42)
model_rf <- caret::train(Oedema_Response~.,
                         data = train_data,
                         method = "rf",
                         metric = "Accuracy",
                         tuneGrid=tunegrid,
                         trControl = trainControl(method = "repeatedcv",
                                                  number =10,
                                                  repeats = 5, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE, 
                                                  allowParallel = TRUE),
                         importance = TRUE,
                         ntree = 500)


summary(model_rf)

plot(model_rf)

test_results <- predict(model_rf, newdata = test_data)
summary(test_results)
test_data$Oedema_Response <- as.factor(test_data$Oedema_Response)
summary(test_data$Oedema_Response)
confusionMatrix(test_results, test_data$Oedema_Response)

final <- data.frame(actual = test_data$Oedema_Response,
                    predict(model_rf, newdata = test_data, type = "prob"))
final$predicted <- 0
final$predicted[final$Positive > final$Negative] <- 'Positive'
final$predicted[final$Positive < final$Negative] <- 'Negative'

final$correct <- 0
final$correct <- ifelse(final$predict == final$actual, 'Correct', 'Incorrect')
summary(final$predicted == final$actual)

imp <- model_rf$finalModel$importance
imp <- as.data.frame(imp)
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
imp_peaks <- varImp(model_rf, scale = TRUE)
plot(imp_peaks, top=20)
imps <- as.data.frame(imp_peaks$importance)
imps$Peak_ID <- rownames(imps)
imps <- imps[,-1]
colnames(imps)[1] <- 'Importance'
imps <- subset(imps,imps$Importance >30)
imps$Peak_ID <- gsub('X','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_ID_HMDB, by='Peak_ID')
imps_hmdb <- with(imps_hmdb,imps_hmdb[order(Importance),])
imps_hmdb_id <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='NA')
imps_hmdb_id <- distinct(imps_hmdb_id, Putative_Metabolite, .keep_all = TRUE)
imps_hmdb_id$Putative_Metabolite <- as.factor(imps_hmdb_id$Putative_Metabolite)

write.csv(imps_hmdb_id, '20210113_AF_mri_oed_FI.csv')

# Plot the annotated peaks from feature importance
ggplot(imps_hmdb_id)+
  geom_col(aes(reorder(Putative_Metabolite, Importance), 
               Importance),
           fill=0x3a5e84,
           colour='black')+
  coord_flip()+
  theme_minimal()+
  labs(y='Relative Importance',
       x='Putative Metabolite')+
  theme(axis.text.y = element_text(size = 8))

### Build linear models for the top peaks from feature selection for each of the disease measurements
mri_melt$Peak_ID <- gsub('X','', mri_melt$Peak_ID)
mri_melt_top <- subset(mri_melt, mri_melt$Peak_ID %in% peak_ID_HMDB$Peak_ID)
# make sure all disease measures are included in the melted df

ints_nested <- mri_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- mri_melt_top %>%
  unnest(cols=c())
identical(mri_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~MRI_Oedema, data = .x)))
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

best_augmented <- top_200_peaks %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p
best_augmented$Peak_ID <- as.numeric(best_augmented$Peak_ID)
best_hmdb <- inner_join(best_augmented, peak_ID_HMDB, by='Peak_ID')
best_hmdb <- subset(best_hmdb, best_hmdb$Putative_Metabolite != 'NA')

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite != 'Citrate')

sig_peaks <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
sig_peaks <- distinct(sig_peaks, HMDB, .keep_all = TRUE)
sig_peaks <- sig_peaks[-1,]

best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

mri_oed <-best_adj_hmdb 
colnames(mri_oed)[17] <-'Disease_Measure'
mri_oed$Disease_Measure_Type <- 'MRI_Oedema'

ggplot(best_adj_hmdb,aes(x = MRI_Erosion, y=Peak_Intensity)) +
  geom_point() + 
  stat_cor(method = "spearman", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  labs(x='ΔMRI Erosion',
       y='ΔPeak Intensity')+
  theme_minimal()+
  xlim(-2,5)

mets_reduce <- function(disease_df){
  disease_met <- disease_df[c(28:29)]
  disease_met <- distinct(disease_met, Putative_Metabolite, .keep_all = TRUE)
}

mri_oed_red <- mets_reduce(mri_oed)
metabolite_counts_mri_oed <- as.data.frame(as.matrix(table(mri_oed_red$Putative_Metabolite, mri_oed_red$Disease_Measure_Type)))
names(metabolite_counts_mri_oed)<- c('Putative_Metabolite', 'Disease_Measure','Frequency')

write.csv(mri_oed_red, '20210118_AF_diff_mri_oed_mets_count.csv')
### Directly investigating differential abundance of metabolites of interest
mri_melt_top$Oedema_Response <- 0
mri_melt_top$Oedema_Response[mri_melt_top$MRI_Oedema >=-2] <- 'Negative'
mri_melt_top$Oedema_Response[mri_melt_top$MRI_Oedema < -2] <- 'Positive'

stat_test <- mri_melt_top %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Oedema_Response) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

write.csv(stat_test, '20210113_AF_diff_oedema_stats.csv')

stat_test%>%
  mutate(Colour = ifelse(p < 0.05, "P.Value < 0.05", "P.Value > 0.05")) %>%
  ggplot(aes(p))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  labs(x='p-value',
       y='Frequency')
