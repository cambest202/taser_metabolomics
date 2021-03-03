## AF baseline analysis and association of metabolites with DAS44. Standardising the approach

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


`%notin%` <- Negate(`%in%`)

setwd('/Users/cameronbest/University/MRC_DTP_in_Precision_Medicine/Project/RA/Taser_data/Analysis/taser_metabolomics')

### files------
resp_ints <- read.csv('20201105_resp_ints.csv', header=TRUE)
resp_diff <- read.csv('20201105_resp_diff.csv', header=TRUE)
peak_IDs <- read.table (file="20200117_Taser_PeakIDs.csv", header=TRUE, row.names=1, sep= "\t")
peak_ID_HMDB <- read.csv(file='peak_IDs_HMDB.csv', header=TRUE, row.names=1)
peak_metadata <- read.csv(file='20210217_neg_peak_metadata.csv', header=TRUE)

sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)

names(sample_sheet)[2] <- 'Sample_Name'
peak_IDs$Peak_ID <- rownames(peak_IDs)
peak_IDs$Peak_ID  <- as.numeric(peak_IDs$Peak_ID )
peak_IDs$moleculeName  <- as.character(peak_IDs$moleculeName )
peak_metadata <- left_join(peak_metadata, peak_IDs, by='Peak_ID')
names(sample_sheet)[2] <- 'Sample_Name'

peak_metadata_IDs <- subset(peak_metadata, peak_metadata$identification != '')
peak_metadata$Peak_ID_test <- peak_metadata$Peak_ID

# Which metabolites at baseline levels are associated with the DAS44 clinical outcome after 18 months of treatment?
resp_ints_A <- subset(resp_ints,resp_ints$Sample == 'A')
resp_A_melt <- resp_ints_A[,c(2,9:1466)]
resp_A_melt <- melt(resp_A_melt)
names(resp_A_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_A_melt$DAS44 <- resp_ints_A$DAS44
resp_A_melt$CRP <- resp_diff$CRP
resp_A_melt$ESR <- resp_diff$ESR
resp_A_melt$HAQ <- resp_diff$HAQ
resp_A_melt$GHVAS <- resp_diff$GHVAS
resp_A_melt$PVAS <- resp_diff$PVAS
resp_A_melt$RAI <- resp_diff$RAI
resp_A_melt$SJC <- resp_diff$SJC

resp_A_melt$Peak_ID <- gsub('X','', resp_A_melt$Peak_ID)

resp_differential <- resp_ints_A[c(8:1466)]
resp_differential_2 <- resp_differential
resp_differential$Response <- as.factor(resp_differential$Response)
rownames(resp_differential) <- paste0(rownames(resp_differential), resp_differential$Response)
resp_differential <- resp_differential[,-1]
resp_differential_t <- t(resp_differential)
rownames(resp_differential_t) <- gsub('X','', rownames(resp_differential_t))
colnames(resp_differential_t)[colnames(resp_differential_t) %like% 'Positive'] <- 'Positive'
colnames(resp_differential_t)[colnames(resp_differential_t) %like% 'Negative'] <- 'Negative'
resp_differential_3 <- resp_differential_2[,-1]
resp_select_peaks <- t(resp_differential_3)
rownames(resp_select_peaks) <- gsub('X','', rownames(resp_select_peaks))
resp_subset <- t(resp_select_peaks)
resp_subset <- as.data.frame(resp_subset)
resp_subset$Response <- resp_differential_2$Response
resp_subset <- resp_subset[,c(ncol(resp_subset),1:(ncol(resp_subset)-1))]
resp_subset$Response <- as.factor(resp_subset$Response)
colnames(resp_subset)[2:1459] <- paste0('X', colnames(resp_subset)[2:1459])

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

AF_limma <- limma_fun(resp_differential_t, 1500, 'Negative', 'Positive')
AF_limma$Sig <- 0
AF_limma$Sig <- ifelse(AF_limma$adj.P.Val <0.05, 'Significant', 'Not Significant') 
AF_limma$Sig_Peaks <- ifelse(AF_limma$P.Value<0.05 & AF_limma$identification != '', AF_limma$Peak_ID, '')

AF_limma_hmdb <- inner_join(AF_limma, peak_ID_HMDB, by='Peak_ID')
AF_limma_hmdb$Sig_Peaks <- ifelse(AF_limma_hmdb$P.Value<0.05, AF_limma_hmdb$Putative_Metabolite, '')
AF_limma_hmdb <- distinct(AF_limma_hmdb, Peak_ID, .keep_all = TRUE)
AF_limma_hmdb <- AF_limma_hmdb[-8,]
AF_limma_hmdb$Sig <- as.factor(AF_limma_hmdb$Sig)

AF_limma_hmdb_sig <- subset(AF_limma_hmdb, AF_limma_hmdb$P.Value < 0.05)

ggplot(data=AF_limma_hmdb, aes(x=logFC, y=-log10(P.Value), 
                               colour=Sig, 
                               group=Sig)) +
  geom_point (alpha=0.7) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Signficance',
        title='Baseline Putative Metabolites \nDifferential Abundance Across DAS44 Response Subgroups') +
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                   box.padding = 0.5,
                   position = position_jitter(seed = 1),
                   arrow = arrow(length = unit(0.0015, "npc"))) +  
  scale_color_brewer(palette = "Set1")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

base_diff_das <- subset(AF_limma_hmdb,AF_limma_hmdb$P.Value < 0.05)
base_diff_das <- subset(base_diff_das,base_diff_das$Putative_Metabolite !='NA')
base_diff_sim <- base_diff_das[,c(1,3,15)]


# Histogram of p-values 
AF_limma%>%
  mutate(Colour = ifelse(P.Value < 0.05, "p-value < 0.05", "p-value > 0.05")) %>%
ggplot(aes(P.Value))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  labs(x='p-value',
       y='Frequency')

library(qvalue)
pi0 <- 2*mean(AF_limma$P.Value > 0.05)
slfdrvals <- lfdr(AF_limma$P.Value, pi0)
qobj <- qvalue(AF_limma$P.Value)
hist(qobj)

# PCA of samples
scaled_intensities <- scale(t(resp_differential_t))
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
set.seed(42)
index <- createDataPartition(resp_subset$Response, p = 0.7, list = FALSE)
train_data <- resp_subset[index, ]
test_data  <- resp_subset[-index, ]
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
set.seed(42)
model_rf <- caret::train(Response~.,
                         data = train_data,
                         method = "rf",
                         metric = "Accuracy",
                         tuneGrid=tuneGrid,
                         trControl = trainControl(method = "repeatedcv",
                                                  number =5,
                                                  repeats = 3, 
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

final <- data.frame(actual = test_data$Response,
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
imps <- subset(imps,imps$Importance >55)
imps$Peak_ID <- gsub('X','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_ID_HMDB, by='Peak_ID')
imps_hmdb <- with(imps_hmdb,imps_hmdb[order(Importance),])
imps_hmdb_id <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='NA')
imps_hmdb_id <- distinct(imps_hmdb_id, Putative_Metabolite, .keep_all = TRUE)
imps_hmdb_id$Putative_Metabolite <- as.factor(imps_hmdb_id$Putative_Metabolite)

#write.csv(imps_hmdb_id, '20210114_AF_base_das_FI.csv')
#imps_hmdb_id <- read.csv('20210114_AF_base_das_FI.csv')

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
# DAS44 ------
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

top_20_peaks <- head(bestest_fit, 20)
top_50_peaks <- head(bestest_fit, 50)
top_100_peaks <- head(bestest_fit, 100)
top_200_peaks <- head(bestest_fit, 200)

best_augmented <- bestest_fit %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
bad_adj <- subset(best_augmented, best_augmented$p.value > 0.5) 
bad_adj$Peak_ID <- as.numeric(bad_adj$Peak_ID)
bad_peaks <- inner_join(bad_adj, peak_metadata, by='Peak_ID')

best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p
best_augmented$Peak_ID <- as.numeric(best_augmented$Peak_ID)
best_hmdb <- inner_join(best_augmented, peak_ID_HMDB, by='Peak_ID')
best_hmdb <- subset(best_hmdb, best_hmdb$Putative_Metabolite != 'NA')

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 

no_ID <- subset(best_augmented_sig, best_augmented_sig$p.value < 0.05) 

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite != 'Citrate')

sig_peaks <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
sig_peaks <- distinct(sig_peaks, HMDB, .keep_all = TRUE)
tab <- sig_peaks[,c(28,29,30)]
tab_flex <- flextable::flextable(tab)

best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)
batches <- sample_sheet[c(2,11)]
batches <- distinct(batches,Sample_Name, .keep_all = TRUE)
batches$Batch <- as.factor(batches$Batch)
best_adj_sample <- inner_join(best_adj_hmdb, batches, by='Sample_Name')
best_adj_sample<- subset(best_adj_sample,best_adj_sample$Putative_Metabolite != 'Fumarate?')
best_adj_sample<- subset(best_adj_sample,best_adj_sample$Putative_Metabolite != 'Traumatic acid')


mets <- c(166:194)

best_adj_hmdb %>%
  #subset(Peak_ID %notin% mets) %>%
  subset(Putative_Metabolite != 'Fumarate?')%>%
ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(size=1,
             #aes(colour= Batch),
             #colour= ifelse(best_adj_hmdb$Peak_Intensity > 17.5, 'red', 'black'),
             alpha=0.7) + 
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
    theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  stat_cor(method = "spearman", 
           vjust=1, hjust=-1.7,
           size=4)+
  labs(x='ΔDAS44',
       y='Peak Intensity',
       title='Negative Ion Mode')+
  theme_minimal()+
  theme(legend.position = 'none')

resp_melt_samples <- subset(resp_melt_top, resp_melt_top$Peak_ID %in% best_adj_hmdb$Peak_ID)
best_adj_hmdb$Sample_Name <- resp_melt_samples$Sample_Name
best_adj_hmdb <- best_adj_hmdb[,c(ncol(best_adj_hmdb),1:(ncol(best_adj_hmdb)-1))]
best_adj_hmdb$Outliers <- 0
best_adj_hmdb$Outliers <- ifelse(best_adj_hmdb$Peak_Intensity > 17 & best_adj_hmdb$Putative_Metabolite != 'Traumatic acid', best_adj_hmdb$Sample_Name, '')
best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID != 351)

DAS_metabolites <- best_adj_hmdb
colnames(DAS_metabolites)[17] <-'Disease_Measure'
DAS_metabolites$Disease_Measure_Type <- 'DAS'

ggplot(best_adj_hmdb,aes(x = DAS44, y=Peak_Intensity)) +
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
  labs(x='ΔDAS44',
       y='ΔPeak Intensity')+
  geom_text_repel(label = best_adj_hmdb$Outliers) + 
  theme_minimal()

# Investigating outliers ------
outlier_samples <- as.data.frame(best_adj_hmdb$Outliers)
names(outlier_samples)[1] <- 'Sample_Name'
outlier_samples <- subset(outlier_samples, outlier_samples !='')
outlier_samples <- distinct(outlier_samples,Outliers, .keep_all = TRUE)

patient_metadata_outliers <- patient_metadata
patient_metadata_outliers$Sample_Name <- rownames(patient_metadata_outliers)
patient_metadata_outliers <- merge(patient_metadata_outliers, outlier_samples)
patient_metadata_outliers$Outlier <- 'Outlier'
patient_metadata_samples <- patient_metadata
patient_metadata_samples$Sample_Name <- rownames(patient_metadata_samples)
patient_metadata_samples$Outlier <- 'Non-Outlier'
patient_metadata_samples<- subset(patient_metadata_samples, patient_metadata_samples$Sample_Name != patient_metadata_outliers$Sample_Name)
patient_metadata_outliers_2 <- rbind.data.frame(patient_metadata_outliers, patient_metadata_samples)

t.test(patient_metadata_outliers_2$Age ~ patient_metadata_outliers_2$Outlier)

patient_metadata_outliers_2%>%
ggplot(aes(y=Age,x=Outlier))+
  geom_boxplot(aes(fill=Outlier))+
  theme_light()+
  labs(x='Age',
       y=NULL)+
  theme(legend.position = 'none')+
  stat_compare_means(aes(x=Outlier, y=Age),
                     method='wilcox')

sex_outliers <- prop.table(table(patient_metadata_outliers_2$Sex, 
                 patient_metadata_outliers_2$Outlier), margin=2)*100

smoker_outliers <- prop.table(table(patient_metadata_outliers_2$Smoking, 
                                 patient_metadata_outliers_2$Outlier), margin=2)*100

table(patient_metadata_outliers_2$Smoking, patient_metadata_outliers_2$Outlier)

prop.test(x=c(7, 32),
          n=c(11,76),
          p = NULL, alternative = "two.sided",
          correct = TRUE)

prop.test(x=c(42.1, 63.6),
          n=c(100,100),
          p = NULL, alternative = "two.sided",
          correct = TRUE)

# CRP------
resp_melt_top <- subset(resp_A_melt, resp_A_melt$Peak_ID %in% imps_hmdb$Peak_ID)

ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~CRP, data = .x)))
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

CRP_metabolites <- best_adj_hmdb
colnames(CRP_metabolites)[17] <-'Disease_Measure'
CRP_metabolites$Disease_Measure_Type <- 'CRP'

sig_peaks <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
sig_peaks <- distinct(sig_peaks, HMDB, .keep_all = TRUE)

best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

ggplot(best_adj_hmdb,aes(x = CRP, y=Peak_Intensity)) +
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
  labs(x='ΔCRP',
       y='ΔPeak Intensity')+
  theme_minimal()

resp_melt_samples <- subset(resp_melt_top, resp_melt_top$Peak_ID %in% best_adj_hmdb$Peak_ID)
best_adj_hmdb$Sample_Name <- resp_melt_samples$Sample_Name
best_adj_hmdb <- best_adj_hmdb[,c(ncol(best_adj_hmdb),1:(ncol(best_adj_hmdb)-1))]
best_adj_hmdb$Outliers <- 0
best_adj_hmdb$Outliers <- ifelse(best_adj_hmdb$Peak_Intensity > 17 & best_adj_hmdb$Putative_Metabolite != 'L-Valine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Theobromine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Thymine'
                                 & best_adj_hmdb$Putative_Metabolite != 'L-1-Pyrroline 3-hydroxy-5-carboxylate', best_adj_hmdb$Sample_Name, '')

ggplot(best_adj_hmdb,aes(x = CRP, y=Peak_Intensity)) +
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
  labs(x='ΔCRP',
       y='Peak Intensity')+
  geom_text_repel(label = best_adj_hmdb$Outliers) + 
  theme_minimal()+
  xlim(-200,20)

outliers_metadata <- flextable(patient_metadata_outliers)%>%
  autofit()%>%
  theme_vanilla()

sample_sheet_outliers <- subset(sample_sheet, sample_sheet$Sample_Name %in% outlier_samples$Sample_Name)
sample_sheet_outliers_A <- subset(sample_sheet_outliers, sample_sheet_outliers$time == 'A')
sample_sheet_outliers_AF <- subset(sample_sheet_outliers, sample_sheet_outliers$time == 'A' | sample_sheet_outliers$time == 'F' )
sample_sheet_diff_outliers <- aggregate(.~Sample_Name, sample_sheet_outliers_AF[2:10], diff, na.rm=TRUE)
names(sample_sheet_diff_outliers)[2:9] <- paste0('Δ', names(sample_sheet_diff_outliers)[2:9])

sample_outliers <- inner_join(sample_sheet_outliers_A, sample_sheet_diff_outliers, by='Sample_Name')
sample_outliers <- sample_outliers[,-1]
outlier_sample_sheet <- flextable(sample_outliers)%>%
  autofit()%>%
  theme_vanilla()

# ESR------
resp_melt_top <- subset(resp_A_melt, resp_A_melt$Peak_ID %in% imps_hmdb$Peak_ID)

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

ESR_metabolites <- best_adj_hmdb
colnames(ESR_metabolites)[17] <-'Disease_Measure'
ESR_metabolites$Disease_Measure_Type <- 'ESR'

sig_peaks <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
sig_peaks <- distinct(sig_peaks, HMDB, .keep_all = TRUE)

best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

ggplot(best_adj_hmdb,aes(x = ESR, y=Peak_Intensity)) +
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
  labs(x='ΔESR',
       y='ΔPeak Intensity')+
  theme_minimal()

resp_melt_samples <- subset(resp_melt_top, resp_melt_top$Peak_ID %in% best_adj_hmdb$Peak_ID)
best_adj_hmdb$Sample_Name <- resp_melt_samples$Sample_Name
best_adj_hmdb <- best_adj_hmdb[,c(ncol(best_adj_hmdb),1:(ncol(best_adj_hmdb)-1))]
best_adj_hmdb$Outliers <- 0
best_adj_hmdb$Outliers <- ifelse(best_adj_hmdb$Peak_Intensity > 17 & best_adj_hmdb$Putative_Metabolite != '_-Ketoisovaleric acid'
                                 & best_adj_hmdb$Putative_Metabolite != 'L-Valine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Thiocysteine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Hypoxanthine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Traumatic acid'
                                 & best_adj_hmdb$Putative_Metabolite != 'Uridine', best_adj_hmdb$Sample_Name, '')

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
       y='Peak Intensity')+
  geom_text_repel(label = best_adj_hmdb$Outliers) + 
  theme_minimal()

# HAQ------
resp_melt_top <- subset(resp_A_melt, resp_A_melt$Peak_ID %in% imps_hmdb$Peak_ID)

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

HAQ_metabolites <- best_adj_hmdb
colnames(HAQ_metabolites)[17] <-'Disease_Measure'
HAQ_metabolites$Disease_Measure_Type <- 'HAQ'

sig_peaks <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
sig_peaks <- distinct(sig_peaks, HMDB, .keep_all = TRUE)

best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

ggplot(best_adj_hmdb,aes(x = HAQ, y=Peak_Intensity)) +
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
  labs(x='ΔHAQ',
       y='Peak Intensity')+
  theme_minimal()

resp_melt_samples <- subset(resp_melt_top, resp_melt_top$Peak_ID %in% best_adj_hmdb$Peak_ID)
best_adj_hmdb$Sample_Name <- resp_melt_samples$Sample_Name
best_adj_hmdb <- best_adj_hmdb[,c(ncol(best_adj_hmdb),1:(ncol(best_adj_hmdb)-1))]
best_adj_hmdb$Outliers <- 0
best_adj_hmdb$Outliers <- ifelse(best_adj_hmdb$Peak_Intensity > 17 & best_adj_hmdb$Putative_Metabolite != '_-Ketoisovaleric acid'
                                 & best_adj_hmdb$Putative_Metabolite != 'L-Valine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Thiocysteine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Hypoxanthine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Traumatic acid'
                                 & best_adj_hmdb$Putative_Metabolite != 'Uridine', best_adj_hmdb$Sample_Name, '')

ggplot(best_adj_hmdb,aes(x = HAQ, y=Peak_Intensity)) +
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
       y='Peak Intensity')+
  geom_text_repel(label = best_adj_hmdb$Outliers) + 
  theme_minimal()

# HAQ------
resp_melt_top <- subset(resp_A_melt, resp_A_melt$Peak_ID %in% imps_hmdb$Peak_ID)

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
sig_peaks <- subset(sig_peaks, sig_peaks$Peak_ID != '115')
best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

ggplot(best_adj_hmdb,aes(x = HAQ, y=Peak_Intensity)) +
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
  labs(x='ΔHAQ',
       y='Peak Intensity')+
  theme_minimal()

resp_melt_samples <- subset(resp_melt_top, resp_melt_top$Peak_ID %in% best_adj_hmdb$Peak_ID)
best_adj_hmdb$Sample_Name <- resp_melt_samples$Sample_Name
best_adj_hmdb <- best_adj_hmdb[,c(ncol(best_adj_hmdb),1:(ncol(best_adj_hmdb)-1))]
best_adj_hmdb$Outliers <- 0
best_adj_hmdb$Outliers <- ifelse(best_adj_hmdb$Peak_Intensity > 17 & best_adj_hmdb$Putative_Metabolite != '_-Ketoisovaleric acid'
                                 & best_adj_hmdb$Putative_Metabolite != 'L-Valine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Thiocysteine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Hypoxanthine'
                                 & best_adj_hmdb$Putative_Metabolite != 'Traumatic acid'
                                 & best_adj_hmdb$Putative_Metabolite != 'Uridine', best_adj_hmdb$Sample_Name, '')

ggplot(best_adj_hmdb,aes(x = HAQ, y=Peak_Intensity)) +
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
       y='Peak Intensity')+
  geom_text_repel(label = best_adj_hmdb$Outliers) + 
  theme_minimal()


### Shared metabolites correlations
mets_reduce <- function(disease_df){
  disease_met <- disease_df[c(28:29)]
  disease_met <- distinct(disease_met, Putative_Metabolite, .keep_all = TRUE)
}

DAS_met <- mets_reduce(DAS_metabolites)
CRP_met <- mets_reduce(CRP_metabolites)
ESR_met <- mets_reduce(ESR_metabolites)
HAQ_met <- mets_reduce(HAQ_metabolites)

names(DAS_met) <- c('Disease_Measure_Type','Putative_Metabolite')
DAS_met$Disease_Measure_Type <- 'DAS44'

shared_mets_reduced <- rbind.data.frame(DAS_met,
                                        CRP_met,
                                        ESR_met,
                                        HAQ_met)

shared_mets_reduced <- read.csv('20210118_AF_base_das.csv')

metabolite_counts <- as.data.frame(as.matrix(table(shared_mets_reduced$Putative_Metabolite, shared_mets_reduced$Disease_Measure_Type)))
names(metabolite_counts)<- c('Putative_Metabolite', 'Disease_Measure','Frequency')

write.csv(shared_mets_reduced, '20210118_AF_base_das.csv')
metabolite_counts%>%
  ggplot()+
  geom_col(aes(reorder(Putative_Metabolite, Frequency),
               Frequency,
               fill=Disease_Measure))+
  theme_pubclean()+
  coord_flip()+
  labs(x='Putative Metabolite',
       y='Frequency')

### Directly investigating differential abundance of metabolites of interest
resp_ints_melt <- resp_ints[c(8:1466)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Response','Peak_ID','Peak_Intensity')

stat_test <- resp_ints_melt %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Response) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


resp_ints_stats <- inner_join(resp_ints_melt, stat_test, by='Peak_ID')
resp_ints_stats$Peak_ID_2 <- resp_ints_stats$Peak_ID
resp_ints_stats$Peak_ID <- gsub('X', '', resp_ints_stats$Peak_ID)
resp_ints_stats$Peak_ID <- as.numeric(resp_ints_stats$Peak_ID)

resp_ints_stats_hmdb <- inner_join(resp_ints_stats, peak_ID_HMDB, by='Peak_ID')
resp_ints_stats_hmdb_id <- subset(resp_ints_stats_hmdb, resp_ints_stats_hmdb$Putative_Metabolite !='NA')  

resp_ints_stats_hmdb_id$Response <- as.factor(resp_ints_stats_hmdb_id$Response)
resp_ints_stats_hmdb_id$Peak_Intensity <- as.numeric(resp_ints_stats_hmdb_id$Peak_Intensity)
resp_ints_stats_hmdb_id$p.adj <- signif(resp_ints_stats_hmdb_id$p.adj,3)
ab <-resp_ints_stats_hmdb_id %>%
  subset(resp_ints_stats_hmdb_id$p.adj < 0.05)
ab <- as.data.frame(ab[c(2,17)])
ab <- distinct(ab, Putative_Metabolite, .keep_all = TRUE)


diff_analysis <- function(peak_ID, position){
  resp_ints_stats_hmdb_id %>%
    subset(Peak_ID == `peak_ID`)%>%
    ggplot(aes(Response, Peak_Intensity,
               fill=Response))+
    theme_classic2()+
    geom_violin()+
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    stat_pvalue_manual(
      resp_ints_stats_hmdb_id[resp_ints_stats_hmdb_id$Peak_ID == `peak_ID`,], 
      y.position = position,
      label = "p.adj")+
    theme(legend.position = 'none')+
    labs(y='Peak Intensity')
}

diff_analysis(196, 25.2)


### Dealing with poor quality peaks
poor_peaks <- as.numeric(c(85, 642, 925))

AF_limma_hmdb_clean <- AF_limma_hmdb[-c(1,4,9),]

AF_limma_hmdb_clean%>%
  ggplot(aes(x=logFC, y=-log10(P.Value), 
             colour=Sig, 
             group=Sig)) +
  geom_point (alpha=0.7) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Signficance',
        title='Differential Putative Metabolites \nBetween Baseline and 18-Months Post Treatment Initiation') +
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                  box.padding =1,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.0015, "npc"))) +  
  theme(axis.text = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
  scale_color_brewer(palette = "Set1")+
  xlim(-2,5)+
  ylim(0,25)

# Histogram of p-values 
AF_limma_hmdb_clean%>%
  mutate(Colour = ifelse(adj.P.Val < 0.05, "adj.P.Val < 0.05", "adj.P.Val > 0.05")) %>%
  ggplot(aes(P.Value))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  scale_fill_brewer(palette = "Set1")+
  labs(x='p-value',
       y='Frequency')


diff_das_peaks <- subset(AF_limma_hmdb_clean, AF_limma_hmdb_clean$adj.P.Val < 0.05 & AF_limma_hmdb_clean$Putative_Metabolite != 'NA')
diff_das_peaks_sim <- diff_das_peaks[,c(1,4,15)]

### Attempt to differentiate patient responses after 18 months using the baseline levels of selected metabolites
note_mets <- c('Testosterone sulfate', 'Tartaric acid', 'Fumarate' ,'Traumatic acid',
               'L-Histidine', 'Thiocysteine', 'L-Arginine', 'L-Serine', '4-Hydroxybenzoic acid',
               'L-Glyceric acid', 'Uric acid', 'L-Valine', 'Uridine', 'Hypoxanthine',
               "Sphingosine 1-phosphate", "3-hydroxybutyrate", "L-Glutamate", 'L-Tyrosine', "Ribulose-5-phosphate")

note_mets_strict <-c('Traumatic acid', 'Thiocysteine', 'L-Glyceric acid', 'L-Valine', 'Uridine', 'Hypoxanthine',
                     "Sphingosine 1-phosphate", "3-hydroxybutyrate", "L-Glutamate", 'L-Tyrosine', "Ribulose-5-phosphate")

peaks_id_mets <- subset(peak_ID_HMDB, peak_ID_HMDB$Putative_Metabolite %in% note_mets)
peaks_id_strict <- subset(peak_ID_HMDB, peak_ID_HMDB$Putative_Metabolite %in% note_mets_strict)

resp_int <- resp_ints_A[,c(8:1466)]
resp_int_t <- resp_int[,-1]
resp_int_t <- as.data.frame(t(resp_int_t))
resp_int_t$Peak_ID <- rownames(resp_int_t)
resp_int_t$Peak_ID <- gsub('X', '', resp_int_t$Peak_ID )
resp_int_t$Peak_ID <- as.numeric(resp_int_t$Peak_ID)

resp_int_mets <- subset(resp_int_t, resp_int_t$Peak_ID %in% peaks_id_mets$Peak_ID)
resp_int_mets <- resp_int_mets[,-64]
resp_int_mets_t <- as.data.frame(t(resp_int_mets))
resp_int_mets_t$Response <- resp_int$Response
resp_int_mets_t <- resp_int_mets_t[,c(ncol(resp_int_mets_t),1:(ncol(resp_int_mets_t)-1))]

resp_int_strict <- subset(resp_int_t, resp_int_t$Peak_ID %in% peaks_id_strict$Peak_ID)
resp_int_strict <- resp_int_strict[,-64]
resp_int_strict_t <- as.data.frame(t(resp_int_strict))
resp_int_strict_t$Response <- resp_int$Response
resp_int_strict_t <- resp_int_strict_t[,c(ncol(resp_int_strict_t),1:(ncol(resp_int_strict_t)-1))]

# Melt
resp_int_mets_melt <- melt(resp_int_mets_t)
names(resp_int_mets_melt) <- c('Response', 'Peak_ID', 'Peak_Intensity')
resp_int_strict_melt <- melt(resp_int_strict_t)
names(resp_int_strict_melt) <- c('Response', 'Peak_ID', 'Peak_Intensity')

resp_int_strict_melt %>%
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

#### Using logistic regression for prediction of patient outcomes. 
resp_int_2 <- resp_int
resp_int_2$Response[resp_int_2$Response=='Positive'] <- 1
resp_int_2$Response[resp_int_2$Response=='Negative'] <- 0
resp_int_2$Response <- as.factor(resp_int_2$Response)

#write.csv(resp_int_2, '20210218_AF_neg_LGR_matrix') 

set.seed(42)
index <- createDataPartition(resp_int_2$Response, p = 0.65, list = FALSE)
train_data <- resp_int_2[index, ]
test_data  <- resp_int_2[-index, ]

model <- glm(Response ~ X243 + X189 + X182 + X172 + X351 + X1028,
             family=binomial(link='logit'),data=train_data)
summary(model)
anova(model, test="Chisq")

#library(pscl)
pR2(model)

fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_data$Response)
print(paste('Accuracy',1-misClasificError))

#library(ROCR)
p <- predict(model,newdata=test_data,type='response')
pr <- prediction(p, test_data$Response)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

# Repeat using fewer variables. Do this to avoid overfitting, selecting only the most significant
model <- glm(Response ~ X115 + X270 + X642 + X1191 + X1361,
             family=binomial(link='logit'),
             data=train_data)
summary(model)
confint(model)
anova(model, test="Chisq") # check the overall effect of the variables on the dependent variable
pR2(model)
fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_data$Response)
print(paste('Accuracy',1-misClasificError))

p <- predict(model,newdata=test_data,type='response')
pr <- prediction(p, test_data$Response)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc # higher the score, the better the model

