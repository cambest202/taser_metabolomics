## AF baseline analysis and association of metabolites with erosion imaging scoring outcomes

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
peak_metadata <- read.csv(file='20210217_neg_peak_metadata.csv', header=TRUE)
sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
mri_erosion <- read.csv('20201202_MRI_mean.csv', header=TRUE, na.strings=c("","NA"))
xray_erosion <- read.csv('20201202_xray_mean.csv', header=TRUE, na.strings=c("","NA"))

names(sample_sheet)[2] <- 'Sample_Name'
peak_IDs$Peak_ID <- rownames(peak_IDs)
peak_IDs$Peak_ID  <- as.numeric(peak_IDs$Peak_ID )
peak_IDs$moleculeName  <- as.character(peak_IDs$moleculeName )
names(sample_sheet)[2] <- 'Sample_Name'

peak_metadata_IDs <- subset(peak_metadata, peak_metadata$identification != '')
peak_metadata$Peak_ID_test <- peak_metadata$Peak_ID

peak_metadata_HMDB <- read.csv('20210217_neg_peak_metadata.csv')
# Which metabolites over the 18 months are associated with the erosion imaging score outcomes after 18 months of treatment?
# Xray and MRI should be done independently
# MRI Erosion
names(mri_erosion)[1] <- 'Sample_Name'
resp_ints_A <- subset(resp_ints,resp_ints$Sample == 'A')
mri_resp <- inner_join(mri_erosion, resp_ints_A, by='Sample_Name')
mri_melt <- mri_resp[,c(1,12:1469)]
mri_melt <- melt(mri_melt)
mri_melt$DAS44 <- mri_melt$DAS44
mri_melt$MRI_Synovitis <- mri_resp$ΔMRI_Synovitis
mri_melt$CRP <- mri_resp$CRP
mri_melt$ESR <- mri_resp$ESR
mri_melt$HAQ <- mri_resp$HAQ
mri_melt$GHVAS <- mri_resp$GHVAS
mri_melt$PVAS <- mri_resp$PVAS
mri_melt$RAI <- mri_resp$RAI
mri_melt$SJC <- mri_resp$SJC
names(mri_melt)[1:3] <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')

mri_resp$Syno_Response <- 0
mri_pre_limma <- mri_resp[, c(3,1470,12:1469)]
mean(mri_pre_limma$ΔMRI_Synovitis)
median(mri_pre_limma$ΔMRI_Synovitis)
histogram(mri_pre_limma$ΔMRI_Synovitis)

mri_pre_limma$Syno_Response[mri_pre_limma$ΔMRI_Synovitis >= -6.5] <- 'Negative'
mri_pre_limma$Syno_Response[mri_pre_limma$ΔMRI_Synovitis < -6.5] <- 'Positive'
mri_pre_limma_2 <- mri_pre_limma

mri_pre_limma %>%
  ggplot(aes(length(Syno_Response),
             fill=Syno_Response))+
  geom_histogram()+
  facet_wrap(~Syno_Response)+
  theme(legend.position = 'none')

mri_pre_limma <- mri_pre_limma[,-1]
rownames(mri_pre_limma) <- paste0(rownames(mri_pre_limma), mri_pre_limma$Syno_Response)
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

mri_limma <- limma_fun(mri_pre_limma_t, 1500, 'Positive', 'Negative')
mri_limma$Sig <- 0
mri_limma$Sig <- ifelse(mri_limma$adj.P.Val <0.05, 'Significant', 'Not significant') 
mri_limma$Sig_Peaks <- ifelse(mri_limma$P.Value<0.05 & mri_limma$identification != '', mri_limma$Peak_ID, '')

mri_limma_hmdb <- inner_join(mri_limma, peak_ID_HMDB, by='Peak_ID')
mri_limma_hmdb <- distinct(mri_limma_hmdb, Peak_ID, .keep_all = TRUE)
mri_limma_hmdb$Putative_Metabolite[mri_limma_hmdb$Putative_Metabolite =='Citrate'] <-NA
mri_limma_hmdb$Putative_Metabolite[mri_limma_hmdb$Peak_ID==424] <- 'Citrate'
mri_limma_hmdb$Sig_Peaks <- ifelse(mri_limma_hmdb$P.Value<0.05, mri_limma_hmdb$Putative_Metabolite, '')

ggplot(data=mri_limma_hmdb, aes(x=logFC, y=-log10(P.Value), 
                                colour=Sig, 
                                group=Sig)) +
  geom_point () +
  theme_light()+
  labs (x='LogFC',
        y='-Log p-value') +
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                  box.padding = 0.5,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.0015, "npc"))) +  
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))#+
  
#ylim(0,2)


# Histogram of p-values 
mri_limma_hmdb%>%
  mutate(Colour = ifelse(P.Value < 0.05, "p-value < 0.05", "p-value > 0.05")) %>%
  ggplot(aes(P.Value))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  labs(x='p-value',
       y='Frequency')
ab <- as.vector(mri_limma$adj.P.Val)
write.table(ab,'qval.test.txt')
library(qvalue)
pi0 <- 2*mean(mri_limma$P.Value > 0.05)
lfdrvals <- lfdr(mri_limma$P.Value, pi0)
qobj <- qvalue(mri_limma$P.Value, fdr.level = 0.1, pi0.method='smoother')
summary(qobj)
hist(qobj)
plot(qobj) #shows The estimated π0 versus the tuning parameter λ;
#The q-values versus the p-values
#The number of significant tests versus each q-value cut-off
#The number of expected false positives versus the number of significant tests

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
mri_resping <- mri_pre_limma_2[,-1]

set.seed(42)
index <- createDataPartition(mri_resping$Syno_Response, p = 0.7, list = FALSE)
train_data <- mri_resping[index, ]
test_data  <- mri_resping[-index, ]
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
set.seed(42)
model_rf <- caret::train(Syno_Response~.,
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
test_data$Syno_Response <- as.factor(test_data$Syno_Response)
summary(test_data$Syno_Response)
confusionMatrix(test_results, test_data$Syno_Response)

final <- data.frame(actual = test_data$Syno_Response,
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
imps <- subset(imps,imps$Importance >40)
imps$Peak_ID <- gsub('X','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_ID_HMDB, by='Peak_ID')
imps_hmdb <- with(imps_hmdb,imps_hmdb[order(Importance),])
imps_hmdb_id <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='NA')
imps_hmdb_id <- distinct(imps_hmdb_id, Putative_Metabolite, .keep_all = TRUE)
imps_hmdb_id$Putative_Metabolite <- as.factor(imps_hmdb_id$Putative_Metabolite)

#write.csv(imps_hmdb_id, '20210114_AF_base_mri_syno_FI.csv')
#imps_hmdb_id <- read.csv('20210114_AF_base_mri_syno_FI.csv')
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
mri_melt_top <- mri_melt
#subset(mri_melt, mri_melt$Peak_ID %in% imps_hmdb$Peak_ID)

mri_melt_top%>%
  subset(Peak_ID =='74')%>%
  ggplot(aes(y=Peak_Intensity, x=MRI_Synovitis))+
  geom_point()+
  theme_light()+
  geom_smooth(method='lm')+
  stat_cor(method='spearman')+
  labs(x='ΔSynovitis after 18 months',
       y='Peak intensity at baseline')

ints_nested <- mri_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- mri_melt_top %>%
  unnest(cols=c())
identical(mri_melt_top, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~MRI_Synovitis, data = .x)))
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
best_augmented$Peak_ID <- as.numeric(best_augmented$Peak_ID)
best_hmdb <- inner_join(best_augmented, peak_ID_HMDB, by='Peak_ID')
best_hmdb <- subset(best_hmdb, best_hmdb$Putative_Metabolite != 'NA')

best_adj <- subset(best_augmented_sig, best_augmented_sig$p.value < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite != 'Citrate')

sig_peaks <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
sig_peaks <- distinct(sig_peaks, HMDB, .keep_all = TRUE)
sig_peaks$ID <- sig_peaks$Peak_ID
best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

syno_metabolites <- best_adj_hmdb
colnames(syno_metabolites)[17] <-'Disease_Measure'
syno_metabolites$Disease_Measure_Type <- 'MRI_Synovitis'
mets_reduce <- function(disease_df){
  disease_met <- disease_df[c(28:29)]
  disease_met <- distinct(disease_met, Putative_Metabolite, .keep_all = TRUE)
}

mri_syno <- mets_reduce(syno_metabolites)

metabolite_counts <- as.data.frame(as.matrix(table(mri_syno$Putative_Metabolite, mri_syno$Disease_Measure_Type)))
names(metabolite_counts)<- c('Putative_Metabolite', 'Disease_Measure','Frequency')

#write.csv(mri_syno, '20210118_AF_base_mri_syno.csv')

best_adj_hmdb%>%
  ggplot(aes(x = MRI_Synovitis, y=Peak_Intensity)) +
  geom_point(size=0.5, alpha=0.7) + 
  geom_smooth(method='lm',
              colour='red')+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  stat_cor(method = "spearman", 
           vjust=1, hjust=0,
           size=3)+
  labs(x='ΔMRI Synovitis',
       y='Peak Intensity',
       title='Negative Ion Mode')+
  theme_minimal()


?geom_smooth
### Directly investigating differential abundance of metabolites of interest
mri_melt_top$Syno_Response <- 0
mri_melt_top$Syno_Response[mri_melt_top$MRI_Synovitis >= -6.5] <- 'Negative'
mri_melt_top$Syno_Response[mri_melt_top$MRI_Synovitis < -6.5] <- 'Positive'

stat_test <- mri_melt_top %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Syno_Response) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

write.csv(stat_test, '20210114_AF_base_mri_syno_stats.csv')

stat_test%>%
  mutate(Colour = ifelse(p < 0.05, "P.Value < 0.05", "P.Value > 0.05")) %>%
  ggplot(aes(p))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  labs(x='p-value',
       y='Frequency')

stat_syno <- read.csv('20210114_AF_base_mri_syno_stats.csv')

stat_test_sig <- subset(stat_test, stat_test$p.adj < 0.05)
mri_melt_top$Peak_ID <- as.numeric(mri_melt_top$Peak_ID)
mri_melt_hmdb <- inner_join(mri_melt_top, peak_ID_HMDB, by='Peak_ID')

mri_melt_hmdb %>%
  subset(Peak_ID %in% stat_test_sig$Peak_ID)%>%
  ggplot()+
  theme_classic2()+
  geom_boxplot(aes(x=Syno_Response, 
                   y=Peak_Intensity,
                   fill=Syno_Response))+
  facet_wrap(~Putative_Metabolite, scale='free')+
  theme(legend.position='none')+
  labs(x= 'MRI Synovitis Outcome',
       y= 'Peak Intensity')+
  stat_compare_means(aes(x=Syno_Response, y=Peak_Intensity),
                     method='wilcox')
  
## Attempt to differentiate Mri synovitis responses using the selected metabolites at baseline
mets <- c('Uric acid', "3-hydroxybutyrate", "Aminoadipic acid", "LysoPC", "Ribulose-5-phosphate")
#mets <- c('Uric acid', "L-Tyrosine", "Ribulose-5-phosphate")

mri_limma_select <- subset(mri_limma_hmdb, mri_limma_hmdb$P.Value< 0.05)
mri_limma_select <- subset(mri_limma_select, mri_limma_select$Putative_Metabolite != 'NA')
  
syno_mets <- mri_resp[,c(3,12:1470)]
syno_mets$Syno_Response[syno_mets$ΔMRI_Synovitis < -6.5] <- 'Positive'
syno_mets$Syno_Response[syno_mets$ΔMRI_Synovitis >= -6.5] <- 'Negative'
syno_mets_2 <- syno_mets[,c(2:1460)]

syno_mets_t <- syno_mets_2[,-1459]
syno_mets_t <- as.data.frame(t(syno_mets_t))
syno_mets_t$Peak_ID <- rownames(syno_mets_t)
syno_mets_t$Peak_ID <- gsub('X', '', syno_mets_t$Peak_ID)
syno_mets_t$Peak_ID <- as.numeric(syno_mets_t$Peak_ID)

peaks_met <- subset(peak_ID_HMDB, peak_ID_HMDB$Putative_Metabolite %in% mets)

syno_mets_select <- subset(syno_mets_t, syno_mets_t$Peak_ID %in% peaks_met$Peak_ID)
syno_mets_select <- syno_mets_select[,-40]
syno_sel <- as.data.frame(t(syno_mets_select))
syno_sel$Response <- syno_mets$Syno_Response

syno_melt <- melt(syno_sel)
names(syno_melt) <- c('Response', 'Peak_ID', 'Peak_Intensity')

stat_test <- syno_melt %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Response) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

ggplot(syno_melt, aes(Response, Peak_Intensity,
             fill=Response))+
  theme_light()+
  geom_violin()+
  geom_jitter(size=0.5, alpha=0.5, colour='purple')+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  stat_compare_means(method= 'wilcox.test',
                     label = "p.format",
                     vjust=8, 
                     hjust=-1)+
  theme(legend.position = 'none')+
  labs(y='Peak Intensity')

#### Using logistic regression for prediction of patient outcomes. 
mri_resping <- mri_pre_limma_2[,-1]
syno_sel <- mri_resping
names(syno_sel)[1]<- 'Response'
syno_sel$Response[syno_sel$Response == 'Positive'] <- 1
syno_sel$Response[syno_sel$Response == 'Negative'] <- 0
syno_sel$Response <- as.numeric(syno_sel$Response)

#write.csv(syno_sel, '20210218_AF_neg_syno.csv')

set.seed(42)
index <- createDataPartition(syno_sel$Response, p = 0.7, list = FALSE)
train_data <- syno_sel[index, ]
test_data  <- syno_sel[-index, ]

model <- glm(Response ~.,family=binomial(link='logit'),data=train_data)
summary(model)
anova(model, test="Chisq")

#library(pscl)
pR2(model)

fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.4,1,0)
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
model <- glm(Response ~ X74 + X162 + X196,family=binomial(link='logit'),data=train_data)
summary(model)
confint(model)
anova(model, test="Chisq") # check the overall effect of the variables on the dependent variable
pR2(model)
fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.6,1,0)
misClasificError <- mean(fitted.results != test_data$Response)
print(paste('Accuracy',1-misClasificError))

p <- predict(model,newdata=test_data,type='response')
pr <- prediction(p, test_data$Response)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc # higher the score, the better the model

train_data %>%
  mutate(prob = ifelse(Response == "1", 1, 0)) %>%
  ggplot(aes(X196, prob)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title = "Logistic Regression Model", 
    x = "Metabolite Abundance",
    y = "Probability of Positive MRI Synovitis Response")+
  theme_light()

