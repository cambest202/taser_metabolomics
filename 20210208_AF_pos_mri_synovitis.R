## AF Differential metabolites and association with MRI Synovitis after 18 months 
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
library(ROCR) 

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
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
mri_erosion <- read.csv('20201202_MRI_mean.csv', header=TRUE, na.strings=c("","NA"))
xray_erosion <- read.csv('20201202_xray_mean.csv', header=TRUE, na.strings=c("","NA"))

names(mri_erosion)[1] <- 'Sample_Name'
### Differential Analysis: Metabolic changes between positive and negative responders for XRay responses ----------
# Associations of Δmetabolite levels with the erosion response
resp_diff_mri <- inner_join(mri_erosion, resp_diff, by='Sample_Name')
resp_diff_mri$Synovitis_Response <- 0
resp_diff_mri <- resp_diff_mri[c(1:4,1601,6:1599)]
summary(resp_diff_mri$ΔMRI_Synovitis)
resp_diff_mri$Synovitis_Response[resp_diff_mri$ΔMRI_Synovitis >= -6.5] <-'Negative'
resp_diff_mri$Synovitis_Response[resp_diff_mri$ΔMRI_Synovitis < -6.5] <-'Positive'

ggplot(resp_diff_mri)+
  geom_bar(aes(x=Synovitis_Response,
               fill=Synovitis_Response))+
  theme_light()+
  theme(legend.position='none')

resp_mri_matrix <- resp_diff_mri[c(5,16:1584)]
rownames(resp_mri_matrix) <- paste0(rownames(resp_mri_matrix), resp_mri_matrix$Synovitis_Response)
resp_mri_matrix <- resp_mri_matrix[,-1]
resp_mri_t <- as.data.frame(t(resp_mri_matrix))
names(resp_mri_t)[names(resp_mri_t) %like% 'Positive'] <- 'Positive'
names(resp_mri_t)[names(resp_mri_t) %like% 'Negative'] <- 'Negative'

#Limma----
Group <- factor(colnames(resp_mri_t), levels = c('Positive', 'Negative'))
design <- model.matrix (~Group)
colnames(design) <- c('Positive', 'PositivevsNegative')
eset <- resp_mri_t
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
toptable$Sig_Peaks <- ifelse(toptable$P.Value<0.05 & toptable$Putative_Metabolite != '', toptable$Putative_Metabolite, '')
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

## P-value histogram---
# Adding local FDR, qvalue and π0 values to p-value histogram
pi0 <- 2*mean(toptable$P.Value > 0.05)
lfdrvals <- lfdr(toptable$P.Value, pi0)
qobj <- qvalue(toptable$P.Value)
hist(qobj)

### Differential Analysis: Metabolic changes between baseline and 18 months ----------
resp_ints_mri <- inner_join(mri_erosion, resp_ints, by='Sample_Name')
resp_ints_mri <- subset(resp_ints_mri, resp_ints_mri$time == 'A')
resp_ints_mri$Synovitis_Response <- 0
resp_ints_mri <- resp_ints_mri[c(1:4,1600,6:1599)]
summary(resp_ints_mri$ΔMRI_Synovitis)
resp_ints_mri$Synovitis_Response[resp_ints_mri$ΔMRI_Synovitis >= -6.5] <-'Negative'
resp_ints_mri$Synovitis_Response[resp_ints_mri$ΔMRI_Synovitis < -6.5] <-'Positive'
resp_base_mri_matrix <- resp_ints_mri[c(5,16:1599)]
rownames(resp_base_mri_matrix) <- paste0(rownames(resp_base_mri_matrix), resp_base_mri_matrix$Synovitis_Response)
resp_base_mri_matrix <- resp_base_mri_matrix[,-1]
resp_base_mri_t <- as.data.frame(t(resp_base_mri_matrix))
names(resp_base_mri_t)[names(resp_base_mri_t) %like% 'Positive'] <- 'Positive'
names(resp_base_mri_t)[names(resp_base_mri_t) %like% 'Negative'] <- 'Negative'

#Limma----
Group_base <- factor(colnames(resp_base_mri_t), levels = c('Positive', 'Negative'))
design_base <- model.matrix (~Group_base)
colnames(design_base) <- c('Positive', 'PositivevsNegative')
eset_base <- resp_base_mri_t
fit_base <- lmFit(eset_base, design_base)
fit_base <- eBayes(fit_base)
toptable_base <- topTable(fit_base, coef = 'PositivevsNegative', adjust = 'BH', number = 1584)
toptable_base <- as.data.frame(toptable_base)
toptable_base$Peak_ID <- rownames(toptable_base)
toptable_base$Peak_ID <- gsub('X', '', toptable_base$Peak_ID)
toptable_base <- toptable_base[,c(ncol(toptable_base),1:(ncol(toptable_base)-1))]
toptable_base <- join(toptable_base, peak_metadata, by = 'Peak_ID')
##-----
toptable_base$Sig <- 0
toptable_base$Sig <- ifelse(toptable_base$adj.P.Val <0.05, 'Significant', 'Not significant') 
toptable_base$Sig_Peaks <- ifelse(toptable_base$P.Value<0.05 & toptable_base$Putative_Metabolite != '', toptable_base$Putative_Metabolite, '')
toptable_base_dist <- distinct(toptable_base, Putative_Metabolite, .keep_all = TRUE)

ggplot(data=toptable_base_dist, aes(x=logFC, y=-log10(P.Value), 
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

## P-value histogram---
# Adding local FDR, qvalue and π0 values to p-value histogram
pi0 <- 2*mean(toptable_base$P.Value > 0.05)
lfdrvals <- lfdr(toptable_base$P.Value, pi0)
qobj <- qvalue(toptable_base$P.Value)
hist(qobj)

# PCA of samples
scaled_intensities <- scale(t(resp_base_mri_t))
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


### Correlations of the differential metabolites with the DAS44
resp_diff_melt <- resp_diff_mri[,c(1,16:1599)]
resp_diff_melt <- melt(resp_diff_melt)
names(resp_diff_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_diff_melt$DAS44 <- resp_diff_mri$DAS44
resp_diff_melt$Peak_ID <- gsub('X','', resp_diff_melt$Peak_ID)
resp_diff_melt$CRP <- resp_diff_mri$CRP
resp_diff_melt$ESR <- resp_diff_mri$ESR
resp_diff_melt$HAQ <- resp_diff_mri$HAQ
resp_diff_melt$GHVAS <- resp_diff_mri$GHVAS
resp_diff_melt$PVAS <- resp_diff_mri$PVAS
resp_diff_melt$RAI <- resp_diff_mri$RAI
resp_diff_melt$SJC <- resp_diff_mri$SJC
resp_diff_melt$MRI_Erosion <- resp_diff_mri$ΔMRI_Erosion
resp_diff_melt$MRI_Synovitis <- resp_diff_mri$ΔMRI_Synovitis
resp_diff_melt$MRI_Oedema <- resp_diff_mri$ΔMRI_Oedema

# MRI Synovitis----
# make sure all disease measures are included in the melted df
ints_nested <- resp_diff_melt %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_diff_melt %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

best_augmented <- bestest_fit %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)
sig_peaks$ID <- sig_peaks$Peak_ID
best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Peak_ID != '168')
ggplot(best_adj_hmdb,aes(x = MRI_Synovitis, y=Peak_Intensity)) +
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
  
  labs(x='ΔMRI Synovitis',
       y='ΔPeak Intensity')+
  theme_minimal()


# MRI Erosion (baseline correlations)----
resp_ints_melt <- resp_ints_mri[,c(1,16:1599)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_ints_melt$DAS44 <- resp_ints_mri$DAS44
resp_ints_melt$Peak_ID <- gsub('X','', resp_ints_melt$Peak_ID)
resp_ints_melt$CRP <- resp_ints_mri$CRP
resp_ints_melt$ESR <- resp_ints_mri$ESR
resp_ints_melt$HAQ <- resp_ints_mri$HAQ
resp_ints_melt$GHVAS <- resp_ints_mri$GHVAS
resp_ints_melt$PVAS <- resp_ints_mri$PVAS
resp_ints_melt$RAI <- resp_ints_mri$RAI
resp_ints_melt$SJC <- resp_ints_mri$SJC
resp_ints_melt$MRI_Erosion <- resp_ints_mri$ΔMRI_Erosion
resp_ints_melt$MRI_Synovitis <- resp_ints_mri$ΔMRI_Synovitis
resp_ints_melt$MRI_Oedema <- resp_ints_mri$ΔMRI_Oedema

ints_nested <- resp_ints_melt %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_ints_melt %>%
  unnest(cols=c())
identical(resp_ints_melt, ints_unnested)
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

best_augmented <- bestest_fit %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)

best_adj_hmdb <- inner_join(best_adj, peak_metadata, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
best_adj_hmdb <- subset(best_adj_hmdb,best_adj_hmdb$Putative_Metabolite!= 'Pyroglutamic acid')
sig_dist <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)
sig_dist$ID <- sig_dist$Peak_ID

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
  stat_cor(method = "pearson", 
           vjust=1, hjust=0,
           size=3)+
  labs(x='ΔMRI Synovitis',
       y='Peak Intensity',
       title='Positive Ion Mode')+
  theme_minimal()

### Machine learning- feature selection
toptable_sig <- subset(toptable, toptable$Sig_Peaks != '')

resp_mri <- resp_ints_mri[,c(5,16:1599)]
resp_mri_lr <- resp_mri[,-1]
resp_mri_lr <- as.data.frame(t(resp_mri_lr))
resp_mri_lr$Peak_ID <- rownames(resp_mri_lr)
resp_mri_lr$Peak_ID <- gsub('X', '', resp_mri_lr$Peak_ID)
#resp_mri_sel <- subset(resp_mri_lr, resp_mri_lr$Peak_ID %in% toptable_sig$Peak_ID)
resp_mri_matrix <- resp_mri_lr[,-40]
resp_mri_matrix <- as.data.frame(t(resp_mri_matrix))
resp_mri_matrix$Response <- resp_ints_mri$Synovitis_Response
resp_mri_matrix_2 <- resp_mri_matrix
resp_mri_matrix$Response[resp_mri_matrix$Response %like% 'Positive'] <- 1
resp_mri_matrix$Response[resp_mri_matrix$Response %like% 'Negative'] <- 0
resp_mri_matrix$Response <- as.numeric(resp_mri_matrix$Response)

write.csv(resp_mri_matrix, '20210218_AF_pos_syno.csv')


set.seed(42)
index <- createDataPartition(resp_mri_matrix$Response, p = 0.75, list = FALSE) # 0.8
train_data <- resp_mri_matrix[index, ]
test_data  <- resp_mri_matrix[-index, ]

model <- glm(Response ~ X396 + X873 + X359,
             family=binomial(link='logit'),
             data=train_data)
summary(model)
anova(model, test="Chisq")

1-pchisq(43.860-35.428, 31-26)

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
auc
# Log regression not successful in picking out metabolites

## Attempt rf feature selection
set.seed(42)
index <- createDataPartition(resp_mri_matrix_2$Response, p = 0.85, list = FALSE) # 0.8
train_data <- resp_mri_matrix_2[index, ]
test_data  <- resp_mri_matrix_2[-index, ]
tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
set.seed(42)

model_rf <- train(Response~.,
                          data = train_data,
                          method = "rf",
                          metric = "Accuracy",
                          tuneGrid=tuneGrid,
                          trControl = trainControl(method = "repeatedcv",
                                                   number =3,
                                                   repeats = 3, 
                                                   savePredictions = TRUE, 
                                                   verboseIter = FALSE, 
                                                   allowParallel = TRUE),
                          importance = TRUE,
                          ntree = 500)

plot(model_rf)
test_results <- predict(model_rf, newdata = test_data)
summary(test_results)
test_data$Response <- as.factor(test_data$Response)
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
imps <- subset(imps,imps$Importance >50)
imps$Peak_ID <- gsub('X','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_metadata, by='Peak_ID')
imps_hmdb <- with(imps_hmdb,imps_hmdb[order(Importance),])
imps_hmdb_id <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='NA')
imps_hmdb_id <- distinct(imps_hmdb_id, Putative_Metabolite, .keep_all = TRUE)
imps_hmdb_id[7,9] <- 'Pyroglutamate'
imps_hmdb_id$Putative_Metabolite <- as.factor(imps_hmdb_id$Putative_Metabolite)

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
  scale_colour_brewer(palette = "Set1")+
  theme(axis.text.y = element_text(size = 8))


