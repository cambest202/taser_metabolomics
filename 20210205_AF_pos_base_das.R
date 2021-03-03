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
peak_metadata <- read.csv(file='20200521_Taser_POS_Peakdata.csv', header=TRUE)

sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)

peak_metadata <- peak_metadata[,-c(1,2)]
names(peak_metadata)[6] <- 'Peak_ID'

# fix disease measures in resp_ints
resp_diff_double <- rbind.data.frame(resp_diff,resp_diff)
resp_diff_double <- with(resp_diff_double,resp_diff_double[order(Sample_Name),])
resp_ints$CRP <- resp_diff_double$CRP
resp_ints$DAS44 <- resp_diff_double$DAS44
resp_ints$ESR <- resp_diff_double$ESR
resp_ints$GHVAS <- resp_diff_double$GHVAS
resp_ints$HAQ <- resp_diff_double$HAQ
resp_ints$PVAS <- resp_diff_double$PVAS
resp_ints$RAI <- resp_diff_double$RAI
resp_ints$SJC <- resp_diff_double$SJC
resp_ints$DAS44_Response[resp_ints$DAS44 < -2.4] <- 'Positive'
resp_ints$DAS44_Response[resp_ints$DAS44 >= -2.4] <- 'Negative'

ggplot(resp_ints)+
  geom_bar(aes(DAS44_Response))

### Differential Analysis: Metabolic changes between positive and negative responders ----------
resp_PN_ints <- subset(resp_ints, resp_ints$time =='A')
resp_PN_ints <- resp_PN_ints[12:1596]

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
toptable$Sig_Peaks <- ifelse(toptable$P.Value<0.05 & toptable$Putative_Metabolite != '', toptable$Putative_Metabolite, '')
toptable_dist <- distinct(toptable, Putative_Metabolite, .keep_all = TRUE)
toptable_AF_sig <- subset(toptable_dist,toptable_dist$P.Value< 0.05)
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

# PCA of samples
scaled_intensities <- scale(t(resp_PN_limma_t))
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
resp_ints_melt <- resp_ints[,c(2,13:1596)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_ints_melt$DAS44 <- resp_diff$DAS44
resp_ints_melt$Peak_ID <- gsub('X','', resp_ints_melt$Peak_ID)
resp_ints_melt$CRP <- resp_diff$CRP
resp_ints_melt$ESR <- resp_diff$ESR
resp_ints_melt$HAQ <- resp_diff$HAQ
resp_ints_melt$GHVAS <- resp_diff$GHVAS
resp_ints_melt$PVAS <- resp_diff$PVAS
resp_ints_melt$RAI <- resp_diff$RAI
resp_ints_melt$SJC <- resp_diff$SJC

resp_ints %>%
  ggplot(aes(X25, X1037))+
  geom_point(size=0.8, alpha=0.7) + 
  stat_cor(method = "spearman", 
           vjust=1, hjust=0,
           size=5)+
  geom_smooth(method='lm',
              colour='red')+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  labs(x='Tyrosine',
       y='Phenylpyruvic acid')+
  theme_minimal()

# DAS44----
# make sure all disease measures are included in the melted df

ints_nested <- resp_ints_melt %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_ints_melt %>%
  unnest(cols=c())
identical(resp_ints_melt, ints_unnested)
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
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
good_augmented_sig <-  subset(best_augmented, best_augmented$p.value< 0.1)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
good_augmented_sig$Peak_ID <- as.numeric(good_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

sig_peaks <- distinct(best_augmented_sig, Peak_ID, .keep_all = TRUE)
peak_metadata_dist <- distinct(peak_metadata, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb <- inner_join(best_augmented_sig, peak_metadata_dist, by='Peak_ID')
good_p_hmdb <- inner_join(good_augmented_sig, peak_metadata_dist, by='Peak_ID')

best_adj_hmdb_dist <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
good_p_hmdb_dist <- distinct(good_p_hmdb, Putative_Metabolite, .keep_all = TRUE)
good_p_hmdb_dist$ID <- good_p_hmdb_dist$Peak_ID

#write.csv(best_adj_hmdb_dist_2, '20210211_AF_base_das_keggs.csv')

plot_aug <- best_adj_hmdb[,c(16,17,31)]
no_mets <- c(6,60,183,186,341,513, 572, 658, 682, 683,1041,1226, 1375,1549)
`%notin%` <- Negate(`%in%`)

best_adj_hmdb_dist <- subset(best_adj_hmdb_dist,best_adj_hmdb_dist$Peak_ID %notin% no_mets)
best_adj_hmdb_dist_2 <-best_adj_hmdb_dist
best_adj_hmdb_dist_2$Peak_ID <- paste0('X',best_adj_hmdb_dist_2$Peak_ID)

best_adj_hmdb%>%
  subset(Peak_ID %notin% no_mets) %>%
  subset(Peak_ID != '1037') %>%
ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(size=0.5, alpha=0.7) + 
  stat_cor(method = "spearman", 
           vjust=1, hjust=0,
           size=3)+
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  labs(x='ΔDAS44',
       y='Peak Intensity',
       title='Positive Ion Mode')+
  theme_minimal()



## Determine the most important features to predict clinical outcomes
# Reduce number of features through those correlating peaks
resp_split <- resp_PN_ints
resp_split <- resp_split[,-1]
resp_split_t <- as.data.frame(t(resp_split))
resp_split_t$Peak_ID <- rownames(resp_split_t)
resp_split_t$Peak_ID <- gsub('X','', resp_split_t$Peak_ID)
#resp_split_done <- subset(resp_split_t, resp_split_t$Peak_ID %in% toptable_AF_sig$Peak_ID)
#resp_split_done <- subset(resp_split_done, resp_split_done$Peak_ID %in% toptable_AF_sig$Peak_ID)

resp_red <- resp_split_t[,-64]
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
model_rf <- caret::train(Response~.,
                         data = train_data,
                         method = "rf",
                         metric = "Accuracy",
                         tuneGrid=tuneGrid,
                         trControl = trainControl(method = "repeatedcv",
                                                  number =15,
                                                  repeats = 5, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE, 
                                                  allowParallel = TRUE),
                         importance = TRUE,
                         ntree = 500)

# num=15, rep=5, ntree=500, split=0.8 Acc= 0.667, K=0.23
# num=15, rep=5, ntree=700, split=0.8 Acc= 0.667, K=0.27

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
imps <- subset(imps,imps$Importance >10)# set 30 for initial selection but plot only > 50
imps$Peak_ID <- gsub('X','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_metadata, by='Peak_ID')
imps_hmdb <- distinct(imps_hmdb, Putative_Metabolite, .keep_all = TRUE)
imps_hmdb <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='')
no_mets <- c(1041,683,186,175,513,183,1375,1103)
imps_hmdb_2 <- subset(imps_hmdb,imps_hmdb$Peak_ID %notin% no_mets)

cor_FI_mets <- subset(best_adj_hmdb_dist_2, best_adj_hmdb_dist_2$Peak_ID %in% imps_hmdb_2$Peak_ID)
cor_FI_mets$Putative_Metabolite

mets_LR <- paste0('X',good_p_hmdb_dist$Peak_ID)


# Plot the annotated peaks from feature importance
ggplot(imps_hmdb_2)+
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

resp_red <- resp_split_t[,-64]
resp_red_t <- as.data.frame(t(resp_red))
resp_red_t$Response <- resp_PN_ints$DAS44_Response
resp_red_t$Response <- as.character(resp_red_t$Response)
resp_red_t$Response[resp_red_t$Response=='Positive'] <- 1
resp_red_t$Response[resp_red_t$Response=='Negative'] <- 0
resp_red_t$Response <- as.integer(resp_red_t$Response)

sig_dist_cor <- subset(best_adj_hmdb, best_adj_hmdb$adj_p < 0.05)
sig_dist_cor <- distinct(sig_dist_cor, Peak_ID, .keep_all = TRUE)
sig_dist_cor$ID <- sig_dist_cor$Peak_ID

write.csv(resp_red_t, '20210218_AF_pos_LGR_matrix') 

set.seed(42)
index <- createDataPartition(resp_red_t$Response, p = 0.75, list = FALSE)
train_data <- resp_red_t[index, ]
test_data  <- resp_red_t[-index, ]

model <- glm(Response ~  X423 + X25  + X1038 + X1037  +
             X1030 + X274 + X484 ,
             family=binomial(link='logit'),
             data=train_data)

#model_mets <- c('X602', 'X5','X6','X168', 'X558','X303','X1549',  'X328')
#model_mets <- gsub('X','', model_mets)

summary(model)
anova(model, test="Chisq")
pR2(model)
fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_data$Response)
print(paste('Accuracy',1-misClasificError))

1-pchisq(60.571-20.583, 44-36)

library(ROCR) 
p <- predict(model,newdata=test_data,type='response')
pr <- prediction(p, test_data$Response)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc


LR_mets <- c(423,25,1038,1037,1030,274,484) 
good_p_hmdb_dist <- subset(good_p_hmdb_dist,good_p_hmdb_dist$Peak_ID == c(423,25,1038,1037,1030,274,484))
good_p_hmdb_dist <- distinct(good_p_hmdb_dist, Putative_Metabolite, .keep_all = TRUE)
good_p_hmdb_dist$Putative_Metabolite

toptable$Sig_Peaks <- ifelse(toptable$Putative_Metabolite == LR_mets, toptable$Putative_Metabolite, '')
toptable_dist <- distinct(toptable, Putative_Metabolite, .keep_all = TRUE)
toptable_AF_sig <- subset(toptable_dist,toptable_dist$P.Value< 0.05)
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