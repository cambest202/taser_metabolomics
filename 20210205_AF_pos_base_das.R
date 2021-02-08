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
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)

names(peak_metadata)[6] <- 'Peak_ID'

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
resp_ints_melt <- resp_ints[,c(2,13:1596)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_ints_melt$DAS44 <- resp_ints$DAS44
resp_ints_melt$Peak_ID <- gsub('X','', resp_ints_melt$Peak_ID)
resp_ints_melt$CRP <- resp_ints$CRP
resp_ints_melt$ESR <- resp_ints$ESR
resp_ints_melt$HAQ <- resp_ints$HAQ
resp_ints_melt$GHVAS <- resp_ints$GHVAS
resp_ints_melt$PVAS <- resp_ints$PVAS
resp_ints_melt$RAI <- resp_ints$RAI
resp_ints_melt$SJC <- resp_ints$SJC

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
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.01) 
sig_peaks <- distinct(best_adj, Peak_ID, .keep_all = TRUE)

peak_metadata_dist <- distinct(peak_metadata, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb <- inner_join(best_adj, peak_metadata_dist, by='Peak_ID')
#best_adj_hmdb <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all = TRUE)
best_adj_hmdb<- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')

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
       y='Peak Intensity')+
  theme_minimal()

## Determine the most important features to predict clinical outcomes
# Reduce number of features through those correlating peaks
resp_split <- resp_PN_ints
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
index <- createDataPartition(resp_red_t$Response, p = 0.85, list = FALSE)
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
                                                  number =10,
                                                  repeats = 10, 
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
index <- createDataPartition(resp_red_t$Response, p = 0.6, list = FALSE)
train_data <- resp_red_t[index, ]
test_data  <- resp_red_t[-index, ]
model <- glm(Response ~.,family=binomial(link='logit'),data=train_data)
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

## refining features in model
