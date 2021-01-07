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

remotes::install_github("csdaw/ggprism")
library(ggprism)


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

# Which metabolites over the 18 months are associated with the DAS44 clinical outcome after 18 months of treatment?
resp_diff_melt <- resp_diff[,c(13:1471)]
resp_diff_2 <- resp_diff[,c(3,14:1471)]
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

resp_differential <- resp_ints[c(3,9:1466)]
resp_differential_2 <- resp_differential
resp_differential$Sample <- as.factor(resp_differential$Sample)
rownames(resp_differential) <- paste0(rownames(resp_differential), resp_differential$Sample)
resp_differential <- resp_differential[,-1]
resp_differential_t <- t(resp_differential)
rownames(resp_differential_t) <- gsub('X','', rownames(resp_differential_t))
colnames(resp_differential_t)[colnames(resp_differential_t) %like% 'A'] <- 'A'
colnames(resp_differential_t)[colnames(resp_differential_t) %like% 'F'] <- 'F'

resp_diff_ML <- resp_diff[,c(11,14:1471)]
colnames(resp_diff_ML)[2:1459] <- gsub('X','', colnames(resp_diff_ML)[2:1459])
resp_diff_ML$DAS44_Response[resp_diff_ML$DAS44_Response %like% 'Good'|resp_diff_ML$DAS44_Response %like% 'Remission'] <- 'Positive'
resp_diff_ML$DAS44_Response[resp_diff_ML$DAS44_Response %like% 'Poor'|resp_diff_ML$DAS44_Response %like% 'Mild'] <- 'Negative'
resp_diff_ML$DAS44_Response <- as.factor(resp_diff_ML$DAS44_Response)
resp_diff_peaks <- resp_diff_ML[,-1]
resp_diff_peaks <- as.data.frame(t(resp_diff_peaks))

resp_diff_ML_2 <- t(resp_diff_peaks)
resp_diff_ML_2 <- cbind.data.frame(resp_diff_ML[1], resp_diff_ML_2)

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

AF_limma <- limma_fun(resp_differential_t, 1500, 'A', 'F')
AF_limma$Sig <- 0
AF_limma$Sig <- ifelse(AF_limma$adj.P.Val <0.05, 1, 0) 
AF_limma$Sig_Peaks <- ifelse(AF_limma$adj.P.Val<0.001 & AF_limma$identification != '', AF_limma$Peak_ID, '')

ggplot(data=AF_limma, aes(x=logFC, y=-log10(P.Value), 
                          colour=Sig, 
                          group=Sig)) +
  geom_point () +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value') +
  geom_label_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                   box.padding = 1,
                   position = position_jitter(seed = 1),
                   arrow = arrow(length = unit(0.015, "npc"))) +  
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Histogram of p-values 
AF_limma%>%
  mutate(Colour = ifelse(adj.P.Val < 0.05, "adj.P.Val < 0.05", "adj.P.Val > 0.05")) %>%
  ggplot(aes(P.Value))+
  geom_histogram(aes(fill=Colour),binwidth=0.02, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  labs(x='p-value',
       y='Frequency')

# PCA of samples
scaled_intensities <- scale(t(resp_differential_t))
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

# Generate the model
# Set training and testing data
#colnames(resp_diff_ML_2)[2:1459] <- paste0('X', colnames(resp_diff_ML_2)[2:1459])
set.seed(42)
index <- createDataPartition(resp_diff_ML_2$DAS44_Response, p = 0.7, list = FALSE)
train_data <- resp_diff_ML_2[index, ]
test_data  <- resp_diff_ML_2[-index, ]

### caret model ----
tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
set.seed(42)
model_rff <- caret::train(DAS44_Response~.,
                          data = train_data,
                          method = "rf",
                          metric = "Accuracy",
                          tuneGrid=tuneGrid,
                          trControl = trainControl(method = "repeatedcv",
                                                   number =10,
                                                   repeats = 5, 
                                                   savePredictions = TRUE, 
                                                   verboseIter = FALSE, 
                                                   allowParallel = TRUE),
                          importance = TRUE,
                          ntree = 300)

##------
customRF <- list(type = "Classification",
                 library = "randomForest",
                 loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                  class = rep("numeric", 2),
                                  label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree=param$ntree)
}

#Predict label
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

#Predict prob
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
# train model
control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)

tunegrid <- expand.grid(.mtry=c(1:sqrt(1458)),.ntree=c(300,500,700))

set.seed(42)
custom <- caret::train(DAS44_Response~., data=train_data, 
                       method=customRF, 
                       metric='Accuracy', 
                       tuneGrid=tunegrid, 
                       trControl=control)

summary(custom)

plot(custom)
stopCluster(cores)

test_results <- predict(custom, newdata = test_data)
summary(test_results)
summary(test_data$DAS44_Response)
confusionMatrix(test_results, test_data$DAS44_Response)

final <- data.frame(actual = test_data$DAS44_Response,
                    predict(custom, newdata = test_data, type = "prob"))
final$predicted <- 0
final$predicted[final$Positive > final$Negative] <- 'Positive'
final$predicted[final$Positive < final$Negative] <- 'Negative'
final$correct <- 0
final$correct <- ifelse(final$predict == final$actual, 'Correct', 'Incorrect')
summary(final$predicted == final$actual)

imp <- custom$finalModel$importance
imp <- as.data.frame(imp)
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
imp_peaks <- varImp(model_rff, scale = TRUE)
plot(imp_peaks, top=20)
imps <- as.data.frame(imp_peaks$importance)
imps$Peak_ID <- rownames(imps)
imps <- imps[,-1]
colnames(imps)[1] <- 'Importance'
imps <- subset(imps,imps$Importance >20)
imps$Peak_ID <- gsub('`','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_ID_HMDB, by='Peak_ID')
imps_hmdb <- with(imps_hmdb,imps_hmdb[order(Importance),])
imps_hmdb_id <- subset(imps_hmdb, imps_hmdb$Putative_Metabolite !='NA')
imps_hmdb_id <- distinct(imps_hmdb_id, Putative_Metabolite, .keep_all = TRUE)
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
  theme(axis.text.y = element_text(size = 8))



### Build linear models for the top peaks from feature selection for each of the disease measurements
# DAS44----
resp_melt_top <- subset(resp_diff_melt, resp_diff_melt$Peak_ID %in% peak_ID_HMDB$Peak_ID)
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
  
  best_augmented <- top_200_peaks %>% 
    mutate(augmented = map(model, ~augment(.x))) %>% 
    unnest(augmented)
  best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
  best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
  adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
  best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)

DAS_metabolites <- best_adj_hmdb
colnames(DAS_metabolites)[17] <-'Disease_Measure'
DAS_metabolites$Disease_Measure_Type <- 'DAS44'

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
  theme_minimal()

# CRP----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

best_augmented <- top_50_peaks %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

CRP_metabolites <- best_adj_hmdb
colnames(CRP_metabolites)[17] <-'Disease_Measure'
CRP_metabolites$Disease_Measure_Type <- 'CRP'

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
  theme_minimal()+
  xlim(-200,10)

# ESR----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

best_augmented <- top_50_peaks %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

ESR_metabolites <- best_adj_hmdb
colnames(ESR_metabolites)[17] <-'Disease_Measure'
ESR_metabolites$Disease_Measure_Type <- 'ESR'

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

# HAQ----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

best_augmented <- top_50_peaks %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

HAQ_metabolites <- best_adj_hmdb
colnames(HAQ_metabolites)[17] <-'Disease_Measure'
HAQ_metabolites$Disease_Measure_Type <- 'HAQ'

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
       y='ΔPeak Intensity')+
  theme_minimal()

# PVAS----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

top_20_peaks <- head(bestest_fit, 20)
top_50_peaks <- head(bestest_fit, 50)
top_100_peaks <- head(bestest_fit, 100)

best_augmented <- top_50_peaks %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented)
best_augmented_sig <- subset(best_augmented, best_augmented$p.value< 0.05)
best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)
adj_p <- p.adjust(best_augmented_sig$p.value, method='BH')
best_augmented_sig$adj_p <- adj_p

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Peak_ID %in% sig_peaks$Peak_ID)

PVAS_metabolites <- best_adj_hmdb
colnames(PVAS_metabolites)[17] <-'Disease_Measure'
PVAS_metabolites$Disease_Measure_Type <- 'PVAS'

ggplot(best_adj_hmdb,aes(x = PVAS, y=Peak_Intensity)) +
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
  labs(x='ΔPVAS',
       y='ΔPeak Intensity')+
  theme_minimal()
# GHVAS ----
ints_nested <- resp_melt_top %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_melt_top %>%
  unnest(cols=c())
identical(resp_diff_melt, ints_unnested)
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

best_adj <- subset(best_augmented_sig, best_augmented_sig$adj_p < 0.05) 
best_adj_hmdb <- inner_join(best_adj, peak_ID_HMDB, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
sig_peaks <- distinct(best_adj_hmdb, Peak_ID, .keep_all = TRUE)

ghvas_metabolites <-best_adj_hmdb 
colnames(ghvas_metabolites)[17] <-'Disease_Measure'
ghvas_metabolites$Disease_Measure_Type <- 'GHVAS'

ggplot(best_adj_hmdb,aes(x = GHVAS, y=Peak_Intensity)) +
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
  labs(x='ΔGHVAS',
       y='ΔPeak Intensity')+
  theme_minimal()

### Find shared metabolites across the above disease measure associations and the feature selection from RF -----
# first, which metabolites appear more than once across the disease measures?
mets_reduce <- function(disease_df){
  disease_met <- disease_df[c(28:29)]
  disease_met <- distinct(disease_met, Putative_Metabolite, .keep_all = TRUE)
}

DAS_met <- mets_reduce(DAS_metabolites)
CRP_met <- mets_reduce(CRP_metabolites)
ESR_met <- mets_reduce(ESR_metabolites)
HAQ_met <- mets_reduce(HAQ_metabolites)
PVAS_met <- mets_reduce(PVAS_metabolites)
GHVAS_met <- mets_reduce(ghvas_metabolites)

shared_mets <- rbind.data.frame(DAS_metabolites,
                                CRP_metabolites,
                                ESR_metabolites,
                                HAQ_metabolites,
                                PVAS_metabolites, 
                                ghvas_metabolites)

shared_mets_reduced <- rbind.data.frame(DAS_met,
                                        CRP_met,
                                        ESR_met,
                                        HAQ_met,
                                        PVAS_met,
                                        GHVAS_met)

metabolite_counts <- as.data.frame(as.matrix(table(shared_mets_reduced$Putative_Metabolite, shared_mets_reduced$Disease_Measure_Type)))
names(metabolite_counts)<- c('Putative_Metabolite', 'Disease_Measure','Frequency')

metabolite_counts%>%
  ggplot()+
  geom_col(aes(reorder(Putative_Metabolite, Frequency),
               Frequency,
               fill=Disease_Measure))+
  theme_pubclean()+
  coord_flip()+
  labs(x='Putative Metabolite',
       y='Frequency')


