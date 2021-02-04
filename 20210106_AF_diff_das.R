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
AF_limma$Sig <- ifelse(AF_limma$adj.P.Val <0.05, 'Significant', 'Not significant') 
AF_limma$Sig_Peaks <- ifelse(AF_limma$adj.P.Val<0.001 & AF_limma$identification != '', AF_limma$Peak_ID, '')

AF_limma_hmdb <- inner_join(AF_limma, peak_ID_HMDB, by='Peak_ID')
AF_limma_hmdb$Sig_Peaks <- ifelse(AF_limma_hmdb$adj.P.Val<0.05, AF_limma_hmdb$Putative_Metabolite, '')

AF_limma_hmdb_dist <- distinct(AF_limma_hmdb, Putative_Metabolite, .keep_all = TRUE)

AF_limma_hmdb_dist$Sig <- as.factor(AF_limma_hmdb_dist$Sig)

ggplot(data=AF_limma_hmdb_dist, aes(x=logFC, y=-log10(P.Value), 
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
  scale_color_brewer(palette = "Set1")+
  ylim(0,17)+
  xlim(-2,5)

# Histogram of p-values 
alpha <- binw <- 0.05 #where α = 0.05
pi0 <- 2*mean(AF_limma$P.Value > 0.05)
sample_size <- nrow(resp_differential_t)
AF_limma%>%
  mutate(Colour = ifelse(adj.P.Val < 0.05, "adj.P.Val < 0.05", "adj.P.Val > 0.05")) %>%
  ggplot(aes(P.Value))+
  geom_histogram(aes(fill=Colour),binwidth=binw, boundary=0, colour='black')+
  theme_pubclean()+
  theme(legend.title = element_blank())+
  geom_hline(yintercept = pi0 * binw * sample_size, col='blue')+
  geom_vline(xintercept=alpha, col='red')+
  labs(x='p-value',
       y='Frequency')

yintercept <- pi0 * binw * sample_size
count(AF_limma$P.Value < 0.025)
FDR_fraction <- pi0 * alpha/mean(AF_limma$P.Value <= alpha)

# Adding local FDR, qvalue and π0 values to p-value histogram
library(qvalue)
pi0 <- 2*mean(AF_limma$P.Value > 0.05)
lfdrvals <- lfdr(AF_limma$P.Value, pi0)
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

library(parallel)
library(doParallel)
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
# train model
control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=5,
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
imp_peaks <- varImp(custom, scale = TRUE)
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

write.csv(imps_hmdb_id, '20210111_AF_diff_das_FI_metabolites.csv')

# Plot the annotated peaks from feature importance
imps_hmdb_id <- read.csv('20210111_AF_diff_das_FI_metabolites.csv')
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

write.csv(shared_mets_reduced, '20210118_AF_diff_das_mets_count.csv')
shared_mets_reduced <- read.csv('20210118_AF_diff_das_mets_count.csv')
metabolite_counts%>%
  ggplot()+
  geom_col(aes(reorder(Putative_Metabolite, Frequency),
               Frequency,
               fill=Disease_Measure))+
  theme_pubclean()+
  coord_flip()+
  labs(x='Putative Metabolite',
       y='Frequency')

### Correlation analysis- which metabolites correlate with each other?
# Dendrogram
resp_diff_cor <- resp_diff[,c(14:1471)]
resp_diff_cor_t <- as.data.frame(t(resp_diff_cor))
resp_diff_cor_t$Peak_ID <- rownames(resp_diff_cor_t)
resp_diff_cor_t$Peak_ID <- gsub('X','', resp_diff_cor_t$Peak_ID)
resp_diff_cor_t$Peak_ID <- as.numeric(resp_diff_cor_t$Peak_ID)

resp_diff_cor_t_hmdb <- inner_join(resp_diff_cor_t,peak_ID_HMDB, by='Peak_ID')
resp_diff_cor_t_hmdb <- subset(resp_diff_cor_t_hmdb, resp_diff_cor_t_hmdb$Putative_Metabolite !='NA')
resp_diff_cor_t_hmdb <- distinct(resp_diff_cor_t_hmdb, Putative_Metabolite, .keep_all = TRUE)

metabolite_list <- rownames(resp_diff_cor_t_hmdb)
rownames(resp_diff_cor_t_hmdb) <- resp_diff_cor_t_hmdb$Putative_Metabolite
resp_diff_cor_t_hmdb <- resp_diff_cor_t_hmdb[1:63]

hc_complete <- hclust(dist(resp_diff_cor_t_hmdb), method='complete')
hc_average <- hclust(dist(resp_diff_cor_t_hmdb), method='average')
hc_single <- hclust(dist(resp_diff_cor_t_hmdb), method='single')

par(mfrow=c(1,3)) 
plot(hc_complete,main='Complete Linkage', xlab='', ylab='', sub='',cex=.9)
plot(hc_average,main='Average Linkage', xlab='', ylab='', sub='',cex=.9)
plot(hc_single,main='Single Linkage', xlab='', ylab='', sub='',cex=.9)
dev.off()
plot(hc_complete,main='Complete Linkage', xlab='', ylab='', sub='',cex=.9,
     )

library(ggdendro)
dendr <- dendro_data(hc_complete, type="rectangle") 
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + 
  scale_y_reverse(expand=c(0.2, 1)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

# Heatmap
resp_diff_cor_heat <- t(resp_diff_cor_t_hmdb)
cormat <- round(cor(resp_diff_cor_heat),2)
melted_cormat <- melt(cormat)

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab") +
theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size =8, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

dev.off()

### Kegg pathway analysis
library(FELLA)
met_list_kegg <- read.csv('metabolite.list.csv', header=TRUE)
met_list_kegg$kegg[met_list_kegg$kegg ==''] <- NA

kegg_list <- as.character(met_list_kegg$kegg)
set.seed(1)
# Filter the hsa01100 overview pathway
graph <- buildGraphFromKEGGREST(organism = "hsa", filter.path = c("01100"))
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(
  keggdata.graph = graph,
  databaseDir = tmpdir,
  internalDir = FALSE,
  matrices = "none",
  normality = "diffusion",
  niter = 100)

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "none"
)

kegg_list<- na.omit(kegg_list)
analysis_taser <- enrich(
  compounds = kegg_list,
  data = fella.data,
  method = "diffusion",
  approx = "normality")

analysis_taser %>%
  getInput %>%
  getName(data = fella.data)

g_taser <-  generateResultsGraph(
  object = analysis_taser,
  data = fella.data,
  method = "diffusion")

tab_taser <- generateResultsTable(
  object = analysis_taser,
  data = fella.data,
  method = "diffusion")

pathways_taser <- tab_taser[tab_taser$Entry.type == "pathway", ]
pathways_taser_flex <- flextable(pathways_taser)%>%
  theme_vanilla()%>%
  autofit(add_w = 0, add_h = 0)%>%
  align(align = 'center', part= 'all')

g_go <- addGOToGraph(graph=g_taser,
                     GOterm = "GO:0005739",
                     godata.options = list(
                       OrgDb = "org.Hs.eg.db", ont = "CC"),
                     mart.options = list(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))


set.seed(42)
plot(g_go,
     data = fella.data,
     vertex.size=2,
     plotLegend = TRUE,
     vertex.color=c("gold","red","green"),
     vertex.frame.color="gray", 
     vertex.label.color="black", 
     vertex.label.cex=.7, 
     vertex.label.dist=50,
     edge.curved=.5,
     edge.arrow.size=2,
     method = "diffusion",
     plimit=20,
     vertex_disjoint_paths(),
     margin=0) 

png("network_taser.png", units="px", width=1000, height=3000, res=300)

analysis.taser <- defineCompounds(compounds=kegg_list, data=fella.data)
analysis.taser<- runDiffusion(object=analysis.taser,
                              data=fella.data,
                              approx="normality")
analysis.taser
nlimit <- 150
vertex.label.cex <- .5

ab <- plot.igraph(
  analysis.taser,
  method='diffusion',
  data=fella.data,
  nlimit=250,
  plimit=25,
  vertex.label.cex= vertex.label.cex)

g <- generateResultsGraph(
object = analysis.taser,
method = "diffusion",
nlimit = nlimit,
data = fella.data)
g
g_go <- addGOToGraph(graph=g,
                     GOterm = "GO:0005739",
                     godata.options = list(
                       OrgDb = "org.Hs.eg.db", ont = "CC"),
                     mart.options = list(biomart = "ensembl", 
                                         dataset = "hsapiens_gene_ensembl"))

set.seed(42)
plot_graph(g_go)

tab.all <- generateResultsTable(
  method='diffusion',
  nlimit=100,
  object=analysis.taser,
  data=fella.data)

tab.enzyme <- generateEnzymesTable(
  method='diffusion',
  nlimit=100,
  object=analysis.taser,
  data=fella.data)

# which peaks from the feature importance are found in the shared metabolites list?
fi_cor_mets <- subset(imps_hmdb_id, imps_hmdb_id$Peak_ID %in% shared_mets$Peak_ID)

ggplot(fi_cor_mets)+
  geom_col(aes(reorder(Putative_Metabolite, Importance), 
               Importance),
           fill=0x3a5e84,
           colour='black')+
  coord_flip()+
  theme_minimal()+
  labs(y='Relative Importance',
       x='Putative Metabolite')+
  theme(axis.text.y = element_text(size = 8))

names(met_list_kegg)[2] <- 'Putative_Metabolite'
met_list_kegg$Putative_Metabolite <- as.factor(met_list_kegg$Putative_Metabolite)
fi_cor_kegg <- inner_join(fi_cor_mets,met_list_kegg, by='Putative_Metabolite')
update_kegg_shared <- fi_cor_kegg
update_kegg_shared$kegg[update_kegg_shared$kegg ==''] <- NA
update_kegg_list <- update_kegg_shared$kegg
#write.csv(update_kegg_shared,'20210108_AF_diff_das_kegg_list.csv')
ab <- read.csv('20210108_AF_diff_das_kegg_list.csv')
analysis_update <- enrich(
  compounds = update_kegg_list,
  data = fella.data,
  method = "diffusion",
  approx = "normality")

analysis_update %>%
  getInput %>%
  getName(data = fella.data)

g_update <-  generateResultsGraph(
  object = analysis_update,
  data = fella.data,
  plot.fun= tkplot,
  method = "diffusion")

tab_taser_2 <- generateResultsTable(
  object = analysis_update,
  data = fella.data,
  method = "diffusion")

update_pathways <- tab_taser_2[tab_taser_2$Entry.type == "pathway", ]
update_pathways_print <- update_pathways[,-2]

pathway_flex <- flextable(update_pathways_print)%>%
  theme_vanilla()%>%
  autofit(add_w = 0, add_h = 0)%>%
  align(align = 'center', part= 'all')


g_go <- addGOToGraph(graph=g_update,
                     GOterm = "GO:0005739",
                     godata.options = list(
                       OrgDb = "org.Hs.eg.db", ont = "CC"),
                     mart.options = list(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))

set.seed(42)
plot_graph(g_update)

tab.all <- generateResultsTable(
  method='diffusion',
  nlimit=100,
  object=analysis_update,
  data=fella.data)

tab.enzyme <- generateEnzymesTable(
  method='diffusion',
  nlimit=100,
  object=analysis_update,
  data=fella.data)

# attempt to split the pathway analysis, plotting only the pathways of interest
# Alanine, aspartate and glutamate metabolism
# Filter the hsa01100 overview pathway
graph <- buildGraphFromKEGGREST(organism = "hsa", filter.path = c("00220"))
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(
  keggdata.graph = graph,
  databaseDir = tmpdir,
  internalDir = FALSE,
  matrices = "none",
  normality = "diffusion",
  niter = 100)

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "none"
)

analysis_update <- enrich(
  compounds = update_kegg_list,
  data = fella.data,
  method = "diffusion",
  approx = "normality")

?enrich
analysis_update %>%
  getInput %>%
  getName(data = fella.data)

g_update <-  generateResultsGraph(
  object = analysis_update,
  data = fella.data,
  method = "diffusion")

tab_taser_2 <- generateResultsTable(
  object = analysis_update,
  data = fella.data,
  method = "diffusion")

update_pathways <- tab_taser_2[tab_taser_2$Entry.type == "pathway", ]
update_pathways_print <- update_pathways[,-2]

pathway_flex <- flextable(update_pathways_print)%>%
  theme_vanilla()%>%
  autofit(add_w = 0, add_h = 0)%>%
  align(align = 'center', part= 'all')


g_go <- addGOToGraph(graph=g_update,
                     GOterm = "GO:0005739",
                     godata.options = list(
                       OrgDb = "org.Hs.eg.db", ont = "CC"),
                     mart.options = list(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))

set.seed(42)
plot_graph(g)

tab.all <- generateResultsTable(
  method='diffusion',
  nlimit=100,
  object=analysis,
  data=fella.data)

tab.enzyme <- generateEnzymesTable(
  method='diffusion',
  nlimit=100,
  object=analysis,
  data=fella.data)

### Directly investigating differential abundance of metabolites of interest
resp_ints_melt <- resp_ints[c(3,9:1466)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Sample','Peak_ID','Peak_Intensity')

stat_test <- resp_ints_melt %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Sample) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


resp_ints_stats <- inner_join(resp_ints_melt, stat_test, by='Peak_ID')
resp_ints_stats$Peak_ID_2 <- resp_ints_stats$Peak_ID
resp_ints_stats$Peak_ID <- gsub('X', '', resp_ints_stats$Peak_ID)
resp_ints_stats$Peak_ID <- as.numeric(resp_ints_stats$Peak_ID)

resp_ints_stats_hmdb <- inner_join(resp_ints_stats, peak_ID_HMDB, by='Peak_ID')
resp_ints_stats_hmdb_id <- subset(resp_ints_stats_hmdb, resp_ints_stats_hmdb$Putative_Metabolite !='NA')  

resp_ints_stats_dist_A <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Sample == 'A')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_dist_F <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Sample == 'F')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_hmdb_comb <- rbind.data.frame(resp_ints_stats_dist_A ,resp_ints_stats_dist_F)

write.csv(resp_ints_stats_hmdb, '20210111_AF_diff_das.csv')
resp_ints_stats_hmdb <- read.csv('20210111_AF_diff_das.csv')
resp_ints_stats_hmdb_comb$Sample <- as.factor(resp_ints_stats_hmdb_comb$Sample)
resp_ints_stats_hmdb_comb$Peak_Intensity <- as.numeric(resp_ints_stats_hmdb_comb$Peak_Intensity)
resp_ints_stats_hmdb_comb$p.adj <- signif(resp_ints_stats_hmdb_comb$p.adj,3)
ab <-resp_ints_stats_hmdb_comb %>%
  subset(resp_ints_stats_hmdb_comb$p.adj < 0.05)
ab <- as.data.frame(ab[c(2,17)])
ab <- distinct(ab, Putative_Metabolite, .keep_all = TRUE)


diff_analysis <- function(peak_ID, position){
  resp_ints_stats_hmdb_id %>%
    subset(Peak_ID == `peak_ID`)%>%
    ggplot(aes(Sample, Peak_Intensity,
               fill=Sample))+
    theme_light()+
    geom_violin()+
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    
    stat_pvalue_manual(
      resp_ints_stats_hmdb_comb[resp_ints_stats_hmdb_comb$Peak_ID == `peak_ID`,], 
      y.position = position,
      label = "p.adj")+
    theme(legend.position = 'none')+
    labs(y='Peak Intensity')
}

  diff_analysis(169, 25.2)

  

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

### Plotting the patient groups into response groups based on potential biomarkers------
# Looking at the changing metabolites and their involvement in responses??
# First plot boxplots side by side to show the change in the mean peak intensity across response groups
kegg_list <- read.csv('metabolite.list.csv', header=TRUE)
peak_ID_HMDB_select <- subset(peak_ID_HMDB, peak_ID_HMDB$Putative_Metabolite %in% kegg_list$x)

resp_diff_box <- resp_diff_2[,-1]
resp_diff_box_t <- as.data.frame(t(resp_diff_box))
resp_diff_box_t$Peak_ID <- rownames(resp_diff_box_t)
resp_diff_box_t <- resp_diff_box_t[,c(ncol(resp_diff_box_t),1:(ncol(resp_diff_box_t)-1))]
resp_diff_box_t$Peak_ID <- gsub('X','', resp_diff_box_t$Peak_ID)
resp_diff_box_t$Peak_ID <- as.numeric(resp_diff_box_t$Peak_ID)

resp_diff_selects <- subset(resp_diff_box_t, resp_diff_box_t$Peak_ID %in% peak_ID_HMDB_select$Peak_ID)
resp_diff_selects_t <- resp_diff_selects[,-1]
resp_diff_selects_t <- as.data.frame(t(resp_diff_selects_t))
resp_diff_selects_t$DAS44 <- resp_diff$DAS44_Response
resp_diff_selects_t$DAS44[resp_diff_selects_t$DAS44 == 'Good' | resp_diff_selects_t$DAS44 == 'Remission']<- 'Positive'
resp_diff_selects_t$DAS44[resp_diff_selects_t$DAS44 == 'Poor' | resp_diff_selects_t$DAS44 == 'Mild']<- 'Negative'

resp_diff_selects_t$Sample <- rownames(resp_diff_selects_t)
resp_melt <- melt(resp_diff_selects_t)
names(resp_melt)[3:4] <- c('Peak_ID', 'Peak_Intensity')

stat_test <- resp_melt %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ DAS44) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Plot showing the mean Δpeak intensities across the selected peaks in positive and negative responders
ggplot(resp_melt,aes(x=DAS44, y=Peak_Intensity))+
  geom_violin(aes(fill=DAS44))+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  theme_minimal()+
  stat_compare_means(method= 'wilcox.test',
                     label = "p.format",
                     vjust=1, 
                     hjust=-1)+
  theme(legend.position= 'none')+
  labs(x='DAS44 Response',
       y='ΔPeak Intensity')

### Attempt above with the baseline levels against the differential response-----
### Directly investigating differential abundance of metabolites of interest
resp_ints_melt <- resp_ints[c(8:1466)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Response','Peak_ID','Peak_Intensity')
resp_ints_melt$Peak_ID <- gsub('X', '', resp_ints_melt$Peak_ID)

resp_ints_melt_select <- subset(resp_ints_melt, resp_ints_melt$Peak_ID %in% peak_ID_HMDB_select$Peak_ID)

stat_test <- resp_ints_melt_select %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Response) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

resp_ints_stats <- inner_join(resp_ints_melt_select, stat_test, by='Peak_ID')
resp_ints_stats$Peak_ID_2 <- resp_ints_stats$Peak_ID
resp_ints_stats$Peak_ID <- gsub('X', '', resp_ints_stats$Peak_ID)
resp_ints_stats$Peak_ID <- as.numeric(resp_ints_stats$Peak_ID)

resp_ints_stats_hmdb <- inner_join(resp_ints_stats, peak_ID_HMDB, by='Peak_ID')
resp_ints_stats_hmdb_id <- subset(resp_ints_stats_hmdb, resp_ints_stats_hmdb$Putative_Metabolite !='NA')  

resp_ints_stats_dist_pos <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Response == 'Positive')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_dist_neg <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Response == 'Negative')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_hmdb_comb <- rbind.data.frame(resp_ints_stats_dist_pos ,resp_ints_stats_dist_neg)
resp_ints_stats_hmdb_comb$Response <- as.factor(resp_ints_stats_hmdb_comb$Response)
resp_ints_stats_hmdb_comb$Peak_Intensity <- as.numeric(resp_ints_stats_hmdb_comb$Peak_Intensity)
resp_ints_stats_hmdb_comb$p.adj <- signif(resp_ints_stats_hmdb_comb$p.adj,3)

resp_ints_stats_hmdb_comb %>%
  ggplot(aes(Response, Peak_Intensity,
             fill=Response))+
  theme_light()+
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  stat_compare_means(method= 't.test',
                     label = "p.format",
                     vjust=1, 
                     hjust=-1)+
  theme(legend.position = 'none')+
  labs(y='Peak Intensity')

# Keep directional metabolites, based on fold change from the AF_limma -----
AF_limma_pos <- subset(AF_limma, AF_limma$logFC >0)
AF_limma_neg <- subset(AF_limma, AF_limma$logFC <0)
peak_ID_HMDB_select_pos <-subset(peak_ID_HMDB_select, peak_ID_HMDB_select$Peak_ID %in% AF_limma_pos$Peak_ID)
peak_ID_HMDB_select_neg <-subset(peak_ID_HMDB_select, peak_ID_HMDB_select$Peak_ID %in% AF_limma_neg$Peak_ID)

resp_ints_melt_select <- subset(resp_ints_melt, resp_ints_melt$Peak_ID %in% peak_ID_HMDB_select_pos$Peak_ID)

stat_test <- resp_ints_melt_select %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Response) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

resp_ints_stats <- inner_join(resp_ints_melt_select, stat_test, by='Peak_ID')
resp_ints_stats$Peak_ID_2 <- resp_ints_stats$Peak_ID
resp_ints_stats$Peak_ID <- gsub('X', '', resp_ints_stats$Peak_ID)
resp_ints_stats$Peak_ID <- as.numeric(resp_ints_stats$Peak_ID)

resp_ints_stats_hmdb <- inner_join(resp_ints_stats, peak_ID_HMDB, by='Peak_ID')
resp_ints_stats_hmdb_id <- subset(resp_ints_stats_hmdb, resp_ints_stats_hmdb$Putative_Metabolite !='NA')  

resp_ints_stats_dist_pos <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Response == 'Positive')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_dist_neg <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Response == 'Negative')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_hmdb_comb <- rbind.data.frame(resp_ints_stats_dist_pos ,resp_ints_stats_dist_neg)
resp_ints_stats_hmdb_comb$Response <- as.factor(resp_ints_stats_hmdb_comb$Response)
resp_ints_stats_hmdb_comb$Peak_Intensity <- as.numeric(resp_ints_stats_hmdb_comb$Peak_Intensity)
resp_ints_stats_hmdb_comb$p.adj <- signif(resp_ints_stats_hmdb_comb$p.adj,3)

resp_ints_stats_hmdb_comb %>%
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
  
## Looking at the change in peak intensity of selected peaks for their change from A to F ------
### Attempt above with the baseline levels against the differential response
### Directly investigating differential abundance of metabolites of interest
resp_ints_melt <- resp_ints[c(3,9:1466)]
resp_ints_melt <- melt(resp_ints_melt)
names(resp_ints_melt) <- c('Sample','Peak_ID','Peak_Intensity')
resp_ints_melt$Peak_ID <- gsub('X', '', resp_ints_melt$Peak_ID)

resp_ints_melt_select <- subset(resp_ints_melt, resp_ints_melt$Peak_ID %in% peak_ID_HMDB_select$Peak_ID)

stat_test <- resp_ints_melt_select %>%
  group_by(Peak_ID) %>%
  rstatix::wilcox_test(Peak_Intensity ~ Sample) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

resp_ints_stats <- inner_join(resp_ints_melt_select, stat_test, by='Peak_ID')
resp_ints_stats$Peak_ID_2 <- resp_ints_stats$Peak_ID
resp_ints_stats$Peak_ID <- gsub('X', '', resp_ints_stats$Peak_ID)
resp_ints_stats$Peak_ID <- as.numeric(resp_ints_stats$Peak_ID)

resp_ints_stats_hmdb <- inner_join(resp_ints_stats, peak_ID_HMDB, by='Peak_ID')
resp_ints_stats_hmdb_id <- subset(resp_ints_stats_hmdb, resp_ints_stats_hmdb$Putative_Metabolite !='NA')  

resp_ints_stats_dist_A <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Sample == 'A')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_dist_ <- subset(resp_ints_stats_hmdb_id, resp_ints_stats_hmdb_id$Sample == 'F')  %>%
  distinct(Peak_ID, .keep_all = TRUE)

resp_ints_stats_hmdb_comb <- rbind.data.frame(resp_ints_stats_dist_A ,resp_ints_stats_dist_F)
resp_ints_stats_hmdb_comb$Sample <- as.factor(resp_ints_stats_hmdb_comb$Sample)
resp_ints_stats_hmdb_comb$Peak_Intensity <- as.numeric(resp_ints_stats_hmdb_comb$Peak_Intensity)
resp_ints_stats_hmdb_comb$p.adj <- signif(resp_ints_stats_hmdb_comb$p.adj,3)
resp_ints_stats_hmdb_padj <- subset(resp_ints_stats_hmdb_comb, resp_ints_stats_hmdb_comb$p.adj < 1e-3)

resp_ints_stats_hmdb_padj %>%
  ggplot(aes(Sample, Peak_Intensity,
             fill=Sample))+
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
resp_int_lr <- resp_diff_ML_2

das_peaks <- c("L-Glutamate", "Hypoxanthine", "L-Valine", "Uridine", "L-Ornithine",
               "Creatinine", "Orotate" ,  "N-Acetylornithine")
das_strict <- c("L-Glutamate", "Hypoxanthine", "L-Histidine", "L-Valine", "Traumatic acid", 'Orotate', "N-Acetylornithine")

peaks_das <- subset(peak_ID_HMDB, peak_ID_HMDB$Putative_Metabolite %in% das_strict)

AF_refined <- subset(AF_limma_hmdb, AF_limma_hmdb$adj.P.Val < 0.01)
AF_refined <- subset(AF_refined, AF_refined$Putative_Metabolite != 'NA')
AF_refined <- distinct(AF_refined, Putative_Metabolite, .keep_all = TRUE)

resp_lr <- resp_int_lr[,-1]
resp_lr_t <- as.data.frame(t(resp_lr))
resp_lr_t$Peak_ID <- rownames(resp_lr_t)
resp_lr_sel <- subset(resp_lr_t, resp_lr_t$Peak_ID %in% peaks_das$Peak_ID)
resp_lr_sel_t <- resp_lr_sel[,-64]
resp_lr_sel_t <- as.data.frame(t(resp_lr_sel_t))
resp_lr_mets<- resp_lr_sel_t
resp_lr_mets$Response <- resp_int_lr$DAS44_Response
resp_lr_mets <- resp_lr_mets[,c(ncol(resp_lr_mets),1:(ncol(resp_lr_mets)-1))]

resp_lr_mets$Response <- as.character(resp_lr_mets$Response)
resp_lr_mets$Response[resp_lr_mets$Response=='Positive'] <- 1
resp_lr_mets$Response[resp_lr_mets$Response=='Negative'] <- 0
resp_lr_mets$Response <- as.numeric(resp_lr_mets$Response)
names(resp_lr_mets)[2:13] <- paste0('X', names(resp_lr_mets)[2:13])

set.seed(42)
index <- createDataPartition(resp_lr_mets$Response, p = 0.8, list = FALSE) # 0.8
train_data <- resp_lr_mets[index, ]
test_data  <- resp_lr_mets[-index, ]

model <- glm(Response ~.,family=binomial(link='logit'),data=train_data)
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
model <- glm(Response ~ X270 + X711 +X738 ,family=binomial(link='logit'),data=train_data)
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


