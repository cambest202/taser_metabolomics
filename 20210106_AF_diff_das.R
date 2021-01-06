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

# Generate the model
# Set training and testing data
set.seed(42)
index <- createDataPartition(resp_diff_ML_2$DAS44_Response, p = 0.85, list = FALSE)
train_data <- resp_diff_ML_2[index, ]
test_data  <- resp_diff_ML_2[-index, ]

# optimal ntrees? ----
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800 )) {
  set.seed(42)
  rf_maxtrees <- train(Response~.,
                       data = train_data,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trainControl(method = "repeatedcv",
                                                number = 10,
                                                repeats = 5, 
                                                savePredictions = TRUE, 
                                                verboseIter = FALSE, 
                                                allowParallel = TRUE),
                       importance = TRUE,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)

### Run the model ----
tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
set.seed(42)
model_rff <- caret::train(DAS44_Response~.,
                          data = train_data,
                          method = "rf",
                          metric = "Accuracy",
                          tuneGrid = tuneGrid,
                          trControl = trainControl(method = "repeatedcv",
                                                   number =5,
                                                   repeats = 3, 
                                                   savePredictions = TRUE, 
                                                   verboseIter = FALSE, 
                                                   allowParallel = TRUE),
                          importance = TRUE,
                          ntree = 500)

plot(model_rff)
print(model_rff)

test_results <- predict(model_rff, newdata = test_data)
summary(test_results)
summary(test_data$Response)
confusionMatrix(test_results, test_data$DAS44_Response)

final <- data.frame(actual = test_data$DAS44_Response,
                    predict(model_rff, newdata = test_data, type = "prob"))
final$predicted <- 0
final$predicted[final$Positive > final$Negative] <- 'Positive'
final$predicted[final$Positive < final$Negative] <- 'Negative'
final$correct <- 0
final$correct <- ifelse(final$predict == final$actual, 'Correct', 'Incorrect')
summary(final$predicted == final$actual)

imp <- model_rff$finalModel$importance
#imp <- imp[order(imp, decreasing = TRUE), ]
imp_peaks <- varImp(model_rff, scale = TRUE)
plot(imp_peaks, top=10)
imps <- as.data.frame(imp_peaks$importance)
imps$Peak_ID <- rownames(imps)
imps <- imps[,-1]
colnames(imps)[1] <- 'Importance'
imps <- subset(imps,imps$Importance >20)
imps$Peak_ID <- gsub('`','', imps$Peak_ID)
imps$Peak_ID <- as.numeric(imps$Peak_ID)

imps_hmdb <- inner_join(imps, peak_ID_HMDB, by='Peak_ID')
imps_hmdb <- with(imps_hmdb,imps_hmdb[order(Importance),])

# Using plotly for feature importance visualisation----
plotly_para <- function(){
  title_feat <- list(
    family = "Arial, sans-serif",
    size = 18,
    color = "black")
  text_feat <- list(
    family = "Old Standard TT, serif",
    size = 14,
    color = "black")
  y_axis <- list(
    title='Peak ID',
    titlefont = title_feat,
    showticklabels = TRUE,
    tickfont = text_feat)
  x_axis <- list(
    title='Importance',
    titlefont = title_feat,
    showticklabels = TRUE,
    tickfont = text_feat)
}

fig <-plot_ly(y=reorder(imps$Peak_ID, imps$Overall),
              x=imps$Overall, 
              data=imps,
              type='bar',
              orientation='h',
              marker = list(color = 'rgba(0, 70, 300, 0.6)',
                            line = list(color = 'rgba(50, 70, 96, 1.0)', width = 1)))
fig <- fig %>% 
  layout(yaxis = y_axis, xaxis = x_axis)
fig

### Build linear models for the top peaks from feature selection for each of the disease measurements
# DAS44
resp_melt_top <- subset(resp_A_melt, resp_A_melt$Peak_ID %in% imps_hmdb$Peak_ID)
# make sure all disease measures are included in the melted df

disease_lm <- function(melted_df, disease_measure){
  ints_nested <- melted_df %>%
    group_by (Peak_ID) %>%
    nest()
  ints_unnested <- ints_nested %>%
    unnest(cols=c())
  identical(melted_df, ints_unnested)
  ints_lm <- ints_nested %>%
    mutate(model = map(data, ~lm(formula = Peak_Intensity~disease_measure, data = .x)))
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
  best_augmented_sig$Peak_ID <- as.numeric(best_augmented_sig$Peak_ID)best_aug_sing <- distinct(best_augmented, Peak_ID, .keep_all = TRUE)
  adj_p <- p.adjust(best_aug_sing$p.value, method='BH')
  best_aug_sing$adj_p <- adj_p
  print(head(best_aug_sing))
}

ggplot(best_augmented_sig, aes(x = DAS44, y=Peak_Intensity)) +
  geom_point() + 
  theme_CHIPS()+
  stat_cor(method = "spearman", 
           vjust=1, hjust=0.1,
           size=4)+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Peak_ID, scales = "free_y")+
  theme(strip.background = element_rect(fill='white', 
                                        size=1.5),
        strip.text.x= element_text(face = "bold.italic",
                                   size=12))+
  labs(x='ΔDAS44',
       y='ΔPeak Intensity')


# CRP


