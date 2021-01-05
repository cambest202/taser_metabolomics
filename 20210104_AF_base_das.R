## AF baseline analysis and association of metabolites with DAS44. Standardising the approach

### ggplot theme------
theme_CHIPS <- function () { 
  theme_bw() %+replace% 
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.background = element_blank(), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      plot.title = element_text(size=16, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
      title = element_text(size = 16, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
      axis.text.y = element_text(size = 14, margin = margin(r = 5),hjust=1,vjust=0.5, family="Arial", face="bold",colour="black"),
      axis.text.x = element_text(size = 14, margin = margin(t = 5),hjust=0.5,vjust=1, family="Arial", face="bold",colour="black"), 
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
peak_metadata <- read.csv(file='20200427_Taser_PeakMetadata.csv', header=TRUE, row.names=1)

sample_sheet = read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)

names(sample_sheet)[2] <- 'Sample_Name'
peak_IDs$Peak_ID <- rownames(peak_IDs)
peak_IDs$Peak_ID  <- as.numeric(peak_IDs$Peak_ID )
peak_IDs$moleculeName  <- as.character(peak_IDs$moleculeName )
colnames(peak_metadata)[6] <- 'Peak_ID'
peak_metadata <- left_join(peak_metadata, peak_IDs)
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
colnames(resp_subset)[2:43] <- paste0('X', colnames(resp_subset)[2:43])

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
AF_limma$Sig <- ifelse(AF_limma$adj.P.Val <0.05, TRUE, FALSE) 
AF_limma$Sig_Peaks <- ifelse(AF_limma$adj.P.Val<0.05 & AF_limma$identification != '', AF_limma$Peak_ID, '')

# Histogram of p-values 
histogram(AF_limma$P.Value) #histogram shows a uniform distribution of p-values for the AF_base_das analysis.
# H0 is that there is no difference in the metabolite levels for the negative and positive response groups in the population
# after multiple testing correction, only a small number would likely show up as being significant.

# Generate the model
# Set training and testing data
set.seed(42)
index <- createDataPartition(resp_subset$Response, p = 0.85, list = FALSE)
train_data <- resp_subset[index, ]
test_data  <- resp_subset[-index, ]

tuneGrid <- expand.grid(.mtry = c(1:sqrt(1458)))
set.seed(42)
model_rf <- caret::train(Response ~ .,
                         data = train_data,
                         method = "rf",
                         preProcess = c("scale", "center"),
                         ntree = 1000,
                         metric='Accuracy',
                         tuneGrid = tuneGrid,
                         trControl =  trainControl(method = "repeatedcv",
                                                   number = 10,
                                                   repeats = 10, 
                                                   savePredictions = TRUE, 
                                                   verboseIter = FALSE, 
                                                   allowParallel = TRUE))
 
plot(model_rf)
print(model_rf)

test_results <- predict(model_rf, newdata = test_data)
summary(test_results)
summary(test_data$Response)
confusionMatrix(test_results, test_data$Response)
 
 
