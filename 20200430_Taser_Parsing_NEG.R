### Parsing TaSER data: Positive and Negative Data

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
## Libraries
library(ggplot2)
library(tidyverse)
library (readr)
library(ggpubr)
library (mosaic)
library (dplyr)
library (data.table)
library(reshape2)
library (gtools)
library (ggstance)
library(plyr)
library(limma)
library(ggrepel)
library(flextable)
library('vsn')


## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
## Set working directory
#setwd ("/Users/cameronbest/University/MRC_DTP_in_Precision_Medicine/Project/Taser_data/Analysis/EDA_May2020/Parsing")

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
## Loading data
negative_peaks <- read.table(file="20160324_TaserAnalysisNeg.txt", sep="\t", quote="", comment.char="", header=TRUE)
sample_sheet <- read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
peak_IDs <- read.table (file="20200117_Taser_PeakIDs.csv", header=TRUE, row.names=1, sep= "\t")

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
## ggplot theme
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

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
### Peak metadata
# negative
neg_peak_data <- negative_peaks[,c(1:3, 384:388)]
colnames(neg_peak_data)[1] <- 'Peak_ID'
neg_peak_data <- with(neg_peak_data, neg_peak_data[order(Peak_ID),])
neg_peak_data$RT <- neg_peak_data$RT/60
write.csv(neg_peak_data, "Peak_data.csv", row.names = FALSE)

## ——————————————————————————————————————————————————
### Peak ID
peak_IDs$Peak_ID <- rownames(peak_IDs)
peak_IDs$Peak_ID <- as.numeric(as.character(peak_IDs$Peak_ID ))

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
### Peak Data- take Log base 2 of the data 
total_peaks_neg <- (negative_peaks[, 4:382])
rownames(total_peaks_neg) <- negative_peaks$id

peaks_neg <- t(total_peaks_neg)
peaks_neg <- as.data.frame(subset(peaks_neg, rownames(peaks_neg) %like% 'A'|rownames(peaks_neg) %like% 'B'))
peaks_neg$Label <- 0
peaks_neg <- peaks_neg[,c(ncol(peaks_neg),1:(ncol(peaks_neg)-1))]
peaks_neg$No.ofSample <- rownames(peaks_neg)
peaks_neg <- peaks_neg[,c(ncol(peaks_neg),1:(ncol(peaks_neg)-1))]
peaks_neg_two <- as.data.frame(do.call('rbind', strsplit(as.character(peaks_neg$No.ofSample),'_',fixed=TRUE)))
peaks_neg <- cbind.data.frame(peaks_neg_two, peaks_neg)
peaks_neg_2 <- peaks_neg[,-1]
peaks_neg_2 <- peaks_neg_2[,-2]
colnames(peaks_neg_2)[1] <- 'No. of Sample'
peaks_neg_2$Label [rownames(peaks_neg_2) %like% 'A'] <- 'A'
peaks_neg_2$Label [rownames(peaks_neg_2) %like% 'B'] <- 'B'

peaks_neg_2$V2 <- 1:153

## Pre-processing
# Remove missing values, non-sample peaks 
neg_peaks <- total_peaks_neg

neg_peaks$Peak_ID <- rownames(neg_peaks)
neg_peaks$Peak_ID <- as.numeric(as.character(neg_peaks$Peak_ID ))
neg_peaks_2 <- with(neg_peaks, neg_peaks[order(Peak_ID),])
neg_peaks_2 <- neg_peaks[, 1:225] # only samples remain

replacezero <- function(x) "[<-"(x, !x|is.na(x), min(x[x>0], na.rm=TRUE)/2)
neg_peaks_2 <- t(neg_peaks_2) #transpose to allow above function
neg_peaks_2 <- apply(neg_peaks_2, 1, replacezero)

## Normalise
fit <- vsn2(neg_peaks_2)                   ## fit
meanSdPlot(fit)       # check fit
negpeaks = predict(fit, newdata=neg_peaks_2)  ## apply fit

# Split into AB, AF and ABF dataframes
total_samples_resp <- negpeaks
colnames(total_samples_resp) <- substring(colnames(total_samples_resp), 8)
total_samples_resp_t <- as.data.frame(t(total_samples_resp))

samples_AB_df <- subset(total_samples_resp_t, rownames(total_samples_resp_t) %like% 'A' | rownames(total_samples_resp_t) %like% 'B' )
write.csv(samples_AB_df, file='20210304_Taser_AB_NEG_PeakIntensities.csv')

samples_AF_df <- subset(total_samples_resp_t, rownames(total_samples_resp_t) %like% 'A' | rownames(total_samples_resp_t) %like% 'F' )
write.csv(samples_AF_df, file='20210304_Taser_AF_NEG_PeakIntensities.csv')

samples_ABF_df <- subset(total_samples_resp_t, rownames(total_samples_resp_t) %like% 'A' | rownames(total_samples_resp_t) %like% 'B' | rownames(total_samples_resp_t) %like% 'F' )
rownames(samples_ABF_df) <- substr(rownames(samples_ABF_df), 1,7)
write.csv(samples_ABF_df, file='20210304_Taser_ABF_NEG_PeakIntensities.csv')



#### Below is prior script, use the above for the recent versions of the dfs
neg_peaks_2 <- t(neg_peaks_2) # transpose again to tidy rownames
neg_peaks_AF <- neg_peaks_2
neg_peaks_ABF <- neg_peaks_2

neg_peaks_2 <- as.data.frame(subset(neg_peaks_2, rownames(neg_peaks_2) %like% 'A'|rownames(neg_peaks_2) %like% 'B'))
neg_peaks_AF <- as.data.frame(subset(neg_peaks_AF, rownames(neg_peaks_AF) %like% 'A'|rownames(neg_peaks_AF) %like% 'F'))
neg_peaks_ABF <- as.data.frame(subset(neg_peaks_ABF, rownames(neg_peaks_ABF) %like% 'A'| rownames(neg_peaks_ABF) %like% 'B'|rownames(neg_peaks_ABF) %like% 'F'))

neg_peaks_sample_list <- list(rownames(neg_peaks_2))
neg_peaks_AF_list <- list(rownames(neg_peaks_AF))
neg_peaks_ABF_list <- list(rownames(neg_peaks_ABF))

neg_peaks_3 <- t(neg_peaks_2)
neg_peaks_AF_t <- t(neg_peaks_AF)
neg_peaks_ABF_t <- t(neg_peaks_ABF)

colnames(neg_peaks_3) <- substring(colnames(neg_peaks_3), 8)
colnames(neg_peaks_AF_t) <- substring(colnames(neg_peaks_AF_t), 8)
colnames(neg_peaks_ABF_t) <- substring(colnames(neg_peaks_ABF_t), 8)

neg_peaks_3 <- neg_peaks_3
neg_peaks_AF <- neg_peaks_AF_t
neg_peaks_ABF <- neg_peaks_ABF_t

#Log transformation
log_neg_peaks <- log(neg_peaks_3)

# Normalisation
fit_AB <- vsn2(neg_peaks_3)                   ## fit
fit_AF <- vsn2(neg_peaks_AF)                   ## fit
fit_ABF <- vsn2(neg_peaks_ABF)                   ## fit

meanSdPlot(fit_ABF)

nnegpeaks_AB = predict(fit_AB, newdata=neg_peaks_3)  ## apply fit
nnegpeaks_AF = predict(fit_AF, newdata=neg_peaks_AF)  ## apply fit
nnegpeaks_ABF = predict(fit_ABF, newdata=neg_peaks_ABF)  ## apply fit

# Reassign rownames as the metabolites
neg_peaks_samples_AB <- as.data.frame(t(nnegpeaks_AB))
neg_peaks_samples_AF <- as.data.frame(t(nnegpeaks_AF))
neg_peaks_samples_ABF <- as.data.frame(t(nnegpeaks_ABF))

# Write csv files for the processed matrix
write.csv(neg_peaks_samples_AB, file='AB_test.csv')
write.csv(neg_peaks_samples_AF, file='AF_test.csv')
write.csv(neg_peaks_samples_ABF, file='BF_test.csv')


# log transformation
log_peaks <- log(nnegpeaks)

### Limma
# Negative
Group_neg <- factor(colnames(log_peaks), levels = c('A', 'B'))
design_neg <- model.matrix (~Group_neg)
colnames(design_neg) <- c('A', 'AvsB')
eset_neg <- log_peaks
fit_neg <- lmFit(eset_neg, design_neg)
fit_neg <- eBayes(fit_neg)
toptable_neg <- topTable(fit_neg, coef = 'AvsB', adjust = 'BH', number=385)
toptable_neg <- as.data.frame(toptable_neg)

# Identify the differentially produced metabolites 
toptable_neg$Peak_ID <- rownames(toptable_neg)
toptable_neg$Peak_ID <- as.numeric(as.character(toptable_neg$Peak_ID))

toptable_neg <- toptable_neg[,c(ncol(toptable_neg),1:(ncol(toptable_neg)-1))]

toptable_neg <- subset(toptable_neg, toptable_neg$Peak_ID %in% peak_IDs$Peak_ID) # remove peaks without standards
toptable_neg <- merge(toptable_neg, peak_IDs)

# Add the m/z and RT
toptable_neg2 <- inner_join(neg_peak_data, toptable_neg, by= 'Peak_ID')
toptable_neg2<- toptable_neg2[,-c(5:8,10,11,14)]
toptable_neg2 <- with(toptable_neg2, toptable_neg2[order(toptable_neg2$adj.P.Val),])

toptops <- head(toptable_neg2, n=10)

#Testing correct peak data
neg_peak_data_test <- subset(neg_peak_data, neg_peak_data$Peak_ID %in% toptable_neg$Peak_ID)

### Results
# Can use the neg_peak_samples in future experiments. Imputed missing values and normalised via VSN

