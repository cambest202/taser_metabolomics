### Parsing TaSER data: Positive and positive Data

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

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
## Loading data
positive_peaks <- read.table(file="20160324_TaserAnalysisPOS.txt", sep="\t", quote="", comment.char="", header=TRUE)
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
# positive
pos_peak_data <- positive_peaks[,c(1:3, 383:388)]
colnames(pos_peak_data)[1] <- 'Peak_ID'
pos_peak_data <- with(pos_peak_data, pos_peak_data[order(Peak_ID),])
pos_peak_data$RT <- pos_peak_data$RT/60
write.csv(pos_peak_data, "Pos_Peak_data.csv", row.names = FALSE)

## ——————————————————————————————————————————————————
### Peak ID
peak_IDs$Peak_ID <- rownames(peak_IDs)
peak_IDs$Peak_ID <- as.numeric(as.character(peak_IDs$Peak_ID ))

## ——————————————————————————————————————————————————
## ——————————————————————————————————————————————————
### Peak Data- take Log base 2 of the data 
total_peaks_pos <- (positive_peaks[, 4:382])
rownames(total_peaks_pos) <- positive_peaks$id

peaks_pos <- t(total_peaks_pos)
peaks_pos <- as.data.frame(subset(peaks_pos, rownames(peaks_pos) %like% 'A'|rownames(peaks_pos) %like% 'B'))
peaks_pos$Label <- 0
peaks_pos <- peaks_pos[,c(ncol(peaks_pos),1:(ncol(peaks_pos)-1))]
peaks_pos$No.ofSample <- rownames(peaks_pos)
peaks_pos <- peaks_pos[,c(ncol(peaks_pos),1:(ncol(peaks_pos)-1))]
peaks_pos_two <- as.data.frame(do.call('rbind', strsplit(as.character(peaks_pos$No.ofSample),'_',fixed=TRUE)))
peaks_pos <- cbind.data.frame(peaks_pos_two, peaks_pos)
peaks_pos_2 <- peaks_pos[,-1]
peaks_pos_2 <- peaks_pos_2[,-2]
colnames(peaks_pos_2)[1] <- 'No. of Sample'
peaks_pos_2$Label [rownames(peaks_pos_2) %like% 'A'] <- 'A'
peaks_pos_2$Label [rownames(peaks_pos_2) %like% 'B'] <- 'B'

peaks_pos_2$V2 <- 1:153

## Pre-processing
# Remove missing values, non-sample peaks 
pos_peaks <- total_peaks_pos

pos_peaks$Peak_ID <- rownames(pos_peaks)
pos_peaks$Peak_ID <- as.numeric(as.character(pos_peaks$Peak_ID ))
pos_peaks_2 <- with(pos_peaks, pos_peaks[order(Peak_ID),])
pos_peaks_2 <- pos_peaks[, 1:225] # only samples remain

replacezero <- function(x) "[<-"(x, !x|is.na(x), min(x[x>0], na.rm=TRUE)/2)
pos_peaks_2 <- t(pos_peaks_2) #transpose to allow above function
pos_peaks_2 <- apply(pos_peaks_2, 1, replacezero)
pos_peaks_2 <- pos_peaks_2[,c(1,2,1587)]

## Normalise
fit <- vsn2(pos_peaks_2)                   ## fit
meanSdPlot(fit)       # check fit
pospeaks = predict(fit, newdata=pos_peaks_2)  ## apply fit

# Split into AB, AF and ABF dataframes
total_samples_resp <- pospeaks
colnames(total_samples_resp) <- substring(colnames(total_samples_resp), 8)
total_samples_resp_t <- as.data.frame(t(total_samples_resp))

samples_AB_df <- subset(total_samples_resp_t, rownames(total_samples_resp_t) %like% 'A' | rownames(total_samples_resp_t) %like% 'B' )
write.csv(samples_AB_df, file='20210308_Taser_AB_POS_PeakIntensities.csv')

samples_AF_df <- subset(total_samples_resp_t, rownames(total_samples_resp_t) %like% 'A' | rownames(total_samples_resp_t) %like% 'F' )
write.csv(samples_AF_df, file='20210308_Taser_AF_POS_PeakIntensities.csv')

samples_ABF_df <- subset(total_samples_resp_t, rownames(total_samples_resp_t) %like% 'A' | rownames(total_samples_resp_t) %like% 'B' | rownames(total_samples_resp_t) %like% 'F' )
rownames(samples_ABF_df) <- substr(rownames(samples_ABF_df), 1,7)
write.csv(samples_ABF_df, file='20210308_Taser_ABF_POS_PeakIntensities.csv')



# Normalisation
fit = vsn2(pos_peaks_3)                   ## fit
pospeaks = predict(fit, newdata=pos_peaks_3)  ## apply fit

# Reassign rownames as the metabolites
pos_peaks_samples <- pospeaks
colnames(pos_peaks_samples) <- colnames(pos_peaks_2)
colnames(pos_peaks_samples) <- substring(colnames(pos_peaks_samples),8,12)

# log transformation
log_peaks <- log(npospeaks)

### Limma
# positive
Group_pos <- factor(colnames(log_peaks), levels = c('A', 'B'))
design_pos <- model.matrix (~Group_pos)
colnames(design_pos) <- c('A', 'AvsB')
eset_pos <- log_peaks
fit_pos <- lmFit(eset_pos, design_pos)
fit_pos <- eBayes(fit_pos)
toptable_pos <- topTable(fit_pos, coef = 'AvsB', adjust = 'BH', number=385)
toptable_pos <- as.data.frame(toptable_pos)

# Identify the differentially produced metabolites 
toptable_pos$Peak_ID <- rownames(toptable_pos)
toptable_pos$Peak_ID <- as.numeric(as.character(toptable_pos$Peak_ID))

toptable_pos <- toptable_pos[,c(ncol(toptable_pos),1:(ncol(toptable_pos)-1))]

toptable_pos <- subset(toptable_pos, toptable_pos$Peak_ID %in% peak_IDs$Peak_ID) # remove peaks without standards
toptable_pos <- merge(toptable_pos, peak_IDs)

# Add the m/z and RT
toptable_pos2 <- inner_join(pos_peak_data, toptable_pos, by= 'Peak_ID')
toptable_pos2<- toptable_pos2[,-c(5:8,10,11,14)]
toptable_pos2 <- with(toptable_pos2, toptable_pos2[order(toptable_pos2$adj.P.Val),])

toptops <- head(toptable_pos2, n=10)

#Testing correct peak data
pos_peak_data_test <- subset(pos_peak_data, pos_peak_data$Peak_ID %in% toptable_pos$Peak_ID)

### Results
# Can use the pos_peak_samples in future experiments. Imputed missing values and normalised via VSN

pos_peaks_samples
write.csv(pos_peaks_samples, file='positive_peak_samples.csv')
