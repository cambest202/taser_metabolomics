### Investigating whether 3 month metabolic changes can be used to predict 18 month outcomes
## May supplement the baseline metabolic signatures of 18 month DAS44 and synovitis outcomes

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
library(pscl)
library(parallel)
library(doParallel)
library(ROCR) 
library(qvalue)

setwd('/Users/cameronbest/University/MRC_DTP_in_Precision_Medicine/Project/RA/Taser_data/Analysis/taser_metabolomics')

### files------
resp_AB <- read.csv(file='20210304_Taser_AB_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
resp_AF <- read.csv(file='20210304_Taser_AF_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
resp_ABF <- read.csv(file='20210304_Taser_ABF_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
peak_IDs <- read.csv(file="20200427_Taser_PeakMetadata.csv", header=TRUE, row.names=1)
peak_metadata <- read.csv(file='20210217_neg_peak_metadata.csv', header=TRUE)
sample_sheet <- read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
novel_peak_data <- read.csv(file='20200430_Taser_NEG_PeakIntensities.csv', header=TRUE, row.names=1)
test <- t(novel_peak_data)
names(peak_IDs)[6] <- 'Peak_ID'
peak_IDs$Peak_ID  <- as.numeric(peak_IDs$Peak_ID )
names(sample_sheet)[2] <- 'Sample_Name'

## Helpful functions
`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)


### Differential peaks across sample time points-----
# Limma function for differential analysis
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
  toptable <- inner_join(toptable, peak_metadata, by = 'Peak_ID')
  }

# Prepare matrices for limma (variables as columns, features to rows)
resp_AB_mat <- as.data.frame((resp_AB))
colnames(resp_AB_mat)[colnames(resp_AB_mat) %like% 'A']<- 'A' ; colnames(resp_AB_mat)[colnames(resp_AB_mat) %like% 'B']<- 'B' 

resp_AF_mat <- as.data.frame((resp_AF))
colnames(resp_AF_mat)[colnames(resp_AF_mat) %like% 'A']<- 'A' ; colnames(resp_AF_mat)[colnames(resp_AF_mat) %like% 'F']<- 'F' 

resp_BF_mat <- as.data.frame((resp_ABF))
resp_BF_mat <- subset(resp_BF_mat, rownames(resp_BF_mat) %notlike% 'A')
colnames(resp_BF_mat)[colnames(resp_BF_mat) %like% 'B']<- 'B' ; colnames(resp_BF_mat)[colnames(resp_BF_mat) %like% 'F']<- 'F' 

# test the same
resp_AB_mat_test <- subset(resp_AB_mat,rownames(resp_AB_mat) %like% 'A')
resp_AF_mat_test <- subset(resp_AF_mat,rownames(resp_AF_mat) %like% 'A')
summary(resp_AB_mat_test$X1 == resp_AF_mat_test$X1)

# format for limma
resp_AB_L <- t(resp_AB_mat)
colnames(resp_AB_L)[colnames(resp_AB_L) %like% 'A'] <- 'A';colnames(resp_AB_L)[colnames(resp_AB_L) %like% 'B'] <- 'B'
rownames(resp_AB_L) <- gsub('X', '', rownames(resp_AB_L))

resp_AF_L <- t(resp_AF_mat)
colnames(resp_AF_L)[colnames(resp_AF_L) %like% 'A'] <- 'A';colnames(resp_AF_L)[colnames(resp_AF_L) %like% 'F'] <- 'F'
rownames(resp_AF_L) <- gsub('X', '', rownames(resp_AF_L))

resp_BF_L <- t(resp_BF_mat)
colnames(resp_BF_L)[colnames(resp_BF_L) %like% 'B'] <- 'B';colnames(resp_BF_L)[colnames(resp_BF_L) %like% 'F'] <- 'F'
rownames(resp_BF_L) <- gsub('X', '', rownames(resp_BF_L))

AB_limma <- limma_fun(resp_AB_L, 1500, 'A', 'B')
AF_limma <- limma_fun(resp_AF_L, 1500, 'A', 'F')
BF_limma <- limma_fun(resp_BF_L, 1500, 'B', 'F')

qvals <- function(limma_table){
  pi0 <- 2*mean(limma_table$P.Value > 0.05)
  lfdrvals <- lfdr(limma_table$P.Value, pi0)
  qobj <- qvalue(limma_table$P.Value)
  hist(qobj)
}
AB_q <- qvals(AB_limma)
AF_q <- qvals(AF_limma)
BF_q <- qvals(BF_limma)
ggarrange(AB_q,AF_q, BF_q,
          labels=c('A', 'B', 'C'))

# Peak Identification Limma
limma_ID <- function(limma_table){
  limma_table_ID <- limma_table
  limma_table_ID$Sig <- 0
  limma_table_ID$Sig <- ifelse(limma_table_ID$adj.P.Val <0.05, '< 0.05', '> 0.05') 
  limma_table_ID$Sig_Peaks <- 0
  limma_table_ID$Sig_Peaks <- ifelse(limma_table_ID$Sig =='< 0.05' & limma_table_ID$identification != '',                                      limma_table_ID$Peak_ID, '')
  limma_table_ID$Sig_Names <-0
  limma_table_ID$Sig_Names <- ifelse(limma_table_ID$Sig =='< 0.05' & limma_table_ID$identification != '',                                      limma_table_ID$Putative_Metabolite, '')
  limma_table_ID%>%
  ggplot(aes(x=logFC, y=-log10(P.Value), 
             colour=Sig, 
             group=Sig)) +
    geom_point (alpha=0.7) +
    theme_minimal() +
    labs (x='LogFC',
          y='-Log p-value',
          colour='Adjusted \np-value',
          title= 'Negative Ion Mode')+
    geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names),
                    box.padding =1,
                    size=2.5,
                    max.overlaps = Inf,
                    position = position_jitter(seed = 1),
                    arrow = arrow(length = unit(0.0015, "npc"))) +  
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))+
    scale_color_brewer(palette = "Set1",direction=-1)
  }
  
AB_ID <- limma_ID(AB_limma)
AB_tab <- AB_ID$data
AF_ID <- limma_ID(AF_limma)
AF_tab <- AF_ID$data
BF_ID <- limma_ID(BF_limma)
BF_tab <- BF_ID$data


### Correlation analysis-----
# Functions ----
cors <- function(melt_table, disease_measure, df_name){
  melt_table$Peak_ID <- as.numeric(melt_table$Peak_ID)
  melt_table <- inner_join(melt_table,peak_metadata, by='Peak_ID')
  ints_nested <- melt_table %>%
    group_by (Peak_ID) %>%
    nest()
  ints_unnested <- melt_table %>%
    unnest(cols=c())
  identical(melt_table, ints_unnested)
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
  best_augmented$adj_p <- p.adjust(best_augmented$p.value, method='BH')
  best_augmented_sig <- subset(best_augmented,best_augmented$adj_p < 0.05)
  best_adj_hmdb <- inner_join(best_augmented_sig, peak_metadata, by='Peak_ID')
  best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
  best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
  gg <-best_adj_hmdb %>%
    subset(Putative_Metabolite != 'Galactonic/Gluconic acid') %>%
    ggplot(aes(x = disease_measure, y=Peak_Intensity)) +
    geom_point(size=1, alpha=0.7) + 
    stat_cor(method = "spearman", 
             vjust=1, hjust=0,
             size=3)+
    geom_smooth(method='lm',
                colour='red')+
    geom_line(aes(y = .fitted), color = "red") +
    facet_wrap(~Putative_Metabolite, scales = "free_y")+
    theme(strip.background = element_rect(fill='white', 
                                          size=1.5),
          strip.text.x= element_text(face = "bold.italic",
                                     size=12))+
    
    labs(x='ΔDAS44',
         y='ΔPeak Intensity',
         title='Negative Ion Mode')+
    theme_minimal()
  print(gg)
  df_name <- return(best_adj_hmdb)
}
sig_dist_peaks <- function(cor_table){
  sig_peaks <- distinct(cor_table, Peak_ID, .keep_all = TRUE)
  sig_peaks <- distinct(sig_peaks, Putative_Metabolite, .keep_all = TRUE)
  sig_peaks <- subset(sig_peaks,sig_peaks$Peak_ID!= 861)
  sig_peaks$ID <- sig_peaks$Peak_ID
  sig_peaks <- sig_peaks[,-c(2,3)]
}

# Prepare datafames
resp_AB_samples <- resp_AB
resp_AB_samples$Sample_Name <- rownames(resp_AB_samples)
#resp_AB_samples$Sample_Name <- substr(resp_AB_samples$Sample_Name, 1,4)

resp_AF_samples <- resp_AF
resp_AF_samples$Sample_Name <- rownames(resp_AF_samples)
#resp_AF_samples$Sample_Name <- substr(resp_AF_samples$Sample_Name, 1,4)

resp_BF_samples <- resp_BF_mat
resp_BF_samples$Sample_Name <- rownames(resp_BF_samples)
#resp_BF_samples$Sample_Name <- substr(resp_BF_samples$Sample_Name, 1,4)

sample_sheet_AB <- subset(sample_sheet, sample_sheet$time == 'A' |sample_sheet$time == 'B' )
sample_sheet_AB$Sample_Name <- paste0(sample_sheet_AB$Sample_Name, sample_sheet_AB$time)
sample_sheet_AF <- subset(sample_sheet, sample_sheet$time == 'A' |sample_sheet$time == 'F' )
sample_sheet_AF$Sample_Name <- paste0(sample_sheet_AF$Sample_Name, sample_sheet_AF$time)
sample_sheet_BF <- subset(sample_sheet, sample_sheet$time == 'B' |sample_sheet$time == 'F' )
sample_sheet_BF$Sample_Name <- paste0(sample_sheet_BF$Sample_Name, sample_sheet_BF$time)

# Complete sample_sheet and complete peak intensities
resp_int_AB <- inner_join(sample_sheet_AB, resp_AB_samples, by='Sample_Name')
resp_int_AB$Sample <- substr(resp_int_AB$Sample_Name, 1,4)
resp_int_AB <- resp_int_AB[,c(ncol(resp_int_AB),1:(ncol(resp_int_AB)-1))]

resp_int_AF <- inner_join(sample_sheet_AF, resp_AF_samples, by='Sample_Name')
resp_int_AF$Sample <- substr(resp_int_AF$Sample_Name, 1,4)
resp_int_AF <- resp_int_AF[,c(ncol(resp_int_AF),1:(ncol(resp_int_AF)-1))]

resp_int_BF <- inner_join(sample_sheet_BF, resp_BF_samples, by='Sample_Name')
resp_int_BF$Sample <- substr(resp_int_BF$Sample_Name, 1,4)
resp_int_BF <- resp_int_BF[,c(ncol(resp_int_BF),1:(ncol(resp_int_BF)-1))]

# Melt the ints samples with peak intensities then add the disease data after
# AB ints
resp_int_AB_2 <- subset(resp_int_AB, resp_int_AB$Sample %in% resp_diff_AB$Sample) ## need to first produce resp_diff_AB which is done below.. 

resp_int_A <- subset(resp_int_AB_2,resp_int_AB_2$time =='A')
resp_int_B <- subset(resp_int_AB_2,resp_int_AB_2$time =='B')

resp_ints_AB_melt <- resp_int_A[,c(1,13:1470)]
resp_ints_AB_melt <- melt(resp_ints_AB_melt)
names(resp_ints_AB_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_ints_AB_melt$End_DAS44 <- resp_int_B$DAS44
resp_ints_AB_melt$ΔDAS44 <- resp_diff_AB$DAS44
resp_ints_AB_melt$Peak_ID <- gsub('X','', resp_ints_AB_melt$Peak_ID)
resp_ints_AB_melt[, c(2:5)] <- map_df(resp_ints_AB_melt[, c(2:5)], as.numeric)

#AF
resp_int_AF_2 <- subset(resp_int_AF, resp_int_AF$Sample %in% resp_diff_AF$Sample) ## need to first produce resp_diff_AB which is done below.. 

resp_int_A_2 <- subset(resp_int_AF_2,resp_int_AF_2$time =='A')
resp_int_F <- subset(resp_int_AF_2,resp_int_AF_2$time =='F')

resp_int_A_2 <- subset(resp_int_A_2,resp_int_A_2$Sample %in% resp_int_F$Sample)

resp_ints_AF_melt <- resp_int_A_2[,c(1,13:1470)]
resp_ints_AF_melt <- melt(resp_ints_AF_melt)
names(resp_ints_AF_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_ints_AF_melt$End_DAS44 <- resp_int_F$DAS44
resp_ints_AF_melt$ΔDAS44 <- resp_diff_AF$DAS44
resp_ints_AF_melt$Peak_ID <- gsub('X','', resp_ints_AF_melt$Peak_ID)
resp_ints_AF_melt[, c(2:5)] <- map_df(resp_ints_AF_melt[, c(2:5)], as.numeric)

#BF
resp_int_BF_2 <- subset(resp_int_BF, resp_int_BF$Sample %in% resp_diff_BF$Sample) ## need to first produce resp_diff_AB which is done below.. 

resp_int_B <- subset(resp_int_BF_2,resp_int_BF_2$time =='B')
resp_int_F_2 <- subset(resp_int_BF_2,resp_int_BF_2$time =='F')

resp_int_B <- subset(resp_int_B,resp_int_B$Sample %in% resp_int_F_2$Sample)

resp_ints_BF_melt <- resp_int_B[,c(1,13:1470)]
resp_ints_BF_melt <- melt(resp_ints_BF_melt)
names(resp_ints_BF_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_ints_BF_melt$End_DAS44 <- resp_int_F_2$DAS44
resp_ints_BF_melt$ΔDAS44 <- resp_diff_BF$DAS44
resp_ints_BF_melt$Peak_ID <- gsub('X','', resp_ints_BF_melt$Peak_ID)
resp_ints_BF_melt[, c(2:5)] <- map_df(resp_ints_BF_melt[, c(2:5)], as.numeric)

### Correlations between baselie peak intensity  and DAS44 changes------
AB_cors_base <- cors(resp_ints_AB_melt, ΔDAS44)
AB_sig_base <- sig_dist_peaks(AB_cors_base) ## no significant correlations

AF_cors_base <- cors(resp_ints_AF_melt, ΔDAS44)
AF_sig_base <- sig_dist_peaks(AF_cors_base) ## no significant correlations

BF_cors_base <- cors(resp_ints_BF_melt, ΔDAS44)
BF_sig_base <- sig_dist_peaks(BF_cors_base) ## no significant correlations

# Differential sample sheet and peak intensities by samples -----
resp_diff_AB <- resp_int_AB[,c(1,4:11,13:1470)]
resp_diff_AB <- aggregate(.~Sample, resp_diff_AB, diff, na.rm=TRUE)
resp_diff_AB <- resp_diff_AB[resp_diff_AB$CRP != 'numeric(0)', ]

resp_diff_AF <- resp_int_AF[,c(1,4:11,13:1470)]
resp_diff_AF <- aggregate(.~Sample, resp_diff_AF, diff, na.rm=TRUE)
resp_diff_AF <- resp_diff_AF[resp_diff_AF$CRP != 'numeric(0)', ]

resp_diff_BF <- resp_int_BF[,c(1,4:11,13:1470)]
resp_diff_BF <- aggregate(.~Sample, resp_diff_BF, diff, na.rm=TRUE)
resp_diff_BF <- resp_diff_BF[resp_diff_BF$CRP != 'numeric(0)', ]

# Double the differential dfs to use alongside the complete peak df
resp_dbl_diff_AB <-rbind.data.frame(resp_diff_AB, resp_diff_AB)
resp_dbl_diff_AB <- with(resp_dbl_diff_AB,resp_dbl_diff_AB[order(Sample),])
resp_dbl_diff_AB[, c(2:1467)] <- map_df(resp_dbl_diff_AB[, c(2:1467)], as.numeric)
resp_dbl_diff_AF <-rbind.data.frame(resp_diff_AF, resp_diff_AF)
resp_dbl_diff_AF <- with(resp_dbl_diff_AF,resp_dbl_diff_AF[order(Sample),])
resp_dbl_diff_AF[, c(2:1467)] <- map_df(resp_dbl_diff_AF[, c(2:1467)], as.numeric)
resp_dbl_diff_BF <-rbind.data.frame(resp_diff_BF, resp_diff_BF)
resp_dbl_diff_BF <- with(resp_dbl_diff_BF,resp_dbl_diff_BF[order(Sample),])
resp_dbl_diff_BF[, c(2:1467)] <- map_df(resp_dbl_diff_BF[, c(2:1467)], as.numeric)

# Melt the samples with peak intensities and add sample sheets after ----
resp_diff_AB_melt <- resp_dbl_diff_AB[,c(1,10:1467)]
resp_diff_AB_melt <- melt(resp_diff_AB_melt)
names(resp_diff_AB_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_diff_AB_melt$DAS44 <- resp_dbl_diff_AB$DAS44
resp_diff_AB_melt$Peak_ID <- gsub('X','', resp_diff_AB_melt$Peak_ID)
resp_diff_AB_melt$CRP <- resp_dbl_diff_AB$CRP
resp_diff_AB_melt$ESR <- resp_dbl_diff_AB$ESR
resp_diff_AB_melt$HAQ <- resp_dbl_diff_AB$HAQ
resp_diff_AB_melt$GHVAS <- resp_dbl_diff_AB$GHVAS
resp_diff_AB_melt$PVAS <- resp_dbl_diff_AB$PVAS
resp_diff_AB_melt$RAI <- resp_dbl_diff_AB$RAI
resp_diff_AB_melt$SJC <- resp_dbl_diff_AB$SJC

resp_diff_AF_melt <- resp_dbl_diff_AF[,c(1,10:1467)]
resp_diff_AF_melt <- melt(resp_diff_AF_melt)
names(resp_diff_AF_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_diff_AF_melt$DAS44 <- resp_dbl_diff_AF$DAS44
resp_diff_AF_melt$Peak_ID <- gsub('X','', resp_diff_AF_melt$Peak_ID)
resp_diff_AF_melt$CRP <- resp_dbl_diff_AF$CRP
resp_diff_AF_melt$ESR <- resp_dbl_diff_AF$ESR
resp_diff_AF_melt$HAQ <- resp_dbl_diff_AF$HAQ
resp_diff_AF_melt$GHVAS <- resp_dbl_diff_AF$GHVAS
resp_diff_AF_melt$PVAS <- resp_dbl_diff_AF$PVAS
resp_diff_AF_melt$RAI <- resp_dbl_diff_AF$RAI
resp_diff_AF_melt$SJC <- resp_dbl_diff_AF$SJC

resp_diff_BF_melt <- resp_dbl_diff_BF[,c(1,10:1467)]
resp_diff_BF_melt <- melt(resp_diff_BF_melt)
names(resp_diff_BF_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_diff_BF_melt$DAS44 <- resp_dbl_diff_BF$DAS44
resp_diff_BF_melt$Peak_ID <- gsub('X','', resp_diff_BF_melt$Peak_ID)
resp_diff_BF_melt$CRP <- resp_dbl_diff_BF$CRP
resp_diff_BF_melt$ESR <- resp_dbl_diff_BF$ESR
resp_diff_BF_melt$HAQ <- resp_dbl_diff_BF$HAQ
resp_diff_BF_melt$GHVAS <- resp_dbl_diff_BF$GHVAS
resp_diff_BF_melt$PVAS <- resp_dbl_diff_BF$PVAS
resp_diff_BF_melt$RAI <- resp_dbl_diff_BF$RAI
resp_diff_BF_melt$SJC <- resp_dbl_diff_BF$SJC

### Correlations between peak intensity changes and DAS44 changes------
AB_cors <- cors(resp_diff_AB_melt, DAS44)
AB_sig <- sig_dist_peaks(AB_cors)

AF_cors <- cors(resp_diff_AF_melt, DAS44)
AF_sig <- sig_dist_peaks(AF_cors)

BF_cors <- cors(resp_diff_BF_melt, DAS44)
BF_sig <- sig_dist_peaks(BF_cors)

## Correlating 3 month metabolite change with 18 month disease change
resp_diff_AB_AFdas <- resp_diff_AB
resp_diff_AB_AFdas <- subset(resp_diff_AB_AFdas, resp_diff_AB_AFdas$Sample %in% resp_diff_AF$Sample)
resp_diff_AF_sub <- subset(resp_diff_AF, resp_diff_AF$Sample %in% resp_diff_AB_AFdas$Sample)

resp_diff_AB_AFdas_melt <- resp_diff_AB_AFdas[,c(1,10:1467)]
resp_diff_AB_AFdas_melt[, c(2:1459)] <- map_df(resp_diff_AB_AFdas_melt[, c(2:1459)], as.numeric)
resp_diff_AB_AFdas_melt <- melt(resp_diff_AB_AFdas_melt)
names(resp_diff_AB_AFdas_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
resp_diff_AB_AFdas_melt$DAS44 <- resp_diff_AF_sub$DAS44
resp_diff_AB_AFdas_melt$DAS44 <- as.numeric(resp_diff_AB_AFdas_melt$DAS44)
resp_diff_AB_AFdas_melt$Peak_ID <- gsub('X','', resp_diff_AB_AFdas_melt$Peak_ID)
resp_diff_AB_AFdas_melt$CRP <- resp_diff_AF_sub$CRP
resp_diff_AB_AFdas_melt$ESR <- resp_diff_AF_sub$ESR
resp_diff_AB_AFdas_melt$HAQ <- resp_diff_AF_sub$HAQ
resp_diff_AB_AFdas_melt$GHVAS <- resp_diff_AF_sub$GHVAS
resp_diff_AB_AFdas_melt$PVAS <- resp_diff_AF_sub$PVAS
resp_diff_AB_AFdas_melt$RAI <- resp_diff_AF_sub$RAI
resp_diff_AB_AFdas_melt$SJC <- resp_diff_AF_sub$SJC

AB_AFdiff_cors <- cors(resp_diff_AB_AFdas_melt, DAS44) # no metabolites correlated here, so no output
AB_AFdiff_sig <- sig_dist_peaks(AB_AFdiff_cors)

### Logistic regression model to predict 18 month DAS44 change
# Use AB_sig peak list to select metabolites to use in the LRM for predicting outcomes
# Outcomes based on EULAR recommendations: ΔDAS44 reduction > 1.2 // End DAS44 < 2.4 is positive result
sample_AF_outcome <- sample_sheet_AF[,c(2,4)]
sample_AF_outcome$Sample_Name <- substr(sample_AF_outcome$Sample_Name, 1,4)
sample_F_outcome <- subset(sample_AF_outcome , rownames(sample_AF_outcome) %like% 'F')
sample_AF_diff <- aggregate(.~Sample_Name, sample_AF_outcome, diff, na.rm= TRUE)
sample_AF_diff <- subset(sample_AF_diff,sample_AF_diff$DAS44 != 'numeric(0)')
sample_AF_diff$DAS44<- as.numeric(sample_AF_diff$DAS44)
sample_AF_diff$End_DAS44 <- 0
names(sample_AF_diff)[1:2] <- c('Sample_Name', 'ΔDAS44')

sample_F_outcome_cut <- subset(sample_F_outcome,sample_F_outcome$Sample_Name %in% sample_AF_diff$Sample_Name)

sample_AF_diff$End_DAS44 <- sample_F_outcome_cut$DAS44
sample_AF_diff$Response <- 0
sample_AF_diff$Response <- with(sample_AF_diff, 
                                ifelse(End_DAS44 <= 1.6, 'Remission',
                                       ifelse(End_DAS44< 2.4 | ΔDAS44 < -1.2, 'Good',
                                              ifelse(End_DAS44 > 3.7 | ΔDAS44 < -0.6, 'Poor', 'Mild'))))

# May need to make the classification of a positive response more strict. 
# Use the end DAS44 < 2.4 as the indication of a good response
sample_AF_diff$Response <- as.factor(sample_AF_diff$Response)
#sample_AF_diff$Response <- with(sample_AF_diff, ifelse(End_DAS44 < 2.4, 1, 0))
sample_AF_diff$Response_Sim <- with(sample_AF_diff, ifelse(ΔDAS44 < -1.2, 1, 0))
sample_AF_diff$Response_Sim <- ifelse(sample_AF_diff$Response =='Remission', 1,0)
sample_AF_diff$Response_Sim <- as.factor(sample_AF_diff$Response_Sim)

summary(sample_AF_diff$Response_Sim)

# Prepare AB peak table
AB_peak_change <- resp_diff_AB[,-c(2:9)]
AB_peak_change[, c(2:1459)] <- map_df(AB_peak_change[, c(2:1459)], as.numeric)
names(AB_peak_change)[1] <-'Sample_Name'
ABF_tab <- inner_join(sample_AF_diff, AB_peak_change, by='Sample_Name')
ABF_prelrm <- ABF_tab[,-c(1:4)]
names(ABF_prelrm)[1] <-'Response'

# Logistic Regression
# Peaks to address: 
# AB mets to predict AF das44 (ABF_prelrm)----
# model_1 X78 + X189 + X23 + X26 +  X274 +  X642 + X421 + X19 + X1062 + X524 + X415 + X162 + X738 +  X555 + X48 + X216
#model_2 X1 + X961 + X787 + X902 + X582 + X770 + X70 + X232 +   X750 + X1165 + X778, 
#model <- glm(Response ~ X1+  X902  + X74+ X270+ X421+ X555+X770   +   X750 + X1165 + X778,family = binomial(link = "logit"),  data = train_data)

# AB mets to predict AB outcomes------

sig_sigs <- subset(sig_peaks, sig_peaks$adj_p < 0.01)
peaks <- sig_sigs[,c(26,32)]
peaks_list <- paste0('X',peaks$ID)

set.seed(42)
index <- createDataPartition(ABF_prelrm$Response, p = 0.7, list = FALSE)
train_data <- ABF_prelrm[index, ]
test_data  <- ABF_prelrm[-index, ]

summary(model)
anova(model, test="Chisq")
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


