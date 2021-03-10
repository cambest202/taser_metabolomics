### Positive Ion Mode 
###Investigating whether 3 month metabolic changes can be used to predict 18 month outcomes
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
library(corrr)
library(ggcorrplot)
library(ape)
library(forcats)


setwd('/Users/cameronbest/University/MRC_DTP_in_Precision_Medicine/Project/RA/Taser_data/Analysis/taser_metabolomics')

### files------
resp_AB <- read.csv(file='20210308_Taser_AB_POS_PeakIntensities.csv', header=TRUE, row.names=1)
resp_AF <- read.csv(file='20210308_Taser_AF_POS_PeakIntensities.csv', header=TRUE, row.names=1)
resp_ABF <- read.csv(file='20210308_Taser_ABF_POS_PeakIntensities.csv', header=TRUE, row.names=1)
peak_metadata <- read.csv(file='20200521_Taser_POS_Peakdata.csv', header=TRUE, row.names=1)
sample_sheet <- read.table (file="20200318_Taser_SampleSheet.csv", header=TRUE, row.names=1)
patient_metadata <- read.csv(file='20190713_Taser_PatientMetadata.csv', header=TRUE, row.names=1)
kegg <- read.csv('20210305_kegg_list.csv')

kegg <- distinct(kegg, Kegg_code, .keep_all = TRUE)
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

AB_limma <- limma_fun(resp_AB_L, 1600, 'A', 'B')
AF_limma <- limma_fun(resp_AF_L, 1600, 'A', 'F')
BF_limma <- limma_fun(resp_BF_L, 1600, 'B', 'F')

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
    subset(Putative_Metabolite != 'Sulfasalazine')%>%
    subset(Putative_Metabolite != 'Galactonic/Gluconic acid')%>%
    ggplot(aes(x=logFC, y=-log10(P.Value), 
               colour=Sig, 
               group=Sig)) +
    geom_point (alpha=0.7) +
    theme_minimal() +
    labs (x='LogFC',
          y='-Log p-value',
          colour='Adjusted \np-value',
          title= 'Positive Ion Mode')+
    geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names),
                    box.padding =1,
                    size=2.5,
                    max.overlaps = Inf,
                    position = position_jitter(seed = 1),
                    arrow = arrow(length = unit(0.0015, "npc"))) +  
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))+
    scale_color_brewer(palette = "Set1",direction=-1)+
    ylim(0,15)
}

AB_ID <- limma_ID(AB_limma)
AB_tab <- AB_ID$data
AF_ID <- limma_ID(AF_limma)
AF_tab <- AF_ID$data
BF_ID <- limma_ID(BF_limma)
BF_tab <- BF_ID$data

ggarrange(AB_ID,AF_ID, BF_ID,
          labels=c('A', 'B', 'C'),
          ncol=3)

### Correlation analysis-----
# Functions ----
cors <- function(melt_table, disease_measure){
  melt_table$Peak_ID <- as.numeric(melt_table$Peak_ID)
  melt_table <- inner_join(melt_table,peak_metadata, by='Peak_ID')
  ints_nested <- melt_table %>%
    group_by (Peak_ID) %>%
    nest()
  ints_unnested <- melt_table %>%
    unnest(cols=c())
  identical(melt_table, ints_unnested)
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
  best_augmented$adj_p <- p.adjust(best_augmented$p.value, method='BH')
  best_augmented_sig <- subset(best_augmented,best_augmented$adj_p < 0.05)
  best_adj_hmdb <- inner_join(best_augmented_sig, peak_metadata, by='Peak_ID')
  best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
  best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
  gg <- best_adj_hmdb %>%
    subset(Putative_Metabolite != 'Galactonic/Gluconic acid') %>%
    subset(Putative_Metabolite != 'GPC') %>%
    subset(Putative_Metabolite != 'LysoPC') %>%
    subset(Peak_ID != '1460') %>%
    ggplot(aes(x = DAS44, y=Peak_Intensity)) +
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
         title='Positive Ion Mode')+
    theme_minimal()
  print(gg)
  df_name <- return(best_adj_hmdb)
}
cors_len <- function(melt_table, disease_measure, df_name){
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
  best_augmented_sig <- subset(best_augmented,best_augmented$adj_p < 0.6)
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

# Prepare datafames -------
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

# Melt the ints samples with peak intensities then add the disease data after -------
# But first do the resp_diff as need this for the subsetting
# Differential sample sheet and peak intensities by samples -----
resp_diff_AB <- resp_int_AB[,c(1,4:11,13:1596)]
resp_diff_AB <- aggregate(.~Sample, resp_diff_AB, diff, na.rm=TRUE)
resp_diff_AB <- resp_diff_AB[resp_diff_AB$CRP != 'numeric(0)', ]

resp_diff_AF <- resp_int_AF[,c(1,4:11,13:1596)]
resp_diff_AF <- aggregate(.~Sample, resp_diff_AF, diff, na.rm=TRUE)
resp_diff_AF <- resp_diff_AF[resp_diff_AF$CRP != 'numeric(0)', ]

resp_diff_BF <- resp_int_BF[,c(1,4:11,13:1596)]
resp_diff_BF <- aggregate(.~Sample, resp_diff_BF, diff, na.rm=TRUE)
resp_diff_BF <- resp_diff_BF[resp_diff_BF$CRP != 'numeric(0)', ]

# Double the differential dfs to use alongside the complete peak df ----------
resp_dbl_diff_AB <-rbind.data.frame(resp_diff_AB, resp_diff_AB)
resp_dbl_diff_AB <- with(resp_dbl_diff_AB,resp_dbl_diff_AB[order(Sample),])
resp_dbl_diff_AB[, c(2:1593)] <- map_df(resp_dbl_diff_AB[, c(2:1467)], as.numeric)
resp_dbl_diff_AF <-rbind.data.frame(resp_diff_AF, resp_diff_AF)
resp_dbl_diff_AF <- with(resp_dbl_diff_AF,resp_dbl_diff_AF[order(Sample),])
resp_dbl_diff_AF[, c(2:1593)] <- map_df(resp_dbl_diff_AF[, c(2:1467)], as.numeric)
resp_dbl_diff_BF <-rbind.data.frame(resp_diff_BF, resp_diff_BF)
resp_dbl_diff_BF <- with(resp_dbl_diff_BF,resp_dbl_diff_BF[order(Sample),])
resp_dbl_diff_BF[, c(2:1593)] <- map_df(resp_dbl_diff_BF[, c(2:1467)], as.numeric)

# AB ints
resp_int_AB_2 <- subset(resp_int_AB, resp_int_AB$Sample %in% resp_diff_AB$Sample) ## need to first produce resp_diff_AB which is done below.. 

resp_int_A <- subset(resp_int_AB_2,resp_int_AB_2$time =='A')
resp_int_B <- subset(resp_int_AB_2,resp_int_AB_2$time =='B')

resp_ints_AB_melt <- resp_int_A[,c(1,13:1596)]
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

resp_ints_AF_melt <- resp_int_A_2[,c(1,13:1596)]
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

resp_ints_BF_melt <- resp_int_B[,c(1,13:1596)]
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

# use more lenient model (cors_len)
AB_cors_base <- cors_len(resp_ints_AB_melt, ΔDAS44)
AB_sig_base <- sig_dist_peaks(AB_cors_base) ## no significant correlations

AF_cors_base <- cors_len(resp_ints_AF_melt, ΔDAS44) #---------
resp_ints_AF_melt$Peak_ID <- as.numeric(resp_ints_AF_melt$Peak_ID)
resp_ints_AF_melt <- inner_join(resp_ints_AF_melt,peak_metadata, by='Peak_ID')
ints_nested <- resp_ints_AF_melt %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_ints_AF_melt %>%
  unnest(cols=c())
identical(resp_ints_AF_melt, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = map(data, ~lm(formula = Peak_Intensity~ΔDAS44, data = .x)))
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
best_augmented_sig <- subset(best_augmented,best_augmented$adj_p < 0.6)
best_adj_hmdb <- inner_join(best_augmented_sig, peak_metadata, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')

best_adj_hmdb %>%
  subset(Putative_Metabolite != 'Galactonic/Gluconic acid') %>%
  subset(Putative_Metabolite != '') %>%
  ggplot(aes(x = ΔDAS44, y=Peak_Intensity)) +
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
  theme_minimal()#---------- #nothing of note here, just including method



# Melt the samples with peak intensities and add sample sheets after ----
resp_diff_AB_melt <- resp_dbl_diff_AB[,c(1,10:1593)]
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

resp_diff_AF_melt <- resp_dbl_diff_AF[,c(1,10:1593)]
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

resp_diff_BF_melt <- resp_dbl_diff_BF[,c(1,10:1593)]
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
# 3 month metabolic change correlating with 3 month DAS44 change
AB_cors <- cors(resp_diff_AB_melt)
AB_sig <- sig_dist_peaks(AB_cors)

# 18 month metabolic change correlating with 18 month DAS44 change
AF_cors <- cors(resp_diff_AF_melt)
AF_sig <- sig_dist_peaks(AF_cors)

# 15 month metabolic change correlating with 15 month DAS44 change
BF_cors <- cors(resp_diff_BF_melt)
BF_sig <- sig_dist_peaks(BF_cors)
BF <- BF_sig[,c(1,26)]

## Correlating 3 month metabolite change with 18 month disease change
resp_diff_AB_AFdas <- resp_diff_AB
resp_diff_AB_AFdas <- subset(resp_diff_AB_AFdas, resp_diff_AB_AFdas$Sample %in% resp_diff_AF$Sample)
resp_diff_AF_sub <- subset(resp_diff_AF, resp_diff_AF$Sample %in% resp_diff_AB_AFdas$Sample)

resp_diff_AB_AFdas_melt <- resp_diff_AB_AFdas[,c(1,10:1593)]
resp_diff_AB_AFdas_melt[, c(2:1585)] <- map_df(resp_diff_AB_AFdas_melt[, c(2:1585)], as.numeric)
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

# 3 month metabolic change correlating with 18 month DAS44 change**** -------
AB_AFdiff_cors <- cors(resp_diff_AB_AFdas_melt, DAS44) # no metabolites correlated here, so no output
AB_AFdiff_sig <- sig_dist_peaks(AB_AFdiff_cors)

# Using a more lenient criteria to obtain possible metabolic changes associated to test in LRM later  --------
resp_diff_AB_AFdas_melt$Peak_ID <- as.numeric(resp_diff_AB_AFdas_melt$Peak_ID)
resp_diff_AB_AFdas_melt <- inner_join(resp_diff_AB_AFdas_melt,peak_metadata, by='Peak_ID')
ints_nested <- resp_diff_AB_AFdas_melt %>%
  group_by (Peak_ID) %>%
  nest()
ints_unnested <- resp_diff_AB_AFdas_melt %>%
  unnest(cols=c())
identical(resp_diff_AB_AFdas_melt, ints_unnested)
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
best_augmented$adj_p <- p.adjust(best_augmented$p.value, method='BH')
best_augmented_sig <- subset(best_augmented,best_augmented$adj_p < 0.1)
best_adj_hmdb <- inner_join(best_augmented_sig, peak_metadata, by='Peak_ID')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
best_dist <- distinct(best_adj_hmdb, Putative_Metabolite, .keep_all=TRUE)
best_dist <- best_dist[,-c(2,3)]

best_adj_hmdb %>%
  subset(Putative_Metabolite != 'Galactonic/Gluconic acid') %>%
  subset(Peak_ID != '861')%>%
  subset(Peak_ID != '101')%>%
  subset(Peak_ID != '415')%>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
  geom_point(size=1, alpha=0.7) + 
  theme_minimal()+
  stat_cor(method = "spearman", 
           vjust=1, hjust=0,
           size=3)+
  geom_smooth(method='lm',
              colour='red')+
  geom_line(aes(y = .fitted), color = "red") +
  facet_wrap(~Putative_Metabolite, scales = "free_y")+
  theme(strip.background = element_rect(fill='white',
                                        colour='white',
                                        size=1.5),
        strip.text.x= element_text(#face = "bold",
          size=12))+
  labs(x='ΔDAS44',
       y='ΔPeak Intensity',
       title='Negative Ion Mode')

### Logistic regression model to predict 18 month DAS44 change
# Use AB_sig peak list to select metabolites to use in the LRM for predicting outcomes
# Outcomes based on EULAR recommendations: ΔDAS44 reduction > 1.2 // End DAS44 < 2.4 is positive result
#AB Sample Sheet -----
sample_AB_outcome <- sample_sheet_AB[,c(2,4)]
sample_AB_outcome$Sample_Name <- substr(sample_AB_outcome$Sample_Name, 1,4)
sample_B_outcome <- subset(sample_AB_outcome , rownames(sample_AB_outcome) %like% 'B')
sample_AB_diff <- aggregate(.~Sample_Name, sample_AB_outcome, diff, na.rm= TRUE)
sample_AB_diff <- subset(sample_AB_diff,sample_AB_diff$DAS44 != 'numeric(0)')
sample_AB_diff$DAS44<- as.numeric(sample_AB_diff$DAS44)
sample_AB_diff$End_DAS44 <- 0
names(sample_AB_diff)[1:2] <- c('Sample_Name', 'ΔDAS44')

sample_B_outcome_cut <- subset(sample_B_outcome,sample_B_outcome$Sample_Name %in% sample_AB_diff$Sample_Name)

sample_AB_diff$End_DAS44 <- sample_B_outcome_cut$DAS44
sample_AB_diff$Response <- 0
sample_AB_diff$Response <- with(sample_AB_diff, 
                                ifelse(End_DAS44 <= 1.6, 'Remission',
                                       ifelse(End_DAS44< 2.4 | ΔDAS44 < -1.2, 'Good',
                                              ifelse(End_DAS44 > 3.7 | ΔDAS44 < -0.6, 'Poor', 'Mild'))))
sample_AB_diff$Response <- ifelse(sample_AB_diff$Response == 'Remission' | sample_AB_diff$Response == 'Good', 1,0)

# Prepare AB peak table
AB_peak_change <- resp_diff_AF[,-c(2:9)]
AB_peak_change[, c(2:1585)] <- map_df(AB_peak_change[, c(2:1585)], as.numeric)
names(AB_peak_change)[1] <-'Sample_Name'
AB_tab <- inner_join(sample_AB_diff, AB_peak_change, by='Sample_Name')
AB_prelrm <- AB_tab[,-c(1:3)]
names(AB_prelrm$Response)[1] <-'Response'
AB_prelrm$Response <- as.factor(AB_prelrm$Response)

#AF Sample Sheet -----
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
sample_AF_diff$Response <- ifelse(sample_AF_diff$Response == 'Remission' | sample_AF_diff$Response == 'Good', 1,0)

# Prepare AF peak tAFle
AF_peak_change <- resp_diff_AF[,-c(2:9)]
AF_peak_change[, c(2:1585)] <- map_df(AF_peak_change[, c(2:1585)], as.numeric)
names(AF_peak_change)[1] <-'Sample_Name'
AF_tAF <- inner_join(sample_AF_diff, AF_peak_change, by='Sample_Name')
AF_prelrm <- AF_tAF[,-c(1:3)]
names(AF_prelrm$Response)[1] <-'Response'
AF_prelrm$Response <- as.factor(AF_prelrm$Response)

#ABF Sample Sheet -----
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

# Use the end DAS44 < 2.4 as the indication of a good response
sample_AF_diff$Response <- as.factor(sample_AF_diff$Response)
#sample_AF_diff$Response <- with(sample_AF_diff, ifelse(End_DAS44 < 2.4, 1, 0))
sample_AF_diff$Response_Sim <- with(sample_AF_diff, ifelse(ΔDAS44 < -1.2, 1, 0))
sample_AF_diff$Response_Sim <- ifelse(sample_AF_diff$Response =='Remission', 1,0)
sample_AF_diff$Response_Sim <- as.factor(sample_AF_diff$Response_Sim)

summary(sample_AF_diff$Response_Sim)

# Prepare AB peak table
AB_peak_change <- resp_diff_AB[,-c(2:9)]
AB_peak_change[, c(2:1585)] <- map_df(AB_peak_change[, c(2:1585)], as.numeric)
names(AB_peak_change)[1] <-'Sample_Name'
ABF_tab <- inner_join(sample_AF_diff, AB_peak_change, by='Sample_Name')
patient_metadata$Sample_Name <- rownames(patient_metadata)
ABF_tab_2 <- inner_join(patient_metadata,ABF_tab, by='Sample_Name')
#ABF_prelrm <- ABF_tab[,-c(1:4)]
ABF_prelrm <-ABF_tab_2
#names(ABF_prelrm)[14] <-'Response'
ABF_prelrm$Response <- ifelse(ABF_prelrm$Response == 'Remission' | ABF_prelrm$Response == 'Good', 1,0)
ABF_prelrm$Response <- as.factor(ABF_prelrm$Response)

# Prepare baseline metabolite


# Differential sample_sheet_AB and baseline sample_sheet_A ---------
sample_sheet_A <- subset(sample_sheet_AB, sample_sheet_AB$time =='A')
sample_sheet_A$Sample_Name <- substr(sample_sheet_A$Sample_Name,1,4)

sample_sheet_AB$Sample_Name <- substr(sample_sheet_AB$Sample_Name,1,4)
sample_sheet_AB_diff <- aggregate(.~Sample_Name, sample_sheet_AB[c(2:10)], diff, na.rm=TRUE)
sample_sheet_AB_diff <- subset(sample_sheet_AB_diff,sample_sheet_AB_diff$CRP !='numeric(0)')
sample_sheet_AB_diff[, c(2:9)] <- map_df(sample_sheet_AB_diff[, c(2:9)], as.numeric)
names(sample_sheet_AB_diff)[2:9]<- paste0('Δ', names(sample_sheet_AB_diff)[2:9])

sample_sheet_AB_compl <- inner_join(sample_sheet_AB_diff,sample_sheet_A, by='Sample_Name')

# Join 3 month Δpeak matrix, patient metadata and adjusted sample_sheet with 3 month and baseline disease measures -------
ABF_prelrm <- inner_join(ABF_prelrm,sample_sheet_AB_compl, by='Sample_Name')

# AB mets to predict AF das44 (ABF_prelrm)----
# model_1 X78 + X189 + X23 + X26 +  X274 +  X642 + X421 + X19 + X1062 + X524 + X415 + X162 + X738 +  X555 + X48 + X216
#model_2 X1 + X961 + X787 + X902 + X582 + X770 + X70 + X232 +   X750 + X1165 + X778, 
#model <- glm(Response ~ X1+  X902  + X74+ X270+ X421+ X555+X770   +   X750 + X1165 + X778,family = binomial(link = "logit"),  data = train_data)
#model <- glm(Response ~ X78  + X26 +  X274 +  X642 + X524 + X415 + X162 + X738 +  X555 + X48 + X216,family = binomial(link = "logit"),  data = train_data)
# best model yet: model <- glm(Response ~ X421 + X26 + X162 +  X555 + X415 + X642 + X216,family = binomial(link = "logit"),  data = train_data)
# AUC 0.545 X738 +  + X642 + X274 +  X216 + X432 +  X70 + X115 + X19,
sig_sigs <- best_dist
peaks_list <- paste0('X',sig_sigs$Peak_ID)
metadata_only <- ABF_prelrm[,-c(15:1472)]
# AF mets to predict AF outcomes------
#AF_prelrm- baseline to 18 month matrix with response
# AB_prelrm- baseline to 3 month matrix with response
# ABF_prelrm- 3 month metabolic change matrix with 18 month response
AF_peaks_cor <- AF_sig$Peak_ID
lrm_output <- function(model){
  fitted.results <- predict(model,newdata=test_data,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misClasificError <- mean(fitted.results != test_data$Response)
  p <- predict(model,newdata=test_data,type='response')
  pr <- prediction(p, test_data$Response)
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  print(summary(model))
  print(paste('Accuracy',1-misClasificError))
  print(auc)
  print(plot(prf))
}
set.seed(42)
index <- createDataPartition(AF_prelrm$Response, p = 0.65, list = FALSE)
train_data <- AF_prelrm[index, ]
test_data  <- AF_prelrm[-index, ]

model <- glm(Response ~ X274 +  X216 + X1062 + X162 ,
             family = binomial(link = "logit"),  data = train_data)

lrm_output(model)

LRM_peaks <- c(274, 216, 1062, 162)
LRM_AB_diff_peaks <- subset(peaks, peaks$ID %in% LRM_peaks)


AF_cors_select <- subset(AF_cors, AF_cors$Peak_ID %in% LRM_peaks)
AF_cors_select %>%
  ggplot(aes(x = DAS44, y=Peak_Intensity)) +
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
       y='ΔPeak Intensity')+
  theme_minimal()




## Plotting the course of metabolite levels during treatment (baseline -- > 3 months --> 18 months) -----
# Look specifically at the metabolites which significantly correlated with DAS44 changes over 18 months
time_course <- resp_ABF
time_course$Sample_Name <- rownames(time_course)

time_melt <- melt(time_course)
names(time_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
time_melt$Time <-0
time_melt$Time[time_melt$Sample_Name %like% 'A'] <- 0
time_melt$Time[time_melt$Sample_Name %like% 'B'] <- 3
time_melt$Time[time_melt$Sample_Name %like% 'F'] <- 18

time_melt$Sample_Name <- substr(time_melt$Sample_Name,1,4)

time_melt_mean <- aggregate(time_melt$Peak_Intensity, FUN=mean, 
                            by=list(time_melt$Peak_ID, time_melt$Time))
names(time_melt_mean) <- c('Peak_ID', 'Time', 'Mean_Peak_Intensity')
time_melt_mean$Peak_ID <- gsub('X', '', time_melt_mean$Peak_ID)
time_melt_mean_sig <- subset(time_melt_mean, time_melt_mean$Peak_ID %in% AF_sig$Peak_ID)

time_melt_mean_sig%>%
  subset(Peak_ID == 421) %>%
  ggplot(aes(x=Time, y=Mean_Peak_Intensity,
             colour=Peak_ID))+
  geom_line()+
  theme_minimal()

###------

## Plotting the correlations of correlating metabolites to DAS44 changes, and with similar metabolites
cor_plots <- function(resp_diff_df, sig_table){
  resp_cor <- resp_diff_df[,-c(1:9)]
  resp_cor <- map_df(resp_cor, as.numeric)
  resp_cor_sig <- as.data.frame(t(resp_cor))
  resp_cor_sig$Peak_ID <- rownames(resp_cor_sig)
  resp_cor_sig$Peak_ID <- gsub('X', '', resp_cor_sig$Peak_ID)
  resp_cor_sig$Peak_ID <- as.numeric(resp_cor_sig$Peak_ID)
  sig_id <- sig_table[,c(1,30)]
  resp_cor_sig_sel <- inner_join(sig_id,resp_cor_sig, by='Peak_ID')
  resp_cor_sig_sel <- subset(resp_cor_sig_sel, resp_cor_sig_sel$Peak_ID != 1548)
  resp_cor_sig_sel <- as.data.frame(resp_cor_sig_sel)
  resp_cor_sig_sel <- distinct(resp_cor_sig_sel, Putative_Metabolite, .keep_all = TRUE )
  rownames(resp_cor_sig_sel) <- resp_cor_sig_sel$Putative_Metabolite
  resp_cor_sig_sel <- resp_cor_sig_sel[,-c(1,2)]
  resp_cor_sig_sel <- t(resp_cor_sig_sel)
  res <- cor(resp_cor_sig_sel)
  heatmap <- ggcorrplot(cor(resp_cor_sig_sel),
                        #outline.col = "white",
                        #method='circle',
                        p.mat = cor_pmat(resp_cor_sig_sel), 
                        insig='blank',
                        #type='lower',
                        hc.order=TRUE)
  correlations_melted <- melt(res) 
  names(correlations_melted) <- c('Metabolite1' ,'Metabolite2', 'Correlation')
  colours <- c("blue","white", "red")
  dd <- as.dist((1-cor(t(res))))
  ab <-plot(hclust(dd, method="complete"),main='', xlab="", sub="")
  hc <- hclust(dd, method="complete")
  colors = c("#0066CC", "#FF0033", "#66CC00", "purple", 'brown')
  clus4 = cutree(hc, 5)
  par(mar = c(1, 1, 1, 1))
  dendrogram <-plot(as.phylo(hc), 
                    type = "phylogram", 
                    tip.color = colors[clus4],
                    cex =1,
                    font =8)
  print(heatmap)
  print(dendrogram)
  sig_id_sim <- sig_id
  print(sig_id_sim)
}

AB_heats <- cor_plots(resp_diff_AB, AB_sig)
AF_heats <- cor_plots(resp_diff_AF, AF_sig)
BF_heats <- cor_plots(resp_diff_BF, BF_sig)

AB_nope <- c(320,333,1548,199,415,730,215,208,218,322,198,414,671,461,938,421,932,667,419, 1173)

AB_dist <- read.csv('20210308_AB_pos_diff_das_cor.csv')
AF_dist <- read.csv('20210308_AB_pos_diff_das_cor.csv')

AB_dist <- AB_dist %>%
  subset(Peak_ID %notin% AB_nope)%>%
  distinct(Putative_Metabolite, .keep_all=TRUE)

AF_dist <- AF_dist %>%
  subset(Peak_ID %notin% AB_nope)%>%
  distinct(Putative_Metabolite, .keep_all=TRUE)

BF_dist <- BF_dist %>%
  subset(Peak_ID %notin% AB_nope)%>%
  distinct(Putative_Metabolite, .keep_all=TRUE)
  

AB_diff_kegg <- AB_sig[,c(30,31)]
AB_diff_kegg$Putative_Metabolite[AB_diff_kegg$Putative_Metabolite == 'Pyroglutamic acid'] <- 'L-1-Pyrroline 3-hydroxy-5-carboxylate'
AB_diff_kegg_ID<- inner_join(kegg,AB_diff_kegg, by='Putative_Metabolite')

### Extending the correlational analyses to investigate further which metabolic pathways are perturbed during treatment
# Which metabolites correlated with a more lenient adj. p-value cut-off?
cors_ext <- function(melt_table, disease_measure){
  melt_table$Peak_ID <- as.numeric(melt_table$Peak_ID)
  melt_table <- inner_join(melt_table,peak_metadata, by='Peak_ID')
  ints_nested <- melt_table %>%
    group_by (Peak_ID) %>%
    nest()
  ints_unnested <- resp_diff_AB_melt %>%
    unnest(cols=c())
  identical(resp_diff_AB_melt, ints_unnested)
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
  best_augmented$adj_p <- p.adjust(best_augmented$p.value, method='BH')
  best_augmented_sig <- subset(best_augmented,best_augmented$adj_p < 0.5)
  best_adj_hmdb <- inner_join(best_augmented_sig, peak_metadata, by='Peak_ID')
  best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='NA')
  best_adj_hmdb <- subset(best_adj_hmdb, best_adj_hmdb$Putative_Metabolite !='')
  gg <- best_adj_hmdb %>%
    subset(Putative_Metabolite != 'Galactonic/Gluconic acid') %>%
    ggplot(aes(x = DAS44, y=Peak_Intensity)) +
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

# between baseline and 3 months?
AB_cors_ext <- cors_ext(resp_diff_AB_melt)
AB_sig_ext <- sig_dist_peaks(AB_cors_ext)
AB_not <- c(297, 25, 6, 189, 101, 692, 564,908, 964, 711, 113, 1278)
AB_sig_sim_dist <- AB_sig_ext[,c(1,26)] %>%
  subset(Peak_ID %notin% AB_not) %>%
  distinct(Putative_Metabolite, .keep_all = TRUE)

# between baseline and 18 months?
AF_cors_ext <- cors_ext(resp_diff_AF_melt)
AF_sig_ext <- sig_dist_peaks(AF_cors_ext)
AF_not <- c(153, 266, 41, 692, 189)
AF_sig_sim_dist <- AF_sig_ext[,c(1,26)] %>%
  subset(Peak_ID %notin% AF_not) %>%
  distinct(Putative_Metabolite, .keep_all = TRUE)

# between 3 months and 18 months?
BF_cors_ext <- cors_ext(resp_diff_BF_melt)
BF_sig_ext <- sig_dist_peaks(BF_cors_ext)
BF_not <- c(297, 25, 6, 189, 1361, 89, 173, 153, 266, 164, 113, 101)
BF_sig_sim_dist <- BF_sig_ext[,c(1,26)] %>%
  subset(Peak_ID %notin% BF_not) %>%
  distinct(Putative_Metabolite, .keep_all = TRUE)

### Which metabolites behave similarly over these periods?
## AB
AB_sig_ext_dist <- AB_sig_ext %>%
  subset(Peak_ID %notin% AB_not) %>%
  distinct(Putative_Metabolite, .keep_all = TRUE)

AB_heats_ext <- cor_plots(resp_diff_AB, AB_sig_ext_dist)

# Kegg ID
AB_ext_kegg <- inner_join(kegg, AB_sig_ext_dist, by='Putative_Metabolite')
AB_ext_kegg_2 <- AB_ext_kegg[,c(2,3)]

## AF
AF_sig_ext_dist <- AF_sig_ext %>%
  subset(Peak_ID %notin% AF_not) %>%
  distinct(Putative_Metabolite, .keep_all = TRUE)
AF_heats_ext <- cor_plots(resp_diff_AF, AF_sig_ext_dist)

# Kegg ID
AF_ext_kegg <- inner_join(kegg, AF_sig_ext_dist, by='Putative_Metabolite')
AF_ext_kegg_2 <- AF_ext_kegg[,c(2,3)]

## BF
BF_sig_ext_dist <- BF_sig_ext %>%
  subset(Peak_ID %notin% BF_not) %>%
  distinct(Putative_Metabolite, .keep_all = TRUE)
BF_heats_ext <- cor_plots(resp_diff_BF, BF_sig_ext_dist)

# Kegg ID
BF_ext_kegg <- inner_join(kegg, BF_sig_ext_dist, by='Putative_Metabolite')
BF_ext_kegg_2 <- BF_ext_kegg[,c(2,3)]

# Changing metabolites across periods-------
ABF_sum <- read.csv('20210308_ABF_pos_diff_das_summ.csv')
ABF_sum <- subset(ABF_sum,ABF_sum$Peak_ID != 'NA')
ABF_sum_code <- ABF_sum
ABF_sum_code$Correlation <- ifelse(ABF_sum_code$Correlation == 'Positive',1,-1)
ABF_sum_code<- ABF_sum_code[-1]

ABF_cast <- dcast(ABF_sum_code, Putative_Metabolite ~Period)
ABF_cast <- ABF_cast[,c(1,2,4,3)]
ABF_cast[is.na(ABF_cast)] <- 0

ABF_melt <- melt(ABF_cast)
names(ABF_melt) <- c('Putative_Metabolite', 'Period', 'Correlation')
ABF_melt$Correlation <- as.factor(ABF_melt$Correlation)
ABF_melt$Cor_2 <- 0
ABF_melt$Cor_2[ABF_melt$Correlation == '1'] <- 1
ABF_melt$Cor_2[ABF_melt$Correlation == ''] <- 0

ABF_melt %>%
  subset(Putative_Metabolite != 'linolenic acid') %>%
  ggplot(aes(x=Period, y= fct_reorder(Putative_Metabolite, Cor_2), fill=Correlation))+
  scale_fill_manual(values = c('#3399FF', '#FFFFFF', '#FF3333'),
                    name = "Correlation of \nΔmetabolite level \nwith ΔDAS44", labels = c("Negative", "None", "Positive"))+
  geom_tile(colour='black')+
  theme_minimal()+
  labs(x='Period of Treatment',
       y='Putative Metabolite')+
  scale_x_discrete(breaks=c("AB", "BF", "AF"),
                   labels=c("Baseline to 3 Months", "3 Months to 18 Months", "Baseline to 18 months"))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

ABF_kegg <- inner_join(kegg, ABF_melt, by='Putative_Metabolite')%>%
  distinct(Putative_Metabolite, .keep_all=TRUE)
ABF_kegg_sim <- ABF_kegg[,c(2,3)]


## Write csv files for combining positive and negative ion modes
AB_dist_2 <-AB_dist
AB_dist_2$Ion_Mode <- 'Positive'
write.csv(AB_dist_2, '20210308_AB_pos_diff_das_cor.csv')

AF_dist_2 <-AF_dist
AF_dist_2$Ion_Mode <- 'Positive'
write.csv(AF_dist_2, '20210308_AF_pos_diff_das_cor.csv')

BF_dist_2 <-BF_dist
BF_dist_2$Ion_Mode <- 'Positive'
write.csv(BF_dist_2, '20210308_BF_pos_diff_das_cor.csv')

# Sample Peak intensities and differential intensities
# AB
write.csv(resp_int_AB_2, '20210309_AB_pos_resp_ints.csv')
write.csv(resp_dbl_diff_AB, '20210309_AB_pos_resp_diff.csv')

# BF
write.csv(resp_int_BF_2, '20210309_BF_pos_resp_ints.csv')
write.csv(resp_dbl_diff_BF, '20210309_BF_pos_resp_diff.csv')

# AF
write.csv(resp_int_AF_2, '20210309_AF_pos_resp_ints.csv')
write.csv(resp_dbl_diff_AF, '20210309_AF_pos_resp_diff.csv')

