### Combining positive and negative ion mode results from 20210308_ABF_pos_initial.R and 20210303_ABG_neg_initial.R
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
#################################################################################################
### files------
peak_metadata_neg <- read.csv(file='20210217_neg_peak_metadata.csv', header=TRUE)
peak_metadata_pos <- read.csv(file='20200521_Taser_POS_Peakdata.csv', header=TRUE, row.names=1)
kegg <- read.csv('20210305_kegg_list.csv')
kegg <- distinct(kegg, Kegg_code, .keep_all = TRUE)
## Helpful functions
`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)

# Combine positive and negative peak metadata
peak_metadata_neg$Peak_ID_Ion <- paste0('Neg_X',peak_metadata_neg$Peak_ID)
peak_metadata_pos$Peak_ID_Ion <- paste0('Pos_X',peak_metadata_pos$Peak_ID)

peak_metadata_neg_2 <- peak_metadata_neg[,c(3,6:10,1,5,11)]
peak_metadata_pos_2 <- peak_metadata_pos[,c(1:6,8,9,10)]
names(peak_metadata_pos_2)

peak_metadata_comb <- rbind.data.frame(peak_metadata_neg_2, peak_metadata_pos_2)

#################################################################################################
# Setting up a list of metabolites to include in the correlation analysis- from 20210303_ABG_neg/pos_initial.R
# function
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

# AB correlating peaks
AB_pos_cor <- read.csv('20210308_AB_pos_diff_das_cor.csv')
AB_neg_cor <- read.csv('20210308_AB_neg_diff_das_cor.csv')

# AF correlating peaks
AF_pos_cor <- read.csv('20210308_AF_pos_diff_das_cor.csv')
AF_neg_cor <- read.csv('20210308_AF_neg_diff_das_cor.csv')

# BF correlating peaks
BF_pos_cor <- read.csv('20210308_BF_pos_diff_das_cor.csv')
BF_neg_cor <- read.csv('20210308_BF_neg_diff_das_cor.csv')

AB_pos_cor <- AB_pos_cor%>%
  subset(Putative_Metabolite != '')%>%
  subset(Putative_Metabolite != 'NA')

AF_pos_cor <- AF_pos_cor%>%
  subset(Putative_Metabolite != '')%>%
  subset(Putative_Metabolite != 'NA')

BF_pos_cor <- BF_pos_cor%>%
  subset(Putative_Metabolite != '')%>%
  subset(Putative_Metabolite != 'NA')

AB_neg_cor <- AB_neg_cor%>%
  subset(Putative_Metabolite != '')%>%
  subset(Putative_Metabolite != 'NA')

AF_neg_cor <- AF_neg_cor%>%
  subset(Putative_Metabolite != '')%>%
  subset(Putative_Metabolite != 'NA')

BF_neg_cor <- BF_neg_cor%>%
  subset(Putative_Metabolite != '')%>%
  subset(Putative_Metabolite != 'NA')


comb_AB <- rbind.data.frame(AB_pos_cor, AB_neg_cor)
comb_AF <- rbind.data.frame(AF_pos_cor, AF_neg_cor)
comb_BF <- rbind.data.frame(BF_pos_cor, BF_neg_cor)
comb_AB <- comb_AB[,-1] ; comb_AF <- comb_AF[,-1] ; comb_BF <- comb_BF[,-1]
comb_AB$Peak_ID_Ion <- paste0(comb_AB$Ion_Mode,'_X', comb_AB$Peak_ID)
comb_AF$Peak_ID_Ion <- paste0(comb_AF$Ion_Mode,'_X', comb_AF$Peak_ID)
comb_BF$Peak_ID_Ion <- paste0(comb_BF$Ion_Mode,'_X', comb_BF$Peak_ID)

# Combine the positive and negative resp_diff dataframes
resp_diff_pos <- read.csv('20210204_AF_resp_diff_pos.csv', header=TRUE)
resp_diff_neg <- read.csv('20201105_resp_diff.csv', header=TRUE)

resp_pos <- resp_diff_pos[,c(13:1596)]
resp_neg <- resp_diff_neg[,c(14:1471)]

rownames(resp_pos) <- paste0(resp_diff_pos$Sample_Name,'_P')
rownames(resp_neg) <- paste0(resp_diff_neg$Sample_Name, '_N')

names(resp_pos) <- paste0('Positive_', names(resp_pos))
names(resp_neg) <- paste0('Negative_', names(resp_neg))

comb_mets <- cbind.data.frame(resp_pos, resp_neg)
comb_cut <- as.data.frame(t(comb_mets))
colnames(comb_cut) <- gsub("(.*)_.*","\\1",colnames(comb_cut)) # remove everyting after '_'
comb_cut$Peak_ID_Ion <- rownames(comb_cut)
comb_cut <- comb_cut[,c(ncol(comb_cut),1:(ncol(comb_cut)-1))]

#################################################################################################
# Generate response/peak intensities matrices based on combined ion mode data for AB, AF and BF periods
# AB
comb_ID_AB <- inner_join(comb_AB, comb_cut, by='Peak_ID_Ion')
rownames(comb_ID_AB) <- paste0(comb_ID_AB$Putative_Metabolite, '_',comb_ID_AB$Ion_Mode)
comb_ID_AB_mat <- comb_ID_AB[,-c(1:4)]
comb_ID_AB_mat <- t(comb_ID_AB_mat)
colnames(comb_ID_AB_mat) <- gsub("(.*)_.*","\\1",colnames(comb_ID_AB_mat)) # remove everyting after '_'

# AF
comb_ID_AF <- inner_join(comb_AF, comb_cut, by='Peak_ID_Ion')
rownames(comb_ID_AF) <- paste0(comb_ID_AF$Putative_Metabolite, '_', comb_ID_AF$Ion_Mode)
comb_ID_AF_mat <- comb_ID_AF[,-c(1:4)]
comb_ID_AF_mat <- t(comb_ID_AF_mat)
colnames(comb_ID_AF_mat) <- gsub("(.*)_.*","\\1",colnames(comb_ID_AF_mat))

#BF
comb_ID_BF <- inner_join(comb_BF, comb_cut[,-64], by='Peak_ID_Ion')
comb_ID_BF$Putative_Metabolite[comb_ID_BF$Putative_Metabolite == 'GPC'] <- ''
comb_ID_BF$Putative_Metabolite[comb_ID_BF$Peak_ID == 1108] <- ''
comb_ID_BF$Putative_Metabolite[comb_ID_BF$Peak_ID == 113] <- ''
comb_ID_BF$Putative_Metabolite[comb_ID_BF$Peak_ID == 1226] <- ''
comb_ID_BF$Putative_Metabolite[comb_ID_BF$Peak_ID == 314] <- ''
comb_ID_BF$Putative_Metabolite[comb_ID_BF$Peak_ID == 1109] <- ''
comb_ID_BF <- subset(comb_ID_BF, comb_ID_BF$Putative_Metabolite != '')
rownames(comb_ID_BF) <- paste0(comb_ID_BF$Putative_Metabolite,'_', comb_ID_BF$Ion_Mode)
comb_ID_BF_mat <- comb_ID_BF[,-c(1:4)]
comb_ID_BF_mat <- t(comb_ID_BF_mat)
colnames(comb_ID_BF_mat) <- gsub("(.*)_.*","\\1",colnames(comb_ID_BF_mat))

#################################################################################################
## Correlation analysis ----

AB_cor_heat <- ggcorrplot(cor(comb_ID_AB_mat),
                          #outline.col = "white",
                          #method='circle',
                          p.mat = cor_pmat(comb_ID_AB_mat), 
                          insig='blank',
                          #type='lower',
                          hc.order=TRUE)

comb_AB_kegg <- as.data.frame(comb_ID_AB[,c(1,2)])%>%
  inner_join(kegg, by='Putative_Metabolite')


AF_cor_heat <- ggcorrplot(cor(comb_ID_AF_mat),
                          #outline.col = "white",
                          #method='circle',
                          p.mat = cor_pmat(comb_ID_AF_mat), 
                          insig='blank',
                          #type='lower',
                          hc.order=TRUE)

comb_AF_kegg <- as.data.frame(comb_ID_AF[,c(1,2)]) %>%
  inner_join(kegg, by='Putative_Metabolite')

BF_cor_heat <- ggcorrplot(cor(comb_ID_BF_mat),
                          #outline.col = "white",
                          #method='circle',
                          p.mat = cor_pmat(comb_ID_BF_mat), 
                          insig='blank',
                          #type='lower',
                          hc.order=TRUE)

comb_BF_kegg <- as.data.frame(comb_ID_BF[,c(1,2)]) %>%
  inner_join(kegg, by='Putative_Metabolite')

#################################################################################################
#Kegg 
# AB
AB_kegg <- comb_AB[,c(1,2)]
AB_kegg <- inner_join(kegg, AB_kegg, by='Putative_Metabolite')

AF_kegg <- comb_AF[,c(1,2)]
AF_kegg <- inner_join(kegg, AF_kegg, by='Putative_Metabolite')

BF_kegg <- comb_BF[,c(1,2)]
BF_kegg <- inner_join(kegg, BF_kegg, by='Putative_Metabolite')

#################################################################################################
### Investigating whether predictions on DAS44 outcomes can be made through the combined positive and negative ion mode data
# files -----
#AB
resp_ints_neg_AB <- read.csv('20210309_AB_neg_resp_ints.csv', header=TRUE)
resp_diff_neg_AB <- read.csv('20210309_AB_neg_resp_diff.csv', header=TRUE)

resp_ints_pos_AB <- read.csv('20210309_AB_pos_resp_ints.csv', header=TRUE)
resp_diff_pos_AB <- read.csv('20210309_AB_pos_resp_diff.csv', header=TRUE)

#BF
resp_ints_neg_BF <- read.csv('20210309_BF_neg_resp_ints.csv', header=TRUE)
resp_diff_neg_BF <- read.csv('20210309_BF_neg_resp_diff.csv', header=TRUE)

resp_ints_pos_BF <- read.csv('20210309_BF_pos_resp_ints.csv', header=TRUE)
resp_diff_pos_BF <- read.csv('20210309_BF_pos_resp_diff.csv', header=TRUE)

# AF
resp_ints_neg_AF <- read.csv('20210309_AF_neg_resp_ints.csv', header=TRUE)
resp_diff_neg_AF <- read.csv('20210309_AF_neg_resp_diff.csv', header=TRUE)

resp_ints_pos_AF <- read.csv('20210309_AF_pos_resp_ints.csv', header=TRUE)
resp_diff_pos_AF <- read.csv('20210309_AF_pos_resp_diff.csv', header=TRUE)

## Combining positive and negative ion modes
# rename peaks for negative and positive ion modes
names(resp_ints_neg_AB) <- ifelse(names(resp_ints_neg_AB) %like% 'X', paste0('Neg_', names(resp_ints_neg_AB)), names(resp_ints_neg_AB))
names(resp_ints_neg_BF) <- ifelse(names(resp_ints_neg_BF) %like% 'X', paste0('Neg_', names(resp_ints_neg_BF)), names(resp_ints_neg_BF))
names(resp_ints_neg_AF) <- ifelse(names(resp_ints_neg_AF) %like% 'X', paste0('Neg_', names(resp_ints_neg_AF)), names(resp_ints_neg_AF))
names(resp_diff_neg_AB) <- ifelse(names(resp_diff_neg_AB) %like% 'X', paste0('Neg_', names(resp_diff_neg_AB)), names(resp_diff_neg_AB))
names(resp_diff_neg_BF) <- ifelse(names(resp_diff_neg_BF) %like% 'X', paste0('Neg_', names(resp_diff_neg_BF)), names(resp_diff_neg_BF))
names(resp_diff_neg_AF) <- ifelse(names(resp_diff_neg_AF) %like% 'X', paste0('Neg_', names(resp_diff_neg_AF)), names(resp_diff_neg_AF))

names(resp_ints_pos_AB) <- ifelse(names(resp_ints_pos_AB) %like% 'X', paste0('Pos_', names(resp_ints_pos_AB)), names(resp_ints_pos_AB))
names(resp_ints_pos_BF) <- ifelse(names(resp_ints_pos_BF) %like% 'X', paste0('Pos_', names(resp_ints_pos_BF)), names(resp_ints_pos_BF))
names(resp_ints_pos_AF) <- ifelse(names(resp_ints_pos_AF) %like% 'X', paste0('Pos_', names(resp_ints_pos_AF)), names(resp_ints_pos_AF))
names(resp_diff_pos_AB) <- ifelse(names(resp_diff_pos_AB) %like% 'X', paste0('Pos_', names(resp_diff_pos_AB)), names(resp_diff_pos_AB))
names(resp_diff_pos_BF) <- ifelse(names(resp_diff_pos_BF) %like% 'X', paste0('Pos_', names(resp_diff_pos_BF)), names(resp_diff_pos_BF))
names(resp_diff_pos_AF) <- ifelse(names(resp_diff_pos_AF) %like% 'X', paste0('Pos_', names(resp_diff_pos_AF)), names(resp_diff_pos_AF))

# combine
AB_ints <- cbind.data.frame(resp_ints_neg_AB, resp_ints_pos_AB[,-c(1:12)])
BF_ints <- cbind.data.frame(resp_ints_neg_BF, resp_ints_pos_BF[,-c(1:12)])
AF_ints <- cbind.data.frame(resp_ints_neg_AF, resp_ints_pos_AF[,-c(1:12)])

AB_diff <- cbind.data.frame(resp_diff_neg_AB, resp_diff_pos_AB[,-c(1:10)])
BF_diff <- cbind.data.frame(resp_diff_neg_BF, resp_diff_pos_BF[,-c(1:10)])
AF_diff <- cbind.data.frame(resp_diff_neg_AF, resp_diff_pos_AF[,-c(1:10)])

## Differential metabolite abundance
peak_list<- AF_ints
peak_list <- peak_list[,-c(1:13)]
peak_list <- as.data.frame(t(peak_list))
peak_list$ID <- rownames(peak_list)
peak_list$Peak_ID <- peak_list$ID
peak_list$Peak_ID <- gsub('Pos_X','', peak_list$Peak_ID)
peak_list$Peak_ID <- gsub('Neg_X','', peak_list$Peak_ID)
peak_list$Peak_ID <- as.numeric(peak_list$ID)

# limma function
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
  toptable$Peak_ID_Ion <- toptable$Peak_ID
  toptable$Peak_ID <- gsub('Neg_X', '', toptable$Peak_ID)
  toptable$Peak_ID <- gsub('Pos_X', '', toptable$Peak_ID)
  toptable$Peak_ID <- as.numeric(as.character(toptable$Peak_ID))
  toptable <- toptable[,c(ncol(toptable),1:(ncol(toptable)-1))]
  toptable <- inner_join(toptable, peak_metadata_comb, by = 'Peak_ID_Ion')
}

qvals <- function(limma_table){
  pi0 <- 2*mean(limma_table$P.Value > 0.05)
  lfdrvals <- lfdr(limma_table$P.Value, pi0)
  qobj <- qvalue(limma_table$P.Value)
  hist(qobj)
}

limma_ID <- function(limma_table){
  limma_table_ID <- limma_table
  limma_table_ID$Sig <- 0
  limma_table_ID$Sig <- ifelse(limma_table_ID$adj.P.Val <0.05, '< 0.05', '> 0.05') 
  limma_table_ID$Sig_Peaks <- 0
  limma_table_ID$Sig_Peaks <- ifelse(limma_table_ID$Sig =='< 0.05' & limma_table_ID$identification != '',limma_table_ID$Peak_ID.x, '')
  limma_table_ID$Sig_Names <-0
  limma_table_ID$Sig_Names <- ifelse(limma_table_ID$Sig =='< 0.05' & limma_table_ID$identification != '',limma_table_ID$Putative_Metabolite, '')
  limma_table_ID%>%
    subset(Putative_Metabolite != 'Sulfasalazine')%>%
    subset(Putative_Metabolite != 'Galactonic/Gluconic acid')%>%
    subset(Putative_Metabolite != 'Valerylglycine')%>%
    ggplot(aes(x=logFC, y=-log10(P.Value), 
               colour=Sig, 
               group=Sig)) +
    geom_point (alpha=0.7) +
    theme_minimal() +
    labs (x='LogFC',
          y='-Log p-value',
          colour='Adjusted \np-value')+
    geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names),
                    box.padding =1,
                    size=4,
                    max.overlaps = Inf,
                    position = position_jitter(seed = 1),
                    arrow = arrow(length = unit(0.0015, "npc"))) +  
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))+
    scale_color_brewer(palette = "Set1",direction=-1)+
    ylim(0,10)+
    xlim(-2.5, 2.5)
}

AF_ints_2 <- AF_ints[,-c(1,2,4:13)]
rownames(AF_ints_2) <- paste0(rownames(AF_ints_2), AF_ints_2$time)
AF_ints_2<- AF_ints_2[,-1]
AF_ints_2 <- AF_ints_2[,!names(AF_ints_2) == "Batch"]
AF_ints_t <- t(AF_ints_2)
colnames(AF_ints_t)[colnames(AF_ints_t) %like% 'A'] <- 'A'
colnames(AF_ints_t)[colnames(AF_ints_t) %like% 'F'] <- 'F'

AF_limma <- limma_fun(AF_ints_t, 3300, 'A', 'F')
AF_q <- qvals(AF_limma)
AF_ID <- limma_ID(AF_limma)
AFtab <- AF_ID$data

# Baseline metabolites to predict 18 month outcomes -----
# Use peak intensities from AB_ints with AF_diff DAS44 outcomes
# Categorise responses
AF_diff$Response <- AF_diff$DAS44
AF_diff <- AF_diff[,c(2:10,3053,11:3052)]
AB_ints_2 <- AB_ints[,c(2,6,14:3056)]

AB_ints_2 <- subset(AB_ints_2, AB_ints_2$Sample %in% AF_diff$Sample)
AF_diff_2 <- subset(AF_diff, AF_diff$Sample %in% AB_ints_2$Sample)

ABF_pred <- cbind.data.frame(AF_diff_2[1:10], AB_ints_2[2:3045])
names(ABF_pred)[c(3,11)]<- c('ΔDAS44', 'End_DAS44')
ABF_pred$Response <- with(ABF_pred, ifelse(End_DAS44 <= 1.6, 'Positive', 'Negative'))
ABF_pred$Sample <- paste0(ABF_pred$Sample, c('A', 'B'))
ABF_pred_2 <- subset(ABF_pred, ABF_pred$Sample %like% 'A')
ABF_pred_3 <- ABF_pred_2
ABF_pred_2 <- ABF_pred_2[,-c(1:9,11)]
ABF_pred_2$Response[ABF_pred_2$Response == 'Positive'] <- 1
ABF_pred_2$Response[ABF_pred_2$Response == 'Negative'] <- 0

ABF_pred_2$Response <- as.factor(ABF_pred_2$Response)

# Same as above but use end DAS44 as the marker of response
ABF_pred_end <- ABF_pred_3
ABF_pred_end$Sample <- gsub('A', '', ABF_pred_end$Sample)
AF_ints_end <- AF_ints
AF_ints_end <- subset(AF_ints_end,AF_ints_end$time == 'F')

AF_ints_end <- subset(AF_ints_end,AF_ints_end$Sample %in% ABF_pred_end$Sample)

ABF_pred_endDAS <- cbind.data.frame(AF_ints_end$DAS44, ABF_pred_end[12:3054])
names(ABF_pred_endDAS)[1] <- 'Response'
ABF_pred_endDAS$Response <- ifelse(ABF_pred_endDAS$Response < 1.6, 1, 0)
ABF_pred_endDAS$Response <- as.factor(ABF_pred_endDAS$Response )

# Logistic regression model -------
comb_AB_peaks <- comb_AB
comb_AB_peaks$Peak_ID_Ion <-  gsub('itive','', comb_AB_peaks$Peak_ID_Ion);comb_AB_peaks$Peak_ID_Ion <-  gsub('ative','', comb_AB_peaks$Peak_ID_Ion)
comb_AB_peaks$Peak_ID_Ion

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
  print(plot(prf), abline(0,1))
}
set.seed(42)
index <- createDataPartition(ABF_pred_endDAS$Response, p = 0.65, list = FALSE)
train_data <- ABF_pred_endDAS[index, ]
test_data  <- ABF_pred_endDAS[-index, ]

#model <- glm(Response ~ Neg_X787 + Pos_X447   + Pos_X304  + Neg_X232   + Neg_X48 + Neg_X642 + Pos_X569,             family = binomial(link = "logit"),  data = train_data)

model <- glm(Response ~ Neg_X787 +  Neg_X48 + Pos_X505+ Pos_X1126+ Pos_X304 +Pos_X569+ Neg_X642+
               Neg_X642 + Pos_X396 + Neg_X421 + Pos_X447 + Neg_X642 + Neg_X1062 + Neg_X770
               ,             
             family = binomial(link = "logit"),  data = train_data)

model <- glm(Response ~ Neg_X787  + Pos_X1126+ Pos_X304 +Pos_X569+
               Neg_X642 + Pos_X396 + Neg_X421 + Pos_X447 + Neg_X1062 
             
             
             ,             
             family = binomial(link = "logit"),  data = train_data)


model <- glm(Response ~ Neg_X787  + Pos_X1126+ Pos_X304+
               Neg_X642 + Neg_X421  + Neg_X1062,
             family = binomial(link = "logit"),  data = train_data)

summary(model)

lrm_output(model)

# Identifying the peaks
peaks_lrm <- c('Neg_X787', 'Pos_X1126', 'Pos_X304', 
               'Neg_X421', 'Neg_X1062', 'Neg_X642')

peaks_lrm_id <- subset(comb_AB_peaks, comb_AB_peaks$Peak_ID_Ion %in% peaks_lrm)

### Baseline correlations with 18 month DAS44 outcomes --------
ABF_cors <- subset(ABF_pred,ABF_pred$Sample %like% 'A') # melt to DAS44, peak intesnities, samples and peaks

grep('Batch', colnames(ABF_cors))

ABF_cors_sim <- ABF_cors[,-c(2:11, 1470)]
ABF_cors_sim$Sample <- gsub('A', '', ABF_cors_sim$Sample)
ABF_melt <- melt(ABF_cors_sim)
names(ABF_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
ABF_melt$DAS44 <- ABF_cors$ΔDAS44


ABF_cor_res <- cors(ABF_melt) # no significant correlations found

# 3 month changes correlating with 18 month DAS44 outcomes ----
ABF_pred$Sample <- gsub('A', '', ABF_pred$Sample);ABF_pred$Sample <- gsub('B', '', ABF_pred$Sample)
ABF_cors_2 <- ABF_pred[,-c(2:11, 1470)]
AB_diff_das <- aggregate(.~Sample, ABF_cors_2, diff, na.rm=TRUE)

AB_AF_pred <- cbind.data.frame(AF_diff_sing[1:10], AB_diff_das[2:3043])

AB_AF_pred_lm <- AB_AF_pred[,-c(2:10)]

AB_AF_melt <- melt(AB_AF_pred_lm)
names(AB_AF_melt) <- c('Sample_Name', 'Peak_ID', 'Peak_Intensity')
AB_AF_melt$DAS44 <- AB_AF_pred$DAS44

AB_Af_cor_res <- cors(AB_AF_melt)# no significant correlations found between baseline and 3 month metabolic change and DAS44 response over 18 months



## Combine correlating metabolites across periods to map changes
comb_AB_2 <- comb_AB
comb_BF_2 <- comb_BF
comb_AF_2 <- comb_AF

comb_AB_2$Period <- 'AB'
comb_BF_2$Period <- 'BF'
comb_AF_2$Period <- 'AF'

comb_all <- rbind.data.frame(comb_AB_2, comb_BF_2, comb_AF_2)
comb_all$Correlation <- 0

#write.csv(comb_all, '20210309_ABF_comb_cor_sum.csv')
comb_all_cor <- read.csv('20210309_ABF_comb_cor_sum.csv', header=TRUE, row.names = 1)

comb_all_cor <- subset(comb_all_cor,comb_all_cor$Peak_ID != 'NA')
comb_all_code <- comb_all_cor
comb_all_code$Correlation[comb_all_code$Correlation == 'Positive'] <-'1'
comb_all_code$Correlation[comb_all_code$Correlation == 'Negative'] <- '-1'
comb_all_code$Correlation <- as.numeric(comb_all_code$Correlation)

ABF_cast <- comb_all_code[,-c(1,3,4)]
ABF_cast <- dcast(ABF_cast, Putative_Metabolite ~Period, fun.aggregate = mean, na.rm = TRUE)
ABF_cast <- ABF_cast[,c(1,2,4,3)]
ABF_cast[is.na(ABF_cast)] <- 0

ABF_melt <- melt(ABF_cast)
names(ABF_melt) <- c('Putative_Metabolite', 'Period', 'Correlation')
ABF_melt$Correlation <- as.factor(ABF_melt$Correlation)

ABF_melt$Period <- factor(ABF_melt$Period, levels=c("AB", "BF", "AF"))

ABF_melt %>%
  subset(Putative_Metabolite !='4-Hydroxybenzoic acid') %>%
  subset(Putative_Metabolite !='GPC') %>%
  ggplot(aes(x=Period, #y= reorder(Putative_Metabolite, Correlation), 
             y= Putative_Metabolite,
             fill=Correlation))+
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

x5p <- as.data.frame(c('1','Xylulose 5-phosphate', 'C00231', 'AB', '-1'))
x5p<- as.data.frame(t(x5p))
names(x5p) <- names(ABF_kegg)
rownames(x5p) <- '1'

ABF_kegg<- rbind.data.frame(ABF_kegg, x5p)
names(ABF_kegg)


