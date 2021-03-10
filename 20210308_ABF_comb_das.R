pos_limma_AF <- read.csv('20210218_AF_pos_diff.csv')
neg_limma_AF <- read.csv('20210218_AF_neg_diff.csv')

names(pos_limma_AF)
names(neg_limma_AF)

pos_limma_AF_cut <- pos_limma_AF[,c(2,3,6,7,9,10,16)]
neg_limma_AF_cut <- neg_limma_AF[,c(2,3,6,7,13,14,12)]

names(pos_limma_AF_cut)
names(neg_limma_AF_cut)

pos_limma_AF_cut$Ion_Mode <- 'Positive'
neg_limma_AF_cut$Ion_Mode <- 'Negative'

comb_limma <- rbind.data.frame(pos_limma_AF_cut, neg_limma_AF_cut)
comb_limma$Sig <- 0
comb_limma$Sig <- ifelse(comb_limma$adj.P.Val <0.05, '< 0.05', '> 0.05') 
comb_limma$Sig_Peaks <- ifelse(comb_limma$adj.P.Val<0.05, comb_limma$Putative_Metabolite, '')
comb_limma_dist <- distinct(comb_limma, Peak_ID, .keep_all = TRUE)
comb_limma_dist <- subset(comb_limma_dist,comb_limma_dist$adj.P.Val < 0.05)
comb_limma_dist <- comb_limma_dist[,-c(9,10)]

#write.csv(comb_limma_dist, '20210223_AF_comb_diff.csv')
das_mets <- c('L-Tyrosine', 'L-Serine', '5-Hydroxytryptophol', 'Linoleyl carnitine', '2-Phenylacetamide', 'N-Acetylleucine'              )

comb_limma$Sig_Peaks_DAS <- ifelse(comb_limma$Putative_Metabolite %in% das_mets, comb_limma$Putative_Metabolite , '')
comb_limma$Sig_Peaks <- 0
comb_limma$Sig_Peaks <- ifelse(comb_limma$adj.P.Val <0.05, comb_limma$Putative_Metabolite, '')

rem_mets <- c(45,75,550,604, 944,1145, 1278, 1331)
`%notin%` <- Negate(`%in%`)

comb_limma%>%
  subset(Peak_ID %notin% rem_mets) %>%
  ggplot(aes(x=logFC, y=-log10(P.Value), 
             colour=Sig, 
             group=Sig)) +
  geom_point (alpha=0.7) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Adjusted \np-value')+
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks),
                  box.padding =1.2,
                  size=3,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.005, "npc"))) +  
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  scale_color_brewer(palette = "Set1",direction=-1)+
  xlim(-2,2)+
  ylim(0,7.5)

# Logistic regression of metabolites for predicting DAS44
pos_resp <- read.csv('20210218_AF_pos_LGR_matrix', row.names=1)
neg_resp <- read.csv('20210218_AF_neg_LGR_matrix', row.names=1)

names(pos_resp)[1:1584] <- paste0('Pos_', names(pos_resp)[1:1584])
names(neg_resp)[2:1459] <- paste0('Neg_', names(neg_resp)[2:1459])

summary(pos_resp$Response == neg_resp$Response)

comb_resp <- cbind.data.frame(neg_resp, pos_resp[,-1585])

set.seed(42)
index <- createDataPartition(comb_resp$Response, p = 0.7, list = FALSE)
train_data <- comb_resp[index, ]
test_data  <- comb_resp[-index, ]

#model <- glm(Response ~ Neg_X243 + Neg_X189 +  Neg_X182 + Neg_X351 + Neg_X1028 +
#            Pos_X423  + Pos_X25 + Pos_X1038 + Pos_X1037 + Pos_X1030 + Pos_X274 + Pos_X484,
#          family=binomial(link='logit'),data=train_data)

model_pos <- glm(Response ~ 
                   Pos_X423  + Pos_X25 + Pos_X1038 + Pos_X1037   + Pos_X484,
                 family=binomial(link='logit'),data=train_data)

model_neg <- glm(Response ~ Neg_X243 + Neg_X189 +   Neg_X351 + Neg_X1028,
             family=binomial(link='logit'),data=train_data)

#model <- glm(Response ~   Neg_X351 + Neg_X1028+ Neg_X589 + Neg_X191 +               Pos_X423  + Pos_X25 + Pos_X1038    + Pos_X484 + Pos_X400,
 #            family=binomial(link='logit'),data=train_data)

model_good <- glm(Response ~   Neg_X1028+  
               Pos_X423  +  Pos_X1038    +  Pos_X400+  Pos_X484,
             family=binomial(link='logit'),data=train_data)


model <- glm(Response ~   Neg_X1028+Pos_X423  +  Pos_X1038    +  Pos_X400+  
                Pos_X107 + Pos_X374,
             family=binomial(link='logit'),data=train_data)

summary(model)
anova(model, test="Chisq")
pR2(model)

fitted.results <- predict(model,newdata=test_data,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test_data$Response)
print(paste('Accuracy',1-misClasificError))

library(ROCR)
p <- predict(model,newdata=test_data,type='response')
pr <- prediction(p, test_data$Response)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

# Plot these metabolites in the volcano plot
pos_limma_AF_cut_2 <- pos_limma_AF_cut
neg_limma_AF_cut_2 <- neg_limma_AF_cut

pos_limma_AF_cut_2$Peak_ID <- paste0('P_',pos_limma_AF_cut_2$Peak_ID)
neg_limma_AF_cut_2$Peak_ID <- paste0('N_',neg_limma_AF_cut_2$Peak_ID)

comb_limma_ann <- rbind.data.frame(pos_limma_AF_cut_2, neg_limma_AF_cut_2)
comb_limma_ann$Sig <- 0
comb_limma_ann$Sig <- ifelse(comb_limma_ann$adj.P.Val <0.05, '< 0.05', '> 0.05') 
comb_limma_ann$Sig_Peaks <- ifelse(comb_limma_ann$adj.P.Val<0.05, comb_limma_ann$Putative_Metabolite, '')
comb_limma_dist_2 <- distinct(comb_limma_ann, Putative_Metabolite, .keep_all = TRUE)

das_mets <- c('P_423', 'P_25' ,'P_1038', 'P_484', 'P_400',
              'N_351', 'N_1028', 'N_589')

comb_limma_das <- comb_limma_ann
comb_limma_das$Sig_Peaks_DAS <- ifelse(comb_limma_ann$Peak_ID %in% das_mets, comb_limma$Putative_Metabolite , '')

comb_limma_das%>%
  ggplot(aes(x=logFC, y=-log10(P.Value), 
             colour=Sig, 
             group=Sig)) +
  geom_point (alpha=0.7) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Adjusted \np-value')+
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Peaks_DAS),
                  box.padding =1.2,
                  size=3,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.005, "npc"))) +  
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  scale_color_brewer(palette = "Set1",direction=-1)+
  xlim(-1.5,1)+
  ylim(0,4)


### Differential analysis
## LRM
# Logistic regression of metabolites for predicting DAS44
pos_resp <- read.csv('20210223_AF_pos_das_resp.csv', row.names=1)
neg_resp <- read.csv('20210223_AF_neg_das_resp.csv', row.names=1)

summary(pos_resp$Response == neg_resp$Response)

comb_resp <- cbind.data.frame(neg_resp, pos_resp[,-1])
comb_resp$Response <- as.factor(comb_resp$Response)
set.seed(12)
index <- createDataPartition(comb_resp$Response, p = 0.7, list = FALSE)
train_data <- comb_resp[index, ]
test_data  <- comb_resp[-index, ]


## Positive peaks : Pos_X37  + Pos_X157 + Pos_X249 + Pos_X322 + Pos_X414 + Pos_X932 + Pos_X1137
## Negative peaks: Neg_X169 + Neg_X172 + Neg_X173 + Neg_X182 + Neg_X189 + Neg_X191 + Neg_X351 + Neg_X589 + Neg_X692 + Neg_X711 + Neg_X964,
model_pos <- glm(Response ~ 
               Pos_X37 +  Pos_X322 + 
                Pos_X1137,
                 family=binomial(link='logit'),data=train_data)

model_neg <- glm(Response ~ 
                   Neg_X692 + Neg_X711 + Neg_X964 + Neg_X1028,
                 family=binomial(link='logit'),data=train_data)

model_old <- glm(Response ~    Pos_X1137 +
                Neg_X711 + Neg_X964 ,
                 family=binomial(link='logit'),data=train_data)

model <- glm(Response ~      
                  Pos_X249 +Pos_X37 +
               Neg_X692 +  Neg_X711,
             family=binomial(link='logit'),data=train_data)


summary(model)
anova(model, test="Chisq")
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


# Logistic regression of metabolites for predicting synovitis
pos_resp <- read.csv('20210223_AF_pos_syno_resp.csv', row.names=1)
neg_resp <- read.csv('20210223_AF_neg_syno_resp.csv', row.names=1)

summary(pos_resp$Response == neg_resp$Response)

comb_resp <- cbind.data.frame(neg_resp, pos_resp[,-1])
comb_resp$Response <- as.factor(comb_resp$Response)
set.seed(42)
index <- createDataPartition(comb_resp$Response, p = 0.72, list = FALSE)
train_data <- comb_resp[index, ]
test_data  <- comb_resp[-index, ]


model <- glm(Response ~ Pos_X57  + Pos_X304 +  Pos_X396 + Pos_X873+
                Neg_X74 + Neg_X270 + Neg_X738 + Neg_X1028,
             family=binomial(link='logit'),data=train_data)

summary(model)
anova(model, test="Chisq")
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
