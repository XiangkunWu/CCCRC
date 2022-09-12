#目的：特征筛选过程：差异基因（分型差异基因/PDX差异记忆你） + CCLE肿瘤特征基因 +  Boruta特征筛选 + LASSO特征筛选
#版本：gemh-22.06.23
#备注：Boruta算法/LASSO算法的输入数据均为GEO表达谱；阈值均可根据实际情况做调参
library(data.table)
library(tidyverse)
library(Boruta)
library(glmnet)
options(stringsAsFactors = FALSE)
#数据读取-------------------------------
TCGA_cluster <- fread("pamr_TCGA_CPTAC.txt",header = T)
TCGA_DE <- rbind(fread("limma_test_result.C1_unique_upexpr_marker.txt") %>%
                   mutate(Compare_Group = "C1_C234"),
                 fread("limma_test_result.C2_unique_upexpr_marker.txt") %>%
                   mutate(Compare_Group = "C2_C134"),
                 fread("limma_test_result.C3_unique_upexpr_marker.txt") %>%
                   mutate(Compare_Group = "C3_C124"),
                 fread("limma_test_result.C4_unique_upexpr_marker.txt") %>%
                   mutate(Compare_Group = "C4_C123"))
TCGA_data <- fread("sva_voom_batch2_T_protein_coding.csv")

feature_Gene <- TCGA_DE$V1
length(feature_Gene)
#Boruta------------------
set.seed(1)
Boruta_input <- TCGA_data %>% filter(V1 %in% feature_Gene) %>% 
  dplyr::select(V1,TCGA_cluster$id) %>%
  column_to_rownames("V1") %>% 
  t() %>% as.matrix()
borutafit <- Boruta(x = Boruta_input, 
                    y = as.factor(TCGA_cluster$IMcluster), # multiclassification
                    doTrace = 2,
                    maxRuns = 50,
                    ntree = 500)
boruta_fea <- attStats(borutafit)
boruta_fea <- rownames(boruta_fea[which(boruta_fea$decision == "Confirmed"),])
#LASSO------------------
LASSO_input <- TCGA_data %>% filter(V1 %in% feature_Gene) %>% 
  dplyr::select(V1,TCGA_cluster$id) %>%
  column_to_rownames("V1") %>% 
  t() %>% as.matrix()
cv.fit <- cv.glmnet(
  x = LASSO_input,
  y = as.factor(TCGA_cluster$IMcluster),
  family = "multinomial",
  alpha = 1,
  nlambda = 100,
  nfolds = 10,
  parallel=TRUE
)
lasso_best <- glmnet(
  x =  LASSO_input,
  y =  as.factor(TCGA_cluster$IMcluster),
  family="multinomial",
  lambda  = cv.fit$lambda.min
)
glmnet_result <- coef(lasso_best)
glmnet_fea <- lapply(1:4, function(i){
  (glmnet_result[[i]] %>%
     as.matrix() %>% data.frame() %>%
     filter(s0 >0) %>% rownames())[-1] %>% return()
}) %>% unlist() %>% unique()
#特征基因筛选结果---------------------------
feature_genelist <- boruta_fea
write.table(feature_Gene,"feature787_GEO.txt",sep = "\t",row.names = F,quote = F)
#存储分析过程-------------------------------
save(feature_genelist,file = "feature_genelist.Rdata")
save(TCGA_GEO_overlap_Genelist,GEO_DEGenelist,PDX_DEGenelist,ccle_filter_1,ccle_filter_2,
     borutafit,boruta_fea,
     cv.fit,lasso_best,glmnet_fea,
     file = "C:/Users/noenemy/Desktop/单基因分类器/feature_process.Rdata")
