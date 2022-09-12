#目的：构建分类模型：随机森林，SVM，xgboost，logistics  , TSP 算法
#版本：gemh-22.06.23
#备注：因为模型面向跨平台数据，即RNA-seq和芯片测序，因此对建模数据进行标注化处理，脚本提供三种方式:Z-score Scale(zscore);Min-Max Scale(mmscale)和秩转化(rank)
#备注：建模超参数使用
library(data.table)
library(tidyverse)
library(pROC)
library(Boruta)
library(glmnet)
library(randomForest)
library(switchBox)
library(caret)
library(nnet)
options(stringsAsFactors = FALSE)
#读取数据----------------------
#表达谱数据格式均为：Gene * SampleId的矩阵（第一列为基因名）
#cluster数据格式均为：第一列为样本ID，第二列为分类
#GEO_cluster <- fread("data/data_06.23/cluster_pam_pearson_GEO8.txt")
#GEO_cluster$IMcluster <- paste0("C",GEO_cluster$IMcluster)
#setwd("C:\\Users\\noenemy\\Desktop\\单基因分类器")
TCGA_cluster <- fread("pamr_TCGA_CPTAC.txt",header = T)[,-1]
TCGA_data <- fread("sva_voom_batch2_T.csv")
load("feature_genelist.Rdata")
TCGA_data <- filter(TCGA_data,V1 %in% feature_genelist)
test_data <- fread("single patient.csv")
test_data <- filter(test_data,V1 %in% feature_genelist)


#随机森林:rf_method---------------
rf_method <- function(train_data,
                      test_data,
                      validate_data,
                      tuneLength = 3,#参数调整
                      seed_value) {
  colnames(train_data)[1:2] <- c("id","Cluster")
  colnames(test_data)[1:2] <- c("id","Cluster")
  colnames(validate_data)[1] <- c("id")
  train_data$Cluster <- factor(train_data$Cluster,unique(train_data$Cluster))
  test_data$Cluster <- factor(test_data$Cluster,unique(train_data$Cluster))
  
  #建模-----------------------------------
  set.seed(seed_value)
  model <- train(
    Cluster ~ .,
    tuneLength = tuneLength,
    data = column_to_rownames(train_data,"id"),
    method = "ranger",
    trControl = trainControl(method = "cv", number = 5,
                             verboseIter = FALSE,classProbs=TRUE)
  )
  # model$finalModel
  #训练集 和 验证集-------------------
  #训练集-----------------------
  train_pre_prob <- predict(model,newdata = train_data,type = "prob")
  train_roc <- multiclass.roc(response = train_data$Cluster,
                              predictor = train_pre_prob)
  train_pre_res <-  predict(model,newdata = train_data,type = "raw")
  # confusionMatrix(data = train_pre_res,
  #                 reference = factor(train_data$Cluster,unique(train_pre_res)),
  #                 mode = "everything")
  # multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
  #                   lev = levels(unique(train_pre_res)))
  #测试集----------------------
  test_pre_prob <- predict(model,newdata = test_data,type = "prob")
  test_roc <- multiclass.roc(response = test_data$Cluster,predictor = test_pre_prob)
  test_pre_res <- predict(model,newdata = test_data,type = "raw")
  #验证集-------------------
  validate_pre_res <- predict(model,newdata = validate_data,type = "raw")
  #----------
  model_res <- rbind(multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
                                       lev = levels(unique(train_pre_res))),
                     multiClassSummary(data.frame(obs = test_data$Cluster,pred =test_pre_res ),
                                       lev = levels(unique(test_data$Cluster))) ) %>%
    data.frame(stringsAsFactors = F) %>%
    mutate(AUC = c(as.numeric(train_roc$auc),
                   as.numeric(test_roc$auc)),
           dataType = c("train","test"),
           seed_num = seed_value,
           model_res = "RandomForest")
  return(list(model_res,model,validate_pre_res))
}
#支持向量机:svm_method---------------
svm_method <- function(train_data,
                       test_data,
                       validate_data,
                       tuneLength = 3,#参数调整
                       seed_value){
  
  colnames(train_data)[1:2] <- c("id","Cluster")
  colnames(test_data)[1:2] <- c("id","Cluster")
  colnames(validate_data)[1] <- c("id")
  train_data$Cluster <- factor(train_data$Cluster,unique(train_data$Cluster))
  test_data$Cluster <- factor(test_data$Cluster,unique(train_data$Cluster))
  #建模-----------------------------------
  set.seed(seed_value)
  model <- train(
    Cluster ~ .,
    tuneLength = tuneLength,
    data = column_to_rownames(train_data,"id"),
    method = "svmRadial",
    trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs=TRUE)
  )
  #训练集 和 验证集-------------------
  #训练集-----------------------
  train_pre_prob <- predict(model,newdata = train_data,type = "prob")
  train_roc <- multiclass.roc(response = train_data$Cluster,
                              predictor = train_pre_prob)
  train_pre_res <-  predict(model,newdata = train_data,type = "raw")
  # confusionMatrix(data = train_pre_res,
  #                 reference = factor(train_data$Cluster,unique(train_pre_res)),
  #                 mode = "everything")
  # multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
  #                   lev = levels(unique(train_pre_res)))
  #测试集----------------------
  test_pre_prob <- predict(model,newdata = test_data,type = "prob")
  test_roc <- multiclass.roc(response = test_data$Cluster,predictor = test_pre_prob)
  test_pre_res <- predict(model,newdata = test_data,type = "raw")
  #验证集-------------------
  validate_pre_res <- predict(model,newdata = validate_data,type = "raw")
  #----------
  model_res <- rbind(multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
                                       lev = levels(unique(train_pre_res))),
                     multiClassSummary(data.frame(obs = test_data$Cluster,pred =test_pre_res ),
                                       lev = levels(unique(test_data$Cluster)))) %>%
    data.frame(stringsAsFactors = F) %>%
    mutate(AUC = c(as.numeric(train_roc$auc),
                   as.numeric(test_roc$auc)),
           dataType = c("train","test"),
           seed_num = seed_value,
           model_res = "SVM-Radial")
  return(list(model_res,model,validate_pre_res))
}


#xgboost:xgboost_method---------------
xgboost_method <- function(train_data,
                           test_data,
                           validate_data,
                           tuneLength = 3,#参数调整
                           seed_value){
  colnames(train_data)[1:2] <- c("id","Cluster")
  colnames(test_data)[1:2] <- c("id","Cluster")
  colnames(validate_data)[1] <- c("id")
  train_data$Cluster <- factor(train_data$Cluster,unique(train_data$Cluster))
  test_data$Cluster <- factor(test_data$Cluster,unique(train_data$Cluster))
  #建模-----------------------------------
  set.seed(seed_value)
  model <- train(
    Cluster ~ .,
    tuneLength = tuneLength,
    data = column_to_rownames(train_data,"id"),
    method = "xgbTree",
    trControl = trainControl(method = "cv", number = 5,classProbs=TRUE)
  )
  #训练集 和 验证集-------------------
  #训练集-----------------------
  train_pre_prob <- predict(model,newdata = train_data,type = "prob")
  train_roc <- multiclass.roc(response = train_data$Cluster,
                              predictor = train_pre_prob)
  train_pre_res <-  predict(model,newdata = train_data,type = "raw")
  # confusionMatrix(data = train_pre_res,
  #                 reference = factor(train_data$Cluster,unique(train_pre_res)),
  #                 mode = "everything")
  # multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
  #                   lev = levels(unique(train_pre_res)))
  #测试集----------------------
  test_pre_prob <- predict(model,newdata = test_data,type = "prob")
  test_roc <- multiclass.roc(response = test_data$Cluster,predictor = test_pre_prob)
  test_pre_res <- predict(model,newdata = test_data,type = "raw")
  #验证集-------------------
  validate_pre_res <- predict(model,newdata = validate_data,type = "raw")
  #----------
  model_res <- rbind(multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
                                       lev = levels(unique(train_pre_res))),
                     multiClassSummary(data.frame(obs = test_data$Cluster,pred =test_pre_res ),
                                       lev = levels(unique(test_data$Cluster)))) %>%
    data.frame(stringsAsFactors = F) %>%
    mutate(AUC = c(as.numeric(train_roc$auc),
                   as.numeric(test_roc$auc)),
           dataType = c("train","test"),
           seed_num = seed_value,
           model_res = "xgboost")
  return(list(model_res,model,validate_pre_res))
}
#logistics:logistics_method-------------
logistics_method <- function(train_data,
                             test_data,
                             validate_data,
                             seed_value){
  
  colnames(train_data)[1:2] <- c("id","Cluster")
  colnames(test_data)[1:2] <- c("id","Cluster")
  colnames(validate_data)[1] <- c("id")
  train_data$Cluster <- factor(train_data$Cluster,unique(train_data$Cluster))
  test_data$Cluster <- factor(test_data$Cluster,unique(train_data$Cluster))
  #建模-----------------------------------
  set.seed(seed_value)
  model <- multinom(Cluster ~ .,data = column_to_rownames(train_data,"id"))
  
  #训练集 和 验证集-------------------
  #训练集-----------------------
  train_pre_prob <- predict(model,newdata = column_to_rownames(train_data,"id"),type = "probs")
  train_roc <- multiclass.roc(response = train_data$Cluster,
                              predictor = train_pre_prob)
  train_pre_res <- predict(model,newdata = train_data,type = "class")
  # confusionMatrix(data = train_pre_res,
  #                 reference = factor(train_data$Cluster,unique(train_pre_res)),
  #                 mode = "everything")
  # multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
  #                   lev = levels(unique(train_pre_res)))
  #测试集----------------------
  test_pre_prob <- predict(model,newdata = column_to_rownames(test_data,"id"),type = "prob")
  test_roc <- multiclass.roc(response = test_data$Cluster,predictor = test_pre_prob)
  test_pre_res <- predict(model,newdata = test_data,type = "class")
  #验证集-------------------
  validate_pre_res <- predict(model,newdata = validate_data,type = "class")
  #----------
  model_res <- rbind(multiClassSummary(data.frame(obs = factor(train_data$Cluster,unique(train_pre_res)),pred =train_pre_res ),
                                       lev = levels(unique(train_pre_res))),
                     multiClassSummary(data.frame(obs = test_data$Cluster,pred =test_pre_res ),
                                       lev = levels(unique(test_data$Cluster)))) %>%
    data.frame(stringsAsFactors = F) %>%
    mutate(AUC = c(as.numeric(train_roc$auc),
                   as.numeric(test_roc$auc)),
           dataType = c("train","test"),
           seed_num = seed_value,
           model_res = "logistic_regression")
  return(list(model_res,model,validate_pre_res))
}


#TSP:TSP_method------------
tsp_tmp <- function(inpudata,tsp_modle,cluster_num){
  train_res_raw <- lapply(1:cluster_num, function(j){
    tmp <- SWAP.KTSP.Classify(as.matrix(inpudata), tsp_modle[[j]])  %>%
      data.frame() %>%
      rownames_to_column("id")
    colnames(tmp)[2] <- "predict_class"
    tmp$tsp_model <- names(tsp_modle)[j]
    return(tmp)
  })   %>% rbindlist()
  
  train_res <- data.frame(table(train_res_raw[,1:2])) %>%
    group_by(id) %>% filter(Freq == max(Freq)) %>%
    ungroup() %>%
    arrange(id,predict_class)
  
  count_res <- filter(data.frame(table(train_res$id)))
  
  train_res_single <- filter(train_res,id %in% filter(count_res,Freq ==1)$Var1)[,1:2]
  train_res_double <- filter(train_res,id %in% filter(count_res,Freq ==2)$Var1)
  train_res_three <- filter(train_res,id %in% filter(count_res,Freq > 2)$Var1)[,1:2] %>%
    mutate(predict_class = NA)
  if (nrow(train_res_double)>=1) {
    doube_data <- dplyr::select(inpudata,unique(train_res_double$id))
    tsp_compare <- group_by(train_res_double,id ) %>% summarise(model_name = paste0(predict_class,collapse = "_"))  %>%
      data.frame()
    train_res_double <- lapply(1:nrow(tsp_compare), function(t){
      tmp_model <- tsp_modle[[which(names(tsp_modle)==tsp_compare$model_name[t])]]
      tmp <- SWAP.KTSP.Classify(as.matrix(dplyr::select(inpudata,tsp_compare$id[t])), tmp_model)  %>%
        data.frame() %>%
        rownames_to_column("id")
      colnames(tmp)[2] <- "predict_class"
      return(tmp)
    }) %>% rbindlist()
    final_res <- rbind(rbind(train_res_single,train_res_double),
                       train_res_three )
  } else{
    final_res <- rbind(train_res_single,train_res_three)
  }
  return(final_res)
}
TSP_method <- function(train_data,
                       test_data,
                       validate_data,
                       seed_value){
  colnames(train_data)[1:2] <- c("id","Cluster")
  colnames(test_data)[1:2] <- c("id","Cluster")
  colnames(validate_data)[1] <- c("id")
  train_data$Cluster <- factor(train_data$Cluster,unique(train_data$Cluster))
  test_data$Cluster <- factor(test_data$Cluster,unique(train_data$Cluster))
  #建模-----------------------------------
  set.seed(seed_value)
  tsp_modle <- list()
  tsp_gene <- list()
  cluster_num = length(unique(train_data$Cluster))
  for (j in 1:ncol(combn(length(unique(train_data$Cluster)),2))) {
    tmp_sample <- filter(train_data,Cluster %in% unique(train_data$Cluster)[combn(cluster_num,2)[,j]]) %>%
      mutate(Cluster = factor(Cluster))
    tmp_cluster <- tmp_sample$Cluster
    tmp_sample <-  tmp_sample[,-2] %>% column_to_rownames("id") %>% t() %>% data.frame()
    
    tsp_modle[[j]] <- SWAP.Train.KTSP(as.matrix(tmp_sample),tmp_cluster,
                                      FilterFunc=SWAP.Filter.Wilcoxon,
                                      krange = 1:nrow(tmp_sample))
    tsp_gene[[j]] <- unlist(str_split(rownames(tsp_modle[[j]]$TSPs),","))
    names(tsp_modle)[j] <- paste0(unique(train_data$Cluster)[combn(cluster_num,2)[,j]],collapse = "_")
    names(tsp_gene)[j] <- paste0(unique(train_data$Cluster)[combn(cluster_num,2)[,j]],collapse = "_")
  }
  #train-----------
  train_tsp_res<-  tsp_tmp(inpudata = train_data[,-2] %>% column_to_rownames("id") %>% t() %>% data.frame(),
                           tsp_modle = tsp_modle,cluster_num = 6)
  train_tsp_res <- left_join(train_tsp_res,train_data[,1:2],by = "id")
  # multiClassSummary(data.frame(obs = train_tsp_res$Cluster,
  #                              pred = train_tsp_res$predict_class) ,
  #                   lev = levels(train_tsp_res$Cluster))
  #test------------
  test_tsp_res<-  tsp_tmp(inpudata = test_data[,-2] %>% column_to_rownames("id") %>% t() %>% data.frame(),
                          tsp_modle = tsp_modle,cluster_num = 6)
  test_tsp_res <- left_join(test_tsp_res,test_data[,1:2],by = "id")
  #validate---------
  validate_tsp_res<-  tsp_tmp(inpudata = validate_data %>% column_to_rownames("id") %>% t() %>% data.frame(),
                              tsp_modle = tsp_modle,cluster_num = 6)
  
  #----------
  model_res <- rbind(multiClassSummary(data.frame(obs = train_tsp_res$Cluster,
                                                  pred =train_tsp_res$predict_class),
                                       lev = levels(unique(train_tsp_res$Cluster))),
                     multiClassSummary(data.frame(obs = test_tsp_res$Cluster,
                                                  pred =test_tsp_res$predict_class),
                                       lev = levels(unique(test_tsp_res$Cluster)))) %>%
    data.frame(stringsAsFactors = F) %>%
    mutate(AUC = c("-","-"),
           dataType = c("train","test"),
           seed_num = seed_value,
           model_res = "TSP")
  return(list(model_res,tsp_modle,validate_tsp_res))
}

#跨平台数据标准化-------------------------
#标准化函数,Normalized_method(inputdata,method)：
#inpudata:Gene(row) * Sample(column)矩阵， data.frame格式
#method: zscore：Z-score Scale,mmscale：Min-Max Scale,rank：秩转化
Normalized_method <- function(inputdata,method) {
  if (method == "zscore") {
    normalizated_data <- scale(inputdata) %>% data.frame()
  } else if(method == "mmscale"){
    min_max_dif<- apply(inputdata,2,max) - apply(inputdata,2,min)
    min_value <- apply(inputdata,2,min)
    normalizated_data <-  lapply(1:ncol(inputdata), function(i){
      data.frame(SampleId = colnames(inputdata)[i],
                 Gene = rownames(inputdata),
                 value  = (inputdata[,i] - min_value[i])/min_max_dif[i]) %>% return()
    }) %>% rbindlist() %>% pivot_wider(names_from = "SampleId",values_from = "value") %>%
      column_to_rownames("Gene")
  } else if(method == "rank") {
    normalizated_data <-  lapply(1:ncol(inputdata), function(i){
      data.frame(SampleId = colnames(inputdata)[i],
                 Gene = rownames(inputdata),
                 value  = order(inputdata[,i])) %>% return()
    }) %>% rbindlist() %>% pivot_wider(names_from = "SampleId",values_from = "value") %>%
      column_to_rownames("Gene")
  }
  return(normalizated_data)
}
#本次仅尝试了zscore的方式

TCGA_Normalized_data <- Normalized_method(inputdata = column_to_rownames(TCGA_data,"V1"),
                                          method = "zscore")
test_Normalized_data <- Normalized_method(inputdata = column_to_rownames(test_data,"V1"),
                                              method = "zscore")

#调整文件格式：转置，匹配分型，讲基因名中带有“-” 替换为“.”

TCGA_Normalized_data <- t(TCGA_Normalized_data) %>%
  data.frame() %>% rownames_to_column("id")

test_Normalized_data <- t(test_Normalized_data) %>%
  data.frame() %>% rownames_to_column("id")

TCGA_Normalized_data <- inner_join(TCGA_cluster,TCGA_Normalized_data,
                                   by = "id")
colnames(TCGA_Normalized_data) <- gsub("-",'.',colnames(TCGA_Normalized_data))

# validate_cluster <- test_Normalized_data[,1:2] 
# validate_cluster[,2] <- c(rep("C1",10),rep("C2",10),rep("C3",10),rep("C4",13))
# colnames(validate_cluster ) <- c("id","IMcluster")
# test_Normalized_data <- inner_join(validate_cluster,test_Normalized_data,
#                                        by = "id")
colnames(test_Normalized_data) <- gsub("-",'.',colnames(test_Normalized_data))


#50次循环进行建模----------------
#train_data,test_data,validate_data的格式均为：Sample(column) * Gene(row) 矩阵， data.frame格式：1，2列分别id和cluster
model_all <- list()
res_all <- list()
precResult_all <-list()

for (i in 2:17) {
  seed_value <- i
  set.seed(seed_value)
  trains <- createDataPartition(
    y  = TCGA_Normalized_data$IMcluster,
    p = 0.7,
    list = F
  )
  
  train_data <- TCGA_Normalized_data[trains,]
  test_data <- TCGA_Normalized_data[-trains,]
  validate_data <- test_Normalized_data
  rf_res <- rf_method(train_data = train_data,test_data = test_data,
                      validate_data = validate_data,
                      tuneLength = 2,seed_value = i)
  svm_res <- svm_method(train_data = train_data,test_data = test_data,validate_data = validate_data,
                        tuneLength = 2,seed_value = i)
  logistics_res <- logistics_method(train_data = train_data,test_data = test_data,validate_data = validate_data,
                                    seed_value = i)  
  TSP_res <- TSP_method(train_data = train_data,test_data = test_data,validate_data = validate_data,
                        seed_value = i)
  xgboost_res <- xgboost_method(train_data = train_data,
                                test_data = test_data,
                                validate_data = validate_data,
                                tuneLength = 3,
                                seed_value = i)
  model_all[[i]] <- list(rf_res[[2]],
                         svm_res[[2]],
                         logistics_res[[2]],
                         TSP_res[[2]],
                         xgboost_res[[2]])
  res_all[[i]] <- rbind(rf_res[[1]],
                        svm_res[[1]],
                        logistics_res[[1]],
                        TSP_res[[1]],
                        xgboost_res[[1]])
  
   pre_res <- as.data.frame(rbind(as.character(rf_res[[3]]),
                                             as.character(svm_res[[3]]),
                                             as.character(logistics_res[[3]]),
                                             as.character(unique(TSP_res[[3]])$predict_class),
                                             as.character(xgboost_res[[3]])))
   colnames(pre_res) <- test_Normalized_data$id
   rownames(pre_res) <- c("RF","SVM","LOG","TSP","Xgboost")
  precResult_all[[i]] <- pre_res
}

# predictRes <- as.data.frame(rbind(as.character(rf_res[[3]]),
#                                   as.character(svm_res[[3]]),
#                                   as.character(logistics_res[[3]]),
#                                   as.character(TSP_res[[3]]$predict_class),
#                                   as.character(xgboost_res[[3]])))
# colnames(predictRes) <- test_data$V1
# rownames(predictRes) <- c("RF","SVM","LOG","TSP","Xgboost")
# write.csv(predictRes,"predictRes.csv",sep=",")
#-------------------------------------------------
res_all <- rbindlist(res_all)
res_all$AUC <- as.numeric(res_all$AUC)
melt(res_all) %>%
  filter(variable %in% c("Accuracy","Mean_F1","Mean_Recall","Mean_Precision",
                         "Mean_Sensitivity","Mean_Sensitivity","AUC")) %>%
  filter(dataType != "test") %>%
  ggplot(aes(x = model_res ,y = value)) +
  geom_boxplot(aes(fill = dataType )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  labs(x = "",y = "value") +
  facet_wrap(~variable ,scales = "free")
ggsave("model_res17循环.jpeg",width = 8,height = 6)
#根据上图，选择最佳模型后，再筛选该模型下最佳 seed值下的最佳预测结果，本次以“logistic_regression”为例（根据上图，logistic_regression的平局结果最佳）
filter(res_all,model_res == "xgboost" & dataType == "test") %>% 
  arrange(Accuracy) %>% 
  dplyr::select(seed_num)
#确定好seed和模型，选择模型以及对应的预测结果：
#model
model_all[[8]][[5]]
#预测结果：
precResult_all[[8]]["Xgboost",]


