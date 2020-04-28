###采用十折交叉验证方法进行逻辑回归
library(caret)
cross_val <- function(df, result_name)
{
  #设定随机数种子
  set.seed(100)
  
  names(df)[1] <- "Y"
  folds <- createFolds(y=df$Y,k=10)
  
  max=0
  num=0
  
  for(i in 1:10)
  {
    fold_test <- df[folds[[i]],]   
    fold_train <- df[-folds[[i]],]
    
    print("*******Number*******")
    #模型拟合
    fold_pre <- glm(Y ~., family=binomial(link='logit'), data=fold_train)
    fold_predict <- predict(fold_pre, type='response', newdata=fold_test)
    
    #计算测试集精确度
    fold_predict <- ifelse(fold_predict>0.5,1,0)
    fold_test$predict <- fold_predict
    fold_error = fold_test$predict - fold_test$Y
    fold_accuracy = (nrow(fold_test)-sum(abs(fold_error)))/nrow(fold_test)
    print(i)
    print("***test_df_accuracy***")
    print(fold_accuracy)
    
    #计算训练集精确度
    fold_predict2 <- predict(fold_pre,type='response',newdata=fold_train)
    fold_predict2 =ifelse(fold_predict2>0.5,1,0)
    fold_train$predict = fold_predict2
    fold_error2 = fold_train$predict - fold_train$Y
    fold_accuracy2 = (nrow(fold_train)-sum(abs(fold_error2)))/nrow(fold_train) 
    print("***train_df_accuracy***")
    print(fold_accuracy2)
    
    
    if(fold_accuracy>max)
    {
      max=fold_accuracy  
      num=i
    }	
  }
  
  ###输出效果最好的一组数据
  print("***Best_result: ***")
  print(max)
  print(num)
  
  
  ###使用最好的结果
  testi <- df[folds[[num]],]
  traini <- df[-folds[[num]],]
  
  prei <- glm(Y ~., family=binomial(link='logit'), data=traini)
  predicti. <- predict.glm(prei,type='response',newdata=testi)
  
  predicti =ifelse(predicti. >0.5,1,0)
  testi$predict = predicti
  errori = testi$predict-testi$Y
  accuracyi = (nrow(testi)-sum(abs(errori)))/nrow(testi) 
  print("****Best_test_df_accuracy:**")
  print(accuracyi)
  
  
  ###画ROC曲线
  library(pROC)
  true_value <- testi$Y
  pdf(paste0("roc_", sub("\\.RData", "", result_name), ".pdf"), width = 7, height = 7)
  plot.roc(true_value, predicti., legacy.axes=TRUE, grid=c(0.1, 0.2),
           print.auc=TRUE, max.auc.polygon=TRUE, print.thres= F, ci = TRUE) 
  dev.off()
  ###保存训练模型
  prei$data = NULL  ## Reduce model size
  save(prei, file = result_name)
}


###逻辑回归数据整理函数
var_scale <- function(df)
{
  df_names <- c("donor","chr","start","end","res","base","nbase","reptime","tfbs","conser","gc","pro","cpg")
  if(!all(names(df)[1:13] == df_names))  break("Error!")
  df$reptime <- scale(df$reptime, center = F)
  df$tfbs <- scale(df$tfbs, center = F)
  df$gc <- scale(df$gc, center = F)
  df$nbase <- as.character(df$nbase)
  return(df)
}

samp_neg_df <- function(df)
{
  df_pos <- df[df$res == 1, ]
  df_neg <- df[df$res == 0, ]
  pos_num <- dim(df_pos)[1]
  neg_num <- dim(df_neg)[1]
  neg_index <- sample(1:neg_num, pos_num)
  samp_neg <- df_neg[neg_index,]
  rbind(df_pos, samp_neg)
}

logit_form <- function(project,file_name, donor = "~/paper/lasso/id_pro_index.tsv")
{	
  ###匹配项目信息
  id_index <- fread(donor)
  mock <- data.frame("icgc_donor_id" = "mock", "project_code" ="mock")
  id_index <- rbind(id_index, mock)
  result <- fread(file_name)
  merge_result <- merge(result, id_index, by.x = "V4", by.y = "icgc_donor_id", all.x = T)
  
  ###选取特定肿瘤
  mut <- merge_result[merge_result$project_code %in% project,]
  
  ###数据归一化
  names(mut)[1:13] <- c("donor","chr","start","end","res","base","nbase","reptime","tfbs","conser","gc","pro","cpg")
  mut <- var_scale(mut)
  
  ###选取等量的对照数据
  mut <- samp_neg_df(mut)
  
  return(as.data.frame(mut))
}


pos_anno <- function(pos_file, anno_file, result_file)
{
  temp1 <- paste("sorted_", basename(pos_file), sep = "")
  temp2 <- paste("sorted_", basename(anno_file), sep = "")
  cmd1 <- paste("sort -k 1,1 -k 2,2n ",pos_file," > ",temp1, sep = "")
  cmd2 <- paste("sort -k 1,1 -k 2,2n ",anno_file," > ",temp2, sep = "")
  system(cmd1)
  system(cmd2)
  cmd3 <- paste("bedtools map -a ",temp1," -b ",temp2," -c 4 -o mean > ",result_file, sep = "")
  system(cmd3)
  file.remove(temp1, temp2)
}




###添加DNase函数
add_dnase <- function(df,file_names)
{
  df$project_code = NULL
  df <- df[,c(2:ncol(df), 1)]
  df_names = colnames(df)
  data.table::fwrite(df, file = "pos.bed", 
                     sep = "\t", col.names = F)
  rm(df); gc()
  pos_anno("pos.bed", file_names, "raw_add.bed")
  
  raw_add <- fread("raw_add.bed", stringsAsFactors = F, data.table = F, header = FALSE)
  col_n <- ncol(raw_add)
  new_value <- raw_add[,col_n]
  raw_add[,col_n][new_value == "."] <- 0
  names(raw_add)[1:length(df_names)] <- df_names
  names(raw_add)[col_n] <- "dnase"
  raw_add$dnase <- as.numeric(raw_add$dnase)
  raw_add$nbase <- as.character(raw_add$nbase)
  
  file.remove("pos.bed", "raw_add.bed")
  return(raw_add)
}


###多重共线性检验
mult_validate <- function(df)
{
  anno <- c("reptime","tfbs","conser","gc","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","H3K9ac","dnase")  
  df <- df[,names(df) %in% anno]
  x <- cor(df)
  
  print("***multiple validata is***")
  print(kappa(x, exact = TRUE))
}

###类型检验
type_validate <- function(df)
{
  print("is the base type character:")
  if(is.character(df$nbase))
  {	
    print("TRUE")
  }else{   
    print("FALSE")
  }
  
  anno <- c("reptime","tfbs","conser","gc","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","H3K9ac","dnase")  
  
  print("is other type is numeric:")
  if(is.numeric(as.matrix(df[,names(df) %in% anno])))
  {
    print("TRUE")
  }else{
    print("FALSE")
  }
}


###整理格式函数
format_trans <- function(df)
{
  anno <- c("res","nbase","pro","cpg","reptime","tfbs","conser","gc","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","H3K9ac","dnase")  
  return(df[ ,names(df) %in% anno])
}





