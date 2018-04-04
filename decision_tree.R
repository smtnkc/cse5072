library(readr)
library(rstudioapi)
library(C50)
library(e1071)
rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
source("preprocess.R")
source("deg_detect.R")

###########################

df_deg <- get_df_deg(df_exp, fc_list, cut_k = 3, top_n = NULL)
df_deg <- get_df_deg(df_exp, fc_list, cut_k = NULL, top_n = 20)

# TRAIN VS TEST
tt <- function(df, print_table, seed_val, tool) {
  set.seed(seed_val) # Change seed to differentiate the samples
  
  s <- sample(nrow(df), 0.7*nrow(df))
  df_train <- df[s,]; 
  df_test <- df[-s,]
  if(tool == "c50") model <- C5.0(df_train[, -1], df_train[, 1])
  else if(tool == "svm") model <- svm(df_train$Group~., df_train)
  
  if(print_table) { print(model); print(summary(model)) }
  
  trainPred <- predict(model, df_train[, -1])
  tab_train <- table(trainPred, df_train[,1])
  if(print_table) print(tab_train)
  train_res <- round(sum(diag(tab_train))/sum(tab_train)*100,2)
  if(print_table) cat("Success on train set: ", train_res, "%\n")
  
  testPred <- predict(model, df_test[, -1])
  tab_test <- table(testPred, df_test[,1])
  if(print_table) print(tab_test)
  test_res <- round(sum(diag(tab_test))/sum(tab_test)*100,2)
  if(print_table) cat("Success on test set: ", test_res, "%\n")
  
  res <- c(train_res, test_res)
  return(res)
}
tt_runner <- tt_runner <- function(df, print_table, seed_begin, nTimes, tool) {
  train_sum <- 0
  test_sum <- 0
  time <- 0
  
  for(i in c(seed_begin:(seed_begin+nTimes-1))) {
    start_time <- Sys.time()
    res <- tt(df, print_table = print_table, i, tool)
    end_time <- Sys.time()
    time_diff <- end_time-start_time
    time <- time + time_diff
    
    train_sum <- train_sum + res[1]
    test_sum <- test_sum + res[2]
    cat("Seed:", i, 
        "\tTrain:", res[1], 
        "\tTest:", res[2], 
        "\tTime:", time_diff, "\n")
  }
  cat("----------------------------------\n")
  cat("Train Avg:", train_sum/nTimes, 
      "\tTest Avg:", test_sum/nTimes, 
      "\tTotal Time:", time, "\n")
}

tt_runner(df_deg, print_table = FALSE, seed_begin = 2000, nTimes = 100, tool = "svm")
tt_runner(df_deg, print_table = FALSE, seed_begin = 2000, nTimes = 100, tool = "c50")

tt_runner(df_all, print_table = FALSE, seed_begin = 5000, nTimes = 1, tool = "c50")

# ALL VS DEG RUNNER
comparative_runner <- function(df1, df2, print_table, nTimes) {
  train_sum1 <- 0
  test_sum1 <- 0
  train_sum2 <- 0
  test_sum2 <- 0
  time1 <- 0
  time2 <- 0
  
  for(i in sample(1000:2000, nTimes, replace = FALSE)) {
    start_time <- Sys.time()
    res1 <- tt(df1, print_table = print_table, i)
    end_time <- Sys.time()
    time_diff <- end_time-start_time
    time1 <- time1 + time_diff
    
    train_sum1 <- train_sum1 + res1[1]
    test_sum1 <- test_sum1 + res1[2]
    cat("Naive -> Seed:", i, 
        "\tTrain:", res1[1], 
        "\tTest:", res1[2], 
        "\tTime:", time_diff, "\n")
    
    start_time <- Sys.time()
    res2 <- tt(df2, print_table = print_table, i)
    end_time <- Sys.time()
    time_diff <- end_time-start_time
    time2 <- time2 + time_diff
    
    train_sum2 <- train_sum2 + res2[1]
    test_sum2 <- test_sum2 + res2[2]
    cat("DEGs  -> Seed:", i, 
        "\tTrain:", res2[1],
        "\tTest:", res2[2], 
        "\tTime:", time_diff, "\n")
  }
  
  cat("----------------------------------\n")
  cat("Naive AVG -> Train:", train_sum1/nTimes, 
      "\tTest:", test_sum1/nTimes,
      "\tTotal Time:", time1, "\n")
  cat("DEGs  AVG -> Train:", train_sum2/nTimes, 
      "\tTest:", test_sum2/nTimes, 
      "\tTotal Time:", time2, "\n")
}
comparative_runner(df_all, df_deg, print_table = FALSE, 10)

# source("simplify_df.R")