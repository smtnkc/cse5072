# Returns a balanced test sample
GetTestSample <- function(seed.val) {
  set.seed(seed.val)
  v.test.sample <- sample(x = c(1:9), size = 2)
  v.test.sample <- c(v.test.sample, sample(x = c(10:15), size = 1))
  v.test.sample <- c(v.test.sample, sample(x = c(16:21), size = 1))
  v.test.sample <- c(v.test.sample, sample(x = c(22:27), size = 1))
  v.test.sample <- c(v.test.sample, sample(x = c(28:35), size = 2))
  return (v.test.sample)
}

# Returns a df of success and execution time of a model on various num of genes
GetSuccessRates <- function(df, list.fc, type, tool, seed.val) {
  
  if(type == "TOP") {
    range <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 
                   100, 120, 150, 200, 250, 500, 1000, 2000, 5000, 10000, 24283)
    cat("n\tTrain\tTest\tTime (sec)\n")
  }
  else if(type == "CUT") {
    range <- c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5)
    cat("k\tTrain\tTest\tTime (sec)\n")
  }
  
  df.stats <- data.frame(n = numeric(length(range)),
                         train = numeric(length(range)),
                         test = numeric(length(range)),
                         time = numeric(length(range)))
  i <- 1
  for(n in range) {
    
    if(type == "TOP")
      df.deg <- GetDfDeg(df, list.fc, cut = NULL, top = n, verbose = 1)
    else if(type == "CUT")
      df.deg <- GetDfDeg(df, list.fc, cut = n, top = NULL, verbose = 1)
    
    v.test.sample <- GetTestSample(seed.val)
    df.train <- df.deg[-v.test.sample, ] 
    df.test <- df.deg[v.test.sample, ]
    
    if(tool == "c50") {
      time.begin <- Sys.time()
      model <- C5.0(df.train[, -1], df.train[, 1])
      time.end <- Sys.time()
    }
    else if(tool == "svm") {
      time.begin <- Sys.time()
      model <- svm(df.train$Group~., df.train)
      time.end <- Sys.time()
    }
    else stop("Enter a valid parameter for tool!")
    
    pred.train <- predict(model, df.train[, -1])
    tab.train <- table(pred.train, df.train[,1])
    success.train <- round(sum(diag(tab.train))/sum(tab.train)*100, 2)
    
    pred.test <- predict(model, df.test[, -1])
    tab.test <- table(pred.test, df.test[, 1])
    success.test <- round(sum(diag(tab.test))/sum(tab.test)*100,2)
    
    df.stats[i, "n"] <- n
    df.stats[i, "train"] <- success.train
    df.stats[i, "test"] <- success.test
    df.stats[i, "time"] <- time.end - time.begin
    
    i <- i + 1
    
    cat(n, "\t", success.train, "\t", 
        success.test, "\t",
        round((time.end - time.begin), 5), "\n")
  }
  return (df.stats)
}

# Returns a simple df which includes only the genes those are splitting attrs
GetDfSplittingAttrs <- function(df, verbose) {
  
  df.deg <- df
  time.begin <- Sys.time()
  model <- C5.0(df.deg[, -1], df.deg[, 1])
  pred <- predict(model, df.deg[, -1])
  time.end <- Sys.time()
  tab <- table(pred, df.deg[, 1])
  success <- round(sum(diag(tab))/sum(tab)*100, 2)
  
  if(verbose > 1) { 
    print(model); 
    print(summary(model)); 
  }
  
  df.imp <- C5imp(model, metric = "usage")  
  df.splitting.attrs <- df.imp[df.imp[,1] > 0, , drop = FALSE]
  v.splitting.attrs <- row.names(df.splitting.attrs)
  
  if(verbose > 0) {
    cat("Splitting Attrs:", v.splitting.attrs, "\n")
    cat("Time (sec):", time.end - time.begin, "\n")
  }
  
  return (df.deg[, v.splitting.attrs])
}

# Returns a vector of success ratios for a model on train and test samples
TT <- function(df, tool, seed.val, verbose) {

  v.test.sample <- GetTestSample(seed.val)
  df.train <- df[-v.test.sample, ] 
  df.test <- df[v.test.sample, ]
  
  if(tool == "c50") 
    model <- C5.0(df.train[, -1], df.train[, 1])
  else if(tool == "svm") 
    model <- svm(df.train$Group~., df.train)
  
  pred.train <- predict(model, df.train[, -1])
  tab.train <- table(pred.train, df.train[,1])
  success.train <- round(sum(diag(tab.train))/sum(tab.train)*100, 2)
  
  pred.test <- predict(model, df.test[, -1])
  tab.test <- table(pred.test, df.test[, 1])
  success.test <- round(sum(diag(tab.test))/sum(tab.test)*100, 2)
  
  if(verbose > 2) {
    print(model)
    print(summary(model))
    print(tab.train)
    cat("Success on train set: ", success.train, "%\n")
    print(tab.test)
    cat("Success on test set: ", success.test, "%\n")
  }
  
  v.success.rates <- c(success.train, success.test)
  return(v.success.rates)
}

# Runs TT function n-times for different seed vals starting from seed.begin
# Returns nothing but prints success ratio and running time for each iteration.
TTRunner <- function(df, tool, seed.begin, n.times, verbose) {
  train.sum <- 0
  test.sum <- 0
  total.time <- 0
  
  if(verbose > 1) cat("Seed")
  if(verbose > 0) cat("\tTrain\tTest\tTime (sec)\n")

  for(i in c(seed.begin:(seed.begin + n.times - 1))) {
    time.begin <- Sys.time()
    v.success.rates <- TT(df, tool, seed.val = i, verbose)
    time.end <- Sys.time()
    time.diff <- time.end - time.begin
    total.time <- total.time + time.diff
    
    train.sum <- train.sum + v.success.rates[1]
    test.sum <- test.sum + v.success.rates[2]
    
    if(verbose > 1) 
      cat(i,"\t", v.success.rates[1], 
          "\t", v.success.rates[2], 
          "\t", round(time.diff,5), "\n")
  }

  if(verbose > 1)
    cat("----------------------------------\n")
  if(verbose > 0) {
    cat("Avg:\t", round(train.sum/n.times, 2), 
      "\t", round(test.sum/n.times,2), 
      "\t", round(total.time/n.times,5), "\n")
  }
}