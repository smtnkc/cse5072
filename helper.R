# Find success ratio and execution time of a model when applied to the top DEG genes
# Models are run for hole set without train-test seperation
# Returns nothing but prints success ratio and execution time for each iteration
SuccessFinder <- function(df, fc.list, type, tool) {
  
  if(type == "TOP") {
    range <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 250, 500)
    cat("n\tSuccess\tTime (sec)\n")
  }
  else if(type == "CUT") {
    range <- c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5)
    cat("k\tSuccess\tTime (sec)\n")
  }
  
  for(n in range) {
    
    if(type == "TOP")
      df.deg <- GetDfDeg(df, fc.list, cut = NULL, top = n, verbose = 1)
    else if(type == "CUT")
      df.deg <- GetDfDeg(df, fc.list, cut = n, top = NULL, verbose = 1)
    
    time.begin <- Sys.time()
    
    if(tool == "c50")
      model <- C5.0(df.deg[, -1], df.deg[, 1])
    else if(tool == "svm")
      model <- svm(df.deg$Group~., df.deg)
    
    pred <- predict(model, df.deg[, -1])
    time.end <- Sys.time()
    tab <- table(pred, df.deg[, 1])
    success <- round(sum(diag(tab))/sum(tab)*100, 2)
    
    cat(n, "\t", success, "\t", round((time.end - time.begin) / length(range), 5), "\n")
  }
}

# Creates a model ...
# Returns a simple df which includes only the genes those have been used as splitting attributes
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
  
  if(verbose > 0) 
    print(v.splitting.attrs)
  
  return (df.deg[, v.splitting.attrs])
}

# Splits data frame to train and test sets and creates a model by using the train set.
# Applies the model to both train and test sets. 
# Returns a vector of success ratio for both train and test sets.
TT <- function(df, train.perc, tool, seed.val, verbose) {
  set.seed(seed.val) # Change seed to differentiate the samples
  
  v.sample <- sample(nrow(df), train.perc*nrow(df))
  df.train <- df[v.sample,] 
  df.test <- df[-v.sample,]

  if(tool == "c50") 
    model <- C5.0(df.train[, -1], df.train[, 1])
  else if(tool == "svm") 
    model <- svm(df.train$Group~., df.train)
  
  if(verbose > 2) { 
    print(model)
    print(summary(model)) 
  }
  
  pred.train <- predict(model, df.train[, -1])
  tab.train <- table(pred.train, df.train[,1])
  
  if(verbose > 2) 
    print(tab.train)
  
  success.train <- round(sum(diag(tab.train))/sum(tab.train)*100, 2)

  if(verbose > 2)
    cat("Success on train set: ", success.train, "%\n")
  
  pred.test <- predict(model, df.test[, -1])
  tab.test <- table(pred.test, df.test[, 1])
  
  if(verbose > 2) 
    print(tab.test)
  
  success.test <- round(sum(diag(tab.test))/sum(tab.test)*100,2)
  
  if(verbose > 2) 
    cat("Success on test set: ", success.test, "%\n")
  
  v.success.rates <- c(success.train, success.test)
  return(v.success.rates)
}

# Runs "tt" function n-times for the seeds [seed.begin, seed.begin + n]
# Returns nothing but prints success ratio and running time for each iteration.
TTRunner <- function(df, train.perc, tool, seed.begin, n.times, verbose) {
  train.sum <- 0
  test.sum <- 0
  total.time <- 0
  
  if(verbose > 1) cat("Seed")
  if(verbose > 0) cat("\tTrain\tTest\tTime (sec)\n")

  for(i in c(seed.begin:(seed.begin + n.times - 1))) {
    time.begin <- Sys.time()
    v.success.rates <- TT(df, train.perc, tool, i, verbose)
    time.end <- Sys.time()
    time.diff <- time.end - time.begin
    total.time <- total.time + time.diff
    
    train.sum <- train.sum + v.success.rates[1]
    test.sum <- test.sum + v.success.rates[2]
    if(verbose > 1) cat(i,"\t", v.success.rates[1], "\t", v.success.rates[2], "\t", round(time.diff,5), "\n")
  }

  if(verbose > 1)
    cat("----------------------------------\n")
  if(verbose > 0) {
    cat("Avg:\t", round(train.sum/n.times, 2), 
      "\t", round(test.sum/n.times,2), 
      "\t", round(total.time/n.times,5), "\n")
  }
}