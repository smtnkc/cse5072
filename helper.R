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

# Returns a simple df which includes only the genes those are splitting attrs
GetDfSplittingAttrs <- function(df, class.type, verbose) {
  
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
  
  df.splitting.attrs <- df.deg[, v.splitting.attrs]
  
  if(class.type == "num")
    df.splitting.attrs <- cbind(Group = df.states[,"Group"], df.splitting.attrs)
  else if(class.type == "char")
    df.splitting.attrs <- cbind(Group = df.states[,"Disease_State"], df.splitting.attrs)
  
  df.splitting.attrs$Group <- factor(df.splitting.attrs$Group)
  
  return (df.splitting.attrs)
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
  else if(tool == "rf")
    model <- randomForest(df.train$Group~., df.train)
  
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

# Runs TT function n-times for different seed values starting from seed.begin
# Returns a vector of stats (avg success and exec time of N iterations)
TTRunner <- function(df, cut, top, list.fc, tool, seed.begin, n.times, class.type, verbose) {
  
  if(is.null(cut) && is.null(top))
    df.deg <- df
  else
    df.deg <- GetDfDeg(df, list.fc, cut, top, class.type, verbose)
  
  train.sum <- 0
  test.sum <- 0
  total.time <- 0
  
  for(i in c(seed.begin:(seed.begin + n.times - 1))) {
    
    time.begin <- Sys.time()
    v.success.rates <- TT(df.deg, tool, seed.val = i, verbose)
    time.end <- Sys.time()
    time.diff <- time.end - time.begin
    total.time <- total.time + time.diff
    
    train.sum <- train.sum + v.success.rates[1]
    test.sum <- test.sum + v.success.rates[2]
  }
  
  if(is.null(cut) && is.null(top)) 
    v.result <- c(-1)
  else 
    v.result <- ifelse(is.null(cut), c(as.integer(top)), c(round(cut,2)))
  
  v.result <- c(v.result,
                as.integer(ncol(df.deg) - 1),
                as.integer(seed.begin),
                round(train.sum/n.times, 2),
                round(test.sum/n.times, 2),
                as.integer(n.times),
                round(total.time/n.times, 5))

  return (v.result)
}

# Runs TTRunner n-times for different num of genes (top-n or cut-k DEGs)
# Returns a df of various stats
GetSizeSuccess <- function(df, list.fc, tool, type, seed.begin, n.times, range, class.type, verbose) {
  
  cat("TOOL:", tool,"\n")
  df.result <- data.frame(top = integer(length(range)),
                          degs = integer(length(range)),
                          seed = integer(length(range)),
                          train = numeric(length(range)),
                          test = numeric(length(range)),
                          ntimes = integer(length(range)),
                          avgtime = numeric(length(range)))
  
  cat("Top\tDEGs\tSeed\tTrain\tTest\tnTimes\tAvgTime\n")
  r = 1
  for(n in range) {
    if(type == "top")
      v.result <- TTRunner(df, cut = NULL, n, list.fc, tool, 
                           seed.begin, n.times, class.type, verbose)
    else if (type == "cut")
      v.result <- TTRunner(df, n, top = NULL, list.fc, tool, 
                           seed.begin, n.times, class.type, verbose)
    
    for(c in c(1:ncol(df.result))) {
      df.result[r,c] <- v.result[c]
      if(verbose > 0)
        cat(v.result[c], "\t")
    }
    cat("\n")
    
    r <- r + 1
  }
  return (df.result)
}

# Runs TTRunner n-times for different seed.begin values
# Returns a df of various stats
GetSeedSuccess <- function(df, list.fc, n.times, tool, seeds, cut, top, class.type, verbose) {

  df.stats <- data.frame(top = integer(length(seeds)),
                         degs = integer(length(seeds)),
                         seed = integer(length(seeds)),
                         train = numeric(length(seeds)),
                         test = numeric(length(seeds)),
                         ntimes = integer(length(seeds)),
                         avgtime = numeric(length(seeds)))
  
  cat("Top\tnDEGs\tSeed\tTrain\tTest\tnTimes\tAvgTime\n")
  
  r <- 1
  for(s in seeds) {
    v.result <- TTRunner(df, cut, top, list.fc, tool, seed.begin = s, 
                         n.times, class.type, verbose)
    
    for(c in c(1:ncol(df.stats))) {
      df.stats[r,c] <- v.result[c]
      if(verbose > 0)
        cat(v.result[c], "\t")
    }
    cat("\n")
    r <- r + 1
  }
  return (df.stats)
}