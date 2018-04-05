GetDfFc <- function(df, target.index) {
  
  df.fc <- data.frame(fc.min = numeric(nrow(df)))
  rownames(df.fc) <- rownames(df)
  
  if(target.index == 1)
    range <- 2:ncol(df)
  else if(target.index == ncol(df)) 
    range <- 1:(ncol(df)-1)
  else
    range <- c(1:(target.index - 1),(target.index + 1):ncol(df))
  
  for(row in 1:nrow(df)) {
    fc.min <- 9999
    for(pair.index in range) {
      fc <- abs(log(df[row, target.index]/df[row, pair.index], base = 2))
      if(fc < fc.min) fc.min <- fc
    }
    df.fc[row, 1] <- round(fc.min,2)
  }
  return(df.fc)
}

GetListFc <- function(df) {
  list.fc <- list()
  for(i in 1:ncol(df))
    list.fc[[i]] <- GetDfFc(df, i)
  return (list.fc)
}

GetDfByTop <- function(df, top) {
  df.top <- as.data.frame(head(df[order(-df[,1]), , drop = FALSE], top))
  return (df.top)
}

GetDfByCut <- function(df, cut) {
  df.cut <- as.data.frame(df[df[, 1] >= cut, , drop = FALSE])
  return (df.cut)
}

GetDegNames <- function(df, cut, top) {

  if(is.null(cut) && is.null(top))
    stop("You must set either one of k or n value!")

  if(!is.null(cut))
    v.deg.names <- rownames(GetDfByCut(df, cut))
  else
    v.deg.names <- rownames(GetDfByTop(df, top))

  return(v.deg.names)
}

GetDfDeg <- function(df, list.fc, cut, top, verbose) {
  v.combined.deg.names <- c()
  
  for(i in 1:length(list.fc)) {
    
    v.deg.names <- GetDegNames(list.fc[[i]], cut, top)
    v.combined.deg.names <- c(v.combined.deg.names, v.deg.names)
    
    if(verbose > 1) 
      cat("DEG count for state", i, "is", length(v.deg.names), "\n")
  }
  
  v.unique.deg.names <- unique(v.combined.deg.names)
  
  if(verbose > 1) {
    cat("Total DEG count is", length(v.combined.deg.names), "\n")
    cat("Unique DEG count is", length(v.unique.deg.names), "\n")
  }

  df.deg <- as.data.frame(t(df[v.unique.deg.names,]))
  df.deg <- cbind(Group = df.states[,"Disease_State"], df.deg)
  
  return (df.deg)
}