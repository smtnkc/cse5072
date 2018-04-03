# If there is any attribute which is NOT differentially
# expressed as against the target attribute, delete that row.
# Set related params to select a cutoff or top N DEGs.
deg_detect <- function(df, target_index, cutoff, topN) {
  
  if(is.null(cutoff) && is.null(topN)) {
    stop("Set either one of cutoff or topN!")
  }
  else {
    if(!is.null(cutoff) && cutoff < 0)
      stop("cutoff must be a positive integer!")
    else if(!is.null(topN) && topN < 0)
      stop("topN must be a positive integer!")
  }
  
  df_minfc <- data.frame(MIN = numeric(nrow(df)))
  #print(dim(df_minfc))
  row.names(df_minfc) <- row.names(df)
  
  if(target_index == 1)
    range <- 2:ncol(df)
  else if(target_index == ncol(df)) 
    range <- 1:(ncol(df)-1)
  else
    range <- c(1:(target_index-1),(target_index+1):ncol(df))
  
  for(row in 1:nrow(df)) {
    minfc <- 999999
    
    for(pair_index in range) {
      fc <- abs(log(df[row, target_index]/df[row, pair_index], base = 2))
      #df_minfc[row, colnames(df)[pair_index]] <- round(fc,2)
      
      if(fc <= minfc) {
        minfc = fc
      }
    }
    
    if(!is.null(cutoff)) {
      if(minfc >= cutoff)
        df_minfc[row, "MIN"] <- round(minfc,2)
      else 
        df_minfc[row, "MIN"] <- NA
    }
    else {
      df_minfc[row, "MIN"] <- round(minfc,2)
    }
  }
  
  if(!is.null(cutoff)) { 
    df_minfc <- na.omit(df_minfc)
  }
  else { df_minfc <- as.data.frame(head(df_minfc[order(-df_minfc$MIN), , drop = FALSE], topN)) }
  
  cat("DEG count for column", target_index, "=", nrow(df_minfc), "\n")
  
  return(df_minfc)
}
get_deg_list <- function(df, cutoff, topN) {
  start_time <- Sys.time()
  df_deg_all <- NULL
  for(i in 1:5) {
    if(i == 1) df_deg_all <- deg_detect(df, target_index = i, cutoff = cutoff, topN = topN)
    else df_deg_all <- rbind(df_deg_all, deg_detect(df, target_index = i, cutoff = cutoff, topN = topN))
  }
  end_time <- Sys.time()
  cat("Getting DEGS takes", end_time-start_time, "sec\n")
  return (row.names(df_deg_all))
}

deg_list <- get_deg_list(df_group, cutoff = DEG_CUTOFF, topN = DEG_TOP_N)

df_deg <- as.data.frame(t(df_exp[deg_list,]))
df_deg <- cbind(Group = df_states[,"Disease_State"], df_deg)
#dtm4deg <- C5.0(df_deg[,-1], df_deg[,1])
#dtm4deg; summary(dtm4deg)
