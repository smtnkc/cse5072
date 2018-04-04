get_fc_table <- function(df, target_index) {
  
  df_fc <- data.frame(minfc = numeric(nrow(df)))
  rownames(df_fc) <- rownames(df)
  
  if(target_index == 1)
    range <- 2:ncol(df)
  else if(target_index == ncol(df)) 
    range <- 1:(ncol(df)-1)
  else
    range <- c(1:(target_index-1),(target_index+1):ncol(df))
  
  for(row in 1:nrow(df)) {
    minfc <- 9999
    for(pair_index in range) {
      fc <- abs(log(df[row, target_index]/df[row, pair_index], base = 2))
      if(fc < minfc) minfc <- fc
    }
    df_fc[row, 1] <- round(minfc,2)
  }
  return(df_fc)
}

generate_fc_list <- function(df) {
  fc_list <- list()
  for(i in 1:ncol(df))
    fc_list[[i]] <- get_fc_table(df, i)
  return (fc_list)
}

fc_list <- generate_fc_list(df_group) # Generates FC values for only once

get_df_top <- function(df, top_n) {
  df_top_n <- as.data.frame(head(df[order(-df[,1]), , drop = FALSE], top_n))
  return (df_top_n)
}

get_df_cut <- function(df, cut_k) {
  df_cut_k <- as.data.frame(df[df[, 1] >= cut_k, , drop = FALSE])
  return (df_cut_k)
}

get_deg_list <- function(df, cut_k, top_n) {

  if(is.null(cut_k) && is.null(top_n))
    stop("You must set either one of k or n value!")

  if(!is.null(cut_k))
    res <- rownames(get_df_cut(df, cut_k))
  else
    res <- rownames(get_df_top(df, top_n))

  return(res)
}

get_df_deg <- function(df, fc_list, cut_k, top_n) {
  combined_deg_list <- c()
  for(i in 1:length(fc_list)) {
    deg_list <- get_deg_list(fc_list[[i]], cut_k, top_n)
    combined_deg_list <- c(combined_deg_list, deg_list)
    cat("DEG count for state", i, "is", length(deg_list), "\n")
  }
  
  unique_deg_list <- unique(combined_deg_list)
  cat("Total DEG count is", length(combined_deg_list), "\n")
  cat("Unique DEG count is", length(unique_deg_list), "\n")
  df_deg <- as.data.frame(t(df[unique_deg_list,]))
  df_deg <- cbind(Group = df_states[,"Disease_State"], df_deg)
  return (df_deg)
}