# Generates a simplified decision tree model by selecting only splitting attributes
simplify_df <- function(df, verbose) {
  model <- C5.0(df[, -1], df[, 1])
  
  if(verbose) { print(model); print(summary(model)) }
  
  df_attr_imp <- C5imp(model, metric = "usage", pct = TRUE) # or metric = "split"
  imp_attr_count <- length(df_attr_imp[df_attr_imp$Overall != 0, ])
  imp_attr_list <- row.names(head(df_attr_imp, imp_attr_count, imp_attr_count))
  
  df_simple <- data.frame(Group = df[, 1])  
  df_simple <- cbind(df_simple, df[, imp_attr_list])
  return(df_simple)
}

#df_simple <- simplify_df(df_all, TRUE)
#dtm_simple <- C5.0(df_simple[, -1], df_simple[, 1])
#plot(dtm_simple)