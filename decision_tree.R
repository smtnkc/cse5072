library(readr)
library(rstudioapi)

rm(list=ls())

setwd(dirname(getSourceEditorContext()$path)) # Set working directory

df_exp  <- as.data.frame(read_csv("data/GSE23561.csv", col_names = TRUE)) # gene expressions
df_states  <- as.data.frame(read_csv("data/states.csv", col_names = TRUE)) # disease states
df_platform  <- as.data.frame(read_csv("data/GPL10775.csv", col_names = TRUE)) # platform

df_platform <- df_platform[, c("ID", "Symbol v12")] # remove unnecessary columns
df_exp <- df_exp[, -1] # remove id_ref column
df_exp <- cbind("Symbol v12" = df_platform[, "Symbol v12"], df_exp) # add the symbols column

# CONVERT EACH INVALID VALUE TO NA:
for(row in 1:nrow(df_exp)) {
  val <- df_exp[row, "Symbol v12"]
  if(!is.null(val) && !is.na(val))
    if(val == "." || val == "")
      df_exp[row, "Symbol v12"] <- NA
}

# REMOVE NA VALUES:
df_exp <- na.omit(df_exp)

# GROUP ROWS BY GENE SYMBOLS:
df_exp <- aggregate(df_exp[, -1], list(df_exp[, 1]), mean)
symbols <- unlist(df_exp[, 1])
df_exp <- round(df_exp[, -1], 2)
row.names(df_exp) <- symbols

################### DECISION TREE ###########################

library(C50)

simplify4plot <- function(df) {
  tree <- C5.0(df[, -1], df[, 1])
  tree
  summary(tree)

  df_attr_imp <- C5imp(tree, metric = "usage", pct = TRUE)
  imp_attr_count <- length(df_attr_imp[df_attr_imp$Overall != 0, ])
  imp_attr_list <- row.names(head(df_attr_imp, imp_attr_count, imp_attr_count))
  
  df_pruned <- data.frame(Group = df[, 1])  
  df_pruned <- cbind(df_pruned, df[, imp_attr_list])
  return(df_pruned)
}

df_trans <- as.data.frame(t(df_exp))
colnames(df_trans) <- sprintf("`%s`", symbols)
df_trans <- cbind(Group = df_states[,"Disease_State"], df_trans)

df_trans_simple <- simplify4plot(df_trans)
tree_simple <- C5.0(df_trans_simple[, -1], df_trans_simple[, 1])
plot(tree_simple)





