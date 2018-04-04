df_exp  <- as.data.frame(read_csv("data/GSE23561.csv", col_names = TRUE)) # gene expressions
df_states  <- as.data.frame(read_csv("data/states.csv", col_names = TRUE)) # disease states
df_platform  <- as.data.frame(read_csv("data/GPL10775.csv", col_names = TRUE)) # platform

df_platform <- df_platform[, c("ID", "Symbol v12")] # remove unnecessary columns
df_exp <- df_exp[, -1] # remove id_ref column
df_exp <- cbind("Symbol v12" = df_platform[, "Symbol v12"], df_exp) # add the symbols column
#df_exp <- cbind("Symbol v12" = sprintf("`%s`", df_platform[, "Symbol v12"]), df_exp)

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

# GROUP COLUMNS BY DISEASE STATES:
df_group <- NULL
df_group <- data.frame(Control = rowMeans(df_exp[ ,1:9]))                 # 9 Control
df_group <- cbind(df_group, data.frame(RA = rowMeans(df_exp[ ,10:15])))   # 6 RA
df_group <- cbind(df_group, data.frame(MetS = rowMeans(df_exp[ ,16:21]))) # 6 MetS
df_group <- cbind(df_group, data.frame(CAD = rowMeans(df_exp[ ,22:27])))  # 6 CAD
df_group <- cbind(df_group, data.frame(T2D = rowMeans(df_exp[ ,28:35])))  # 8 T2D
df_group <- round(df_group, 2)
row.names(df_group) <- symbols

df_all <- as.data.frame(t(df_exp))
df_all <- cbind(Group = df_states[,"Disease_State"], df_all)