df.exp  <- as.data.frame(read_csv("data/GSE23561.csv", col_names = TRUE)) # gene expressions
df.states  <- as.data.frame(read_csv("data/states.csv", col_names = TRUE)) # disease states
df.platform  <- as.data.frame(read_csv("data/GPL10775.csv", col_names = TRUE)) # platform

df.platform <- df.platform[, c("ID", "Symbol v12")] # remove unnecessary columns
df.exp <- df.exp[, -1] # remove id_ref column
df.exp <- cbind("Symbol v12" = df.platform[, "Symbol v12"], df.exp) # add the symbols column
#df.exp <- cbind("Symbol v12" = sprintf("`%s`", df.platform[, "Symbol v12"]), df.exp)

FilterBadData <- function(df) {
  for(row in 1:nrow(df)) {
    val <- df[row, "Symbol v12"]
    if(!is.null(val) && !is.na(val))
      if(val == "." || val == "")
        df[row, "Symbol v12"] <- NA
  }
  df <- na.omit(df)
  return (df)
}

GroupByGeneSymbols <- function(df) {
  df <- aggregate(df[, -1], list(df[, 1]), mean)
  v.symbols <- unlist(df[, 1])
  df <- round(df[, -1], 2)
  row.names(df) <- v.symbols
  return(df)
}

GroupByDiseaseStates <- function(df) {
  df.group <- NULL
  df.group <- data.frame(Control = rowMeans(df.exp[ ,1:9]))             # 9 Control
  df.group <- cbind(df.group, data.frame(RA = rowMeans(df[ ,10:15])))   # 6 RA
  df.group <- cbind(df.group, data.frame(MetS = rowMeans(df[ ,16:21]))) # 6 MetS
  df.group <- cbind(df.group, data.frame(CAD = rowMeans(df[ ,22:27])))  # 6 CAD
  df.group <- cbind(df.group, data.frame(T2D = rowMeans(df[ ,28:35])))  # 8 T2D
  df.group <- round(df.group, 2)
  row.names(df.group) <- row.names(df)
  return(df.group)
}

df.exp <- FilterBadData(df.exp)
df.exp <- GroupByGeneSymbols(df.exp)
df.grouped <- GroupByDiseaseStates(df.exp)

df.all <- as.data.frame(t(df.exp))
df.all <- cbind(Group = df.states[,"Disease_State"], df.all)