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