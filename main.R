library(readr)
library(rstudioapi)
library(C50)
library(e1071) # For SVM
rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

source("preprocess.R")
source("deg_detect.R")
source("helper.R") # requires deg_detect.R

df.exp  <- as.data.frame(read_csv("data/GSE23561.csv"))
df.states  <- as.data.frame(read_csv("data/states.csv"))
df.platform  <- as.data.frame(read_csv("data/GPL10775.csv"))
df.platform <- df.platform[, c("ID", "Symbol v12")]
df.exp <- df.exp[, -1]
df.exp <- cbind("Symbol v12" = df.platform[, "Symbol v12"], df.exp)
df.exp <- FilterBadData(df.exp)
df.exp <- GroupByGeneSymbols(df.exp)
df.exp.tr <- as.data.frame(t(df.exp))
df.exp.tr <- cbind(Group = df.states[,"Disease_State"], df.exp.tr)
df.grouped <- GroupByDiseaseStates(df.exp)
list.fc <- GetListFc(df.grouped)

######## STATS ##########

df.result.size <- GetSizeSuccess(
  df = df.exp,
  list.fc = list.fc,
  tool = "c50",
  type = "top",
  seed.begin = 3000,
  n.time = 10,
  range = c(40,50,60,70,80,90,100,110,120,150,200),
  #range = c(5,10,25,50,75,100,150,250,500,1000),
  verbose = 1
)

df.deg <- GetDfDeg(df = df.exp, list.fc = list.fc, cut = NULL, top = 70, verbose = 2)
df.splitting.attrs <- GetDfSplittingAttrs(df.deg, verbose = 1)
df.splitting.attrs <- cbind(Group = df.states[,"Disease_State"], df.splitting.attrs)

df.result.seed.deg <- GetSeedSuccess(
  df.deg,
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 1, 
  seeds = c(3000:3020),
  cut = NULL, 
  top = NULL, 
  verbose = 1
)

df.result.seed.exp <- GetSeedSuccess(
  df.exp.tr,
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 1, 
  seeds = c(3000:3020),
  cut = NULL, 
  top = NULL, 
  verbose = 1
)

df.result.seed.comb <- cbind(
  df.result.seed.exp[, c(3,4,5,7)],
  df.result.seed.deg[, c(4,5,7)])
colnames(df.result.seed.comb)[2] <- "train.all"
colnames(df.result.seed.comb)[3] <- "test.all"
colnames(df.result.seed.comb)[4] <- "time.all"
colnames(df.result.seed.comb)[5] <- "train.70"
colnames(df.result.seed.comb)[6] <- "test.70"
colnames(df.result.seed.comb)[7] <- "time.70"
write.csv(df.result.seed.comb, file = "24283vs70.csv")
