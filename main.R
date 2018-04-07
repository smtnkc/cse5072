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
###################################################################

# df.result.success <- GetSuccessRates(
#   df.exp, 
#   list.fc, 
#   type = "TOP", 
#   tool = "c50", 
#   seed.val = 1000
# )

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

df.result.seed <- GetSeedSuccess(
  df.splitting.attrs, 
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 10, 
  seeds = c(1000,2000,3000,4000,5000,6000,7000,8000,9000),
  cut = NULL, 
  top = NULL, 
  verbose = 1
)

