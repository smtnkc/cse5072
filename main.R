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

df.stats <- GetSuccessRates(df.exp, list.fc, type = "TOP", tool = "c50", seed.val = 1000)

df.deg <- GetDfDeg(df = df.exp, list.fc = list.fc, cut = NULL, top = 20, verbose = 2)
# df.splitting.attrs <- GetDfSplittingAttrs(df.deg, verbose = 1)
# df.splitting.attrs <- cbind(Group = df.states[,"Disease_State"], df.splitting.attrs)

TTRunner(#df = df.exp.tr,
         df = df.deg,
          tool = "c50", 
          seed.begin = 2000,
          n.times = 100,
          verbose = 2)

# TT(df.deg, "c50", 6, 3)
