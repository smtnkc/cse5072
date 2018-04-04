library(readr)
library(rstudioapi)
library(C50)
library(e1071) # For SVM
rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

source("preprocess.R")
source("deg_detect.R")
source("helper.R")

SuccessFinder(df.exp, list.fc, type = "TOP", tool = "c50")
df.deg <- GetDfDeg(df = df.exp, list.fc = list.fc, cut = NULL, top = 100, verbose = 2)
df.splitting.attrs <- GetDfSplittingAttrs(df.deg, verbose = 1)
df.splitting.attrs <- cbind(Group = df.states[,"Disease_State"], df.splitting.attrs)

TTRunner(df = df.splitting.attrs,
          train.perc = 0.8, 
          tool = "c50", 
          seed.begin = 5000,
          n.times = 50,
          verbose = 2)
