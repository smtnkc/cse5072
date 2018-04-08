library(readr)
library(rstudioapi)
library(C50)
library(e1071) # SVM
library(randomForest)
library(ggplot2)
library(tidyr)

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))

source("preprocess.R")
source("degDetect.R")
source("helper.R")

df.exp  <- as.data.frame(read_csv("data/GSE23561.csv"))
df.states  <- as.data.frame(read_csv("data/states.csv"))
df.platform  <- as.data.frame(read_csv("data/GPL10775.csv"))
df.platform <- df.platform[, c("ID", "Symbol v12")]
df.exp <- df.exp[, -1]
df.exp <- cbind("Symbol v12" = df.platform[, "Symbol v12"], df.exp)
df.exp <- FilterBadData(df.exp)
df.exp[, "Symbol v12"] <- make.names(df.exp[, "Symbol v12"])
df.exp <- GroupByGeneSymbols(df.exp)
df.exp.tr <- as.data.frame(t(df.exp))
df.exp.tr <- cbind(Group = df.states[,"Disease_State"], df.exp.tr)
df.grouped <- GroupByDiseaseStates(df.exp)
list.fc <- GetListFc(df.grouped)
palette <- c("#12a8e5", "#23366a", "#f44a26", "#e4c19f", "#747a74")