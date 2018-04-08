
############## ALL SEED SUCCESS ###########################
df.seed.all.c50 <- GetSeedSuccess(
  df.exp.tr,
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 1, 
  seeds = c(1000,2000,3000,4000,5000),
  cut = NULL, 
  top = NULL,
  class.type = "char",
  verbose = 1)

############## DEG SEED SUCCESS ##########################
df.deg <- GetDfDeg(df = df.exp, list.fc = list.fc, cut = NULL, top = 120, class.type = "char", verbose = 2)

df.seed.all.c50 <- GetSeedSuccess(
  df.exp.tr,
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 10, 
  seeds = c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
  cut = NULL, 
  top = NULL,
  class.type = "char",
  verbose = 1
)

df.seed.deg.c50 <- GetSeedSuccess(
  df.deg,
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 50, 
  seeds = c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
  cut = NULL, 
  top = NULL,
  class.type = "char",
  verbose = 1
)

df.seed.deg.svm <- GetSeedSuccess(
  df.deg,
  tool = "svm", 
  list.fc = list.fc, 
  n.times = 50, 
  seeds = c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
  cut = NULL, 
  top = NULL,
  class.type = "char",
  verbose = 1
)

df.seed.deg.rf <- GetSeedSuccess(
  df.deg,
  tool = "rf", 
  list.fc = list.fc, 
  n.times = 50, 
  seeds = c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
  cut = NULL, 
  top = NULL,
  class.type = "char",
  verbose = 1
)
#write.csv(df.seed.all.c50, file = "../RESULTS/C50_all_tn10.csv", row.names = FALSE)
#write.csv(df.seed.deg.c50, file = "../RESULTS/C50_deg560_tn10.csv", row.names = FALSE)
#write.csv(df.seed.deg.svm, file = "../RESULTS/SVM_deg560_tn10.csv", row.names = FALSE)
#write.csv(df.seed.deg.rf, file = "../RESULTS/RF_deg560_tn10.csv", row.names = FALSE)






############## SPLITTING SEED SUCCESS ###########################
df.splitting.attrs <- GetDfSplittingAttrs(df.deg, class.type = "char", verbose = 1)
 df.splitting.attrs <- GetDfSplittingAttrs(df.exp.tr, class.type = "num", verbose = 1)
 plot(C5.0(df.splitting.attrs[, -1], df.splitting.attrs[, 1]))
df.result.seed.split <- GetSeedSuccess(
  df.splitting.attrs,
  tool = "c50", 
  list.fc = list.fc, 
  n.times = 1, 
  seeds = c(3000:3020),
  cut = NULL, 
  top = NULL,
  class.type = "char",
  verbose = 1)
df.result.seed.comb  <- as.data.frame(read_csv("../RESULTS/24283vs70.csv"))
df.result.seed.comb <- cbind(
  df.result.seed.comb,
  df.result.seed.split[, c(4,5,7)])
colnames(df.result.seed.comb)[8] <- "train.split"
colnames(df.result.seed.comb)[9] <- "test.split"
colnames(df.result.seed.comb)[10] <- "time.split"
write.csv(df.result.seed.comb, file = "../RESULTS/24283vs70.csv", row.names = FALSE)




############## COMPARE C50-SVM-RF ###############

df.result.size.c50 <- GetSizeSuccess(
  df = df.exp,
  list.fc = list.fc,
  tool = "c50",
  type = "top",
  seed.begin = 5000,
  n.time = 100,
  range = c(10,20,30,40,50,60,70,80,90,100,110,120,150,200,300,400,500,1000),
  class.type = "char",
  verbose = 1)

# write.csv(df.result.size.c50, file = "../RESULTS/C50_tn100_sb5000.csv", row.names = FALSE)

df.result.size.svm <- GetSizeSuccess(
  df = df.exp,
  list.fc = list.fc,
  tool = "svm",
  type = "top",
  seed.begin = 5000,
  n.time = 100,
  range = c(10,20,30,40,50,60,70,80,90,100,110,120,150,200,300,400,500,1000),
  class.type = "char",
  verbose = 1)

#write.csv(df.result.size.svm, file = "../RESULTS/SVM_tn100_sb5000.csv", row.names = FALSE)

df.result.size.rf <- GetSizeSuccess(
  df = df.exp,
  list.fc = list.fc,
  tool = "rf",
  type = "top",
  seed.begin = 3000,
  n.time = 10,
  range = c(10,20,30,40,50,60,70,80,90,100,110,120,150,200,300,400,500,1000),
  class.type = "char",
  verbose = 1)

# write.csv(df.result.size.rf, file = "../RESULTS/RF_tn10_sb5000.csv", row.names = FALSE)

df.result.time <- as.data.frame(cbind(
  top = df.result.size.c50[,"top"], 
  degs = df.result.size.c50[,"degs"],
  time.c50 = df.result.size.c50[,"avgtime"], 
  time.svm = df.result.size.svm[,"avgtime"], 
  time.rf = df.result.size.rf[,"avgtime"]))

write.csv(df.result.time, file = "../RESULTS/C50_RF_SVM.csv", row.names = FALSE)
