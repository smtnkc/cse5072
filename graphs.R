# Class Histogram
ggplot(df.states, aes(x=Disease_State)) +
  geom_bar(fill = palette) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Disease State", y = "Count") 

# Data Distribution 
df.grouped.long <- gather(df.grouped)
colnames(df.grouped.long)[1] <- "Disease State"
ggplot(df.grouped.long, aes(x = value, color = `Disease State`)) + 
  ggtitle("Data Distribution of GSE23561 Dataset") +
  scale_colour_manual(values = palette) +
  geom_freqpoly(size = 1) +
  xlim(c(0,10)) +
  theme(axis.line = element_line(color = "black")) +
  labs(x = "Gene Expression", y = "Count") 

# C50 size vs execution time
ggplot(df.result.size.c50, aes(degs)) + 
  geom_smooth(aes(y = avgtime), size = 1, method = "loess") +
  labs(x = "Number of genes", y = "Execution Time (sec)") +
  xlim(c(0, 25000))








# C50 size vs accuracy
df.temp  <- as.data.frame(read_csv("../RESULTS/C50_tn100_sb5000.csv"))
ggplot(df.temp, aes(top)) +
  ggtitle("Average accuracy of C50 where seed values are in [5000:5099]") +
  geom_smooth(aes(y = test, color = "test"), size = 1, show.legend=TRUE, method = 'loess') +
  geom_smooth(aes(y = train, color = "train"), size = 1, show.legend=TRUE, method = 'loess') +
  scale_color_manual(labels=c("Test", "Train"), values=palette[1:2]) +
  labs(color="Sample", x = "N value for DegDetect", y = "Success %") +
  theme(axis.line = element_line(color = "black"))

# C50 vs SVM vs RF accuracy in different number of genes
df.c50 <- as.data.frame(read_csv("../RESULTS/C50_tn100_sb5000.csv"))
df.svm <- as.data.frame(read_csv("../RESULTS/SVM_tn100_sb5000.csv"))
df.rf <- as.data.frame(read_csv("../RESULTS/RF_tn10_sb5000.csv"))
ggplot(df.c50, aes(top)) +
  ggtitle("Average accuracy of C50 vs SVM vs RF where seed values are in [5000:5099]") +
  geom_smooth(aes(y = df.c50$test, color = "C50"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_smooth(aes(y = df.rf$test, color = "RF"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_smooth(aes(y = df.svm$test, color = "SVM"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  scale_color_manual(labels=c("C50", "RF", "SVM"), values=palette[c(1,2,3)]) +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700,800,900,1000)) +
  labs(color="TOOL", x = "N value for DegDetect", y = "Success %") +
  theme(axis.line = element_line(color = "black"))


# C50 vs SVM vs RF accuracy in different seeds (TOPN = 120, DEGs = 560)
df.c50 <- as.data.frame(read_csv("../RESULTS/C50_deg560_tn10.csv"))
df.svm <- as.data.frame(read_csv("../RESULTS/SVM_deg560_tn10.csv"))
df.rf <- as.data.frame(read_csv("../RESULTS/RF_deg560_tn10.csv"))
ggplot(df.c50, aes(seed)) +
  ggtitle("Average accuracy of C50 vs SVM vs RF (560 unique DEGs for TOP-120 cut)") +
  geom_smooth(aes(y = df.c50$test, color = "C50"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_smooth(aes(y = df.rf$test, color = "RF"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_smooth(aes(y = df.svm$test, color = "SVM"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  scale_color_manual(labels=c("C50", "RF", "SVM"), values=palette[c(1,2,3)]) +
  #scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700,800,900,1000)) +
  labs(color="TOOL", x = "Seed Begin", y = "Success %") +
  theme(axis.line = element_line(color = "black"))




# C50 vs SVM vs RF time in different seeds (TOPN = 120, DEGs = 560)
df.c50 <- as.data.frame(read_csv("../RESULTS/C50_deg560_tn10.csv"))
df.svm <- as.data.frame(read_csv("../RESULTS/SVM_deg560_tn10.csv"))
df.rf <- as.data.frame(read_csv("../RESULTS/RF_deg560_tn10.csv"))
ggplot(df.c50, aes(seed)) +
  ggtitle("Average execution time of C50 vs SVM vs RF (560 unique DEGs for TOP-120 cut)") +
  geom_line(aes(y = df.c50$avgtime, color = "C50"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_line(aes(y = df.rf$avgtime, color = "RF"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_line(aes(y = df.svm$avgtime, color = "SVM"), size = 1, show.legend=TRUE, method = 'loess', se = FALSE) +
  scale_color_manual(labels=c("C50", "RF", "SVM"), values=palette[c(1,2,3)]) +
  #scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700,800,900,1000)) +
  labs(color="TOOL", x = "Seed Begin", y = "Avg. Execution Time (sec)") +
  theme(axis.line = element_line(color = "black"))






# Only 23283
df.all <- as.data.frame(read_csv("../RESULTS/C50_all_tn10.csv"))
ggplot(df.all, aes(seed)) +
  ggtitle("Accuracy of C50 over 24283 genes (without DegDetect)") +
  geom_line(aes(y = df.all$train, color = "train.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = df.all$test, color = "test.all"), size = 1, show.legend=TRUE) +
  scale_color_manual(labels=c("test.all", "train.all"), values=palette[c(2,1)]) +
  labs(color="Sample", x = "Seed", y = "Success %") +
  theme(axis.line = element_line(color = "black"))


# 23283 vs 120
df.c50 <- as.data.frame(read_csv("../RESULTS/C50_deg560_tn10.csv"))
df.all <- as.data.frame(read_csv("../RESULTS/C50_all_tn10.csv"))
ggplot(df.c50, aes(seed)) +
  ggtitle("Average accuracy of C50 over 24283 genes vs TOP-120 cut") +
  geom_line(aes(y = df.all$train, color = "train.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = df.c50$train, color = "train.120"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = df.all$test, color = "test.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = df.c50$test, color = "test.120"), size = 1, show.legend=TRUE) +
  scale_color_manual(labels=c("test.120", "test.all", "train.120", "train.all"), 
                     values=palette[c(3,2,4,1)]) +
  labs(color="Sample", x = "Seed", y = "Success %") +
  theme(axis.line = element_line(color = "black"))


# 23283 vs 70 vs 5
df.temp  <- as.data.frame(read_csv("../RESULTS/24283vs70.csv"))
ggplot(df.temp, aes(seed)) +
  ggtitle("Accuracy of C50 over 24283 vs 70 vs 5 genes") +
  #geom_line(aes(y = train.all, color = "train.all"), size = 0.5, show.legend=TRUE) +
  #geom_line(aes(y = train.70, color = "train.70"), size = 0.5, show.legend=TRUE) +
  #geom_line(aes(y = train.split, color = "train.split"), size = 0.5, show.legend=TRUE) +
  geom_line(aes(y = test.all, color = "test.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = test.70, color = "test.70"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = test.split, color = "test.split"), size = 1, show.legend=TRUE) +
  scale_color_manual(labels=c("test.70", "test.all", "test.split"),
                              #,"train.split", "train.70", "train.all"), 
                     values=palette[c(1,2,3)]) +
  labs(color="Sample", x = "Seed", y = "Success %") +
  theme(axis.line = element_line(color = "black"))

##################################################################

ggplot(df.result.time, aes(degs)) +
  ggtitle("Average time to create a model with C50 vs SVM vs RF") +
  geom_smooth(aes(y = time.c50, color = "C50"), show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_smooth(aes(y = time.rf, color = "RF"), show.legend=TRUE, method = 'loess', se = FALSE) +
  geom_smooth(aes(y = time.svm, color = "SVM"), show.legend=TRUE, method = 'loess', se = FALSE) +
  scale_color_manual(labels=c("C50", "RF", "SVM"), values=palette[c(1,2,3)]) +
  scale_y_continuous(breaks = c(0,10,20,30,40,50)) +
  labs(color="Tool", x = "Number of Genes", y = "Execution Time (sec)") +
  theme(axis.line = element_line(color = "black"))



##################### SCATTER of SUCCESS
df.c50 <- as.data.frame(read_csv("../RESULTS/C50_deg560_tn10.csv"))
df.svm <- as.data.frame(read_csv("../RESULTS/SVM_deg560_tn10.csv"))
df.rf <- as.data.frame(read_csv("../RESULTS/RF_deg560_tn10.csv"))
df.all <- as.data.frame(read_csv("../RESULTS/C50_all_tn10.csv"))
ggplot() +
  geom_point(aes(x=mean(df.rf$train), y=mean(df.rf$test), color = "RF.560"), size=3) +
  geom_point(aes(x=mean(df.all$train), y=mean(df.all$test), color = "C50.ALL"), size=3) +
  geom_point(aes(x=mean(df.c50$train), y=mean(df.c50$test), color = "C50.560"), size=3)+
  geom_point(aes(x=mean(df.svm$train), y=mean(df.svm$test), color = "SVM.560"), size=3)+
  scale_color_manual(labels=c("C50.560", "C50.ALL", "RF.560", "SVM.560"), 
                     values=palette[c(1,4,2,3)]) +
  scale_y_continuous(breaks = c(50,55,60,65,70,75,80)) +
  labs(color="Tool", x = "Train Accuracy", y = "Test Accuracy") +
  theme(axis.line = element_line(color = "black"))







