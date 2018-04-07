library(ggplot2)
library(tidyr)

palette <- c("#12a8e5", "#23366a", "#f44a26", "#e4c19f", "#747A74")

df.result.seed.comb  <- as.data.frame(read_csv("../RESULTS/24283vs70.csv"))

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
  scale_colour_manual(values = palette) +
  geom_freqpoly() +
  xlim(c(0,10)) +
  theme(axis.line = element_line(color = "black")) +
  labs(x = "Gene Expression", y = "Count") 

# Stats
ggplot(df.stats, aes(n)) + 
  geom_smooth(aes(y = time), method = "loess") +
  labs(x = "Number of genes", y = "Execution Time (sec)") +
  xlim(c(0, 25000))


# Results
ggplot(df.result.size, aes(top)) +
  ggtitle("Average results for 100 different sampling where seed = [5000:5099]") +
  geom_smooth(aes(y = test, color = "test"), show.legend=TRUE, method = 'loess') +
  geom_smooth(aes(y = train, color = "train"), show.legend=TRUE, method = 'loess') +
  scale_color_manual(labels=c("Test", "Train"), values=palette[1:2]) +
  labs(color="Sample", x = "Number of genes", y = "Success %") +
  theme(axis.line = element_line(color = "black"))



# 23283 vs 70
ggplot(df.result.seed.comb, aes(seed)) +
  ggtitle("Accuracy over 24283 genes vs top 70 selected") +
  geom_line(aes(y = train.all, color = "train.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = train.70, color = "train.70"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = test.all, color = "test.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = test.70, color = "test.70"), size = 1, show.legend=TRUE) +
  scale_color_manual(labels=c("test.70", "test.all", "train.70", "train.all"), 
                     values=palette[c(3,2,4,1)]) +
  labs(color="Sample", x = "Seed", y = "Success %") +
  theme(axis.line = element_line(color = "black"))

# Only 23283
ggplot(df.result.seed.comb, aes(seed)) +
  ggtitle("Accuracy over 24283 genes (without DegDetect)") +
  geom_line(aes(y = train.all, color = "train.all"), size = 1, show.legend=TRUE) +
  geom_line(aes(y = test.all, color = "test.all"), size = 1, show.legend=TRUE) +
  scale_color_manual(labels=c("test.all", "train.all"), values=palette[c(2,1)]) +
  labs(color="Sample", x = "Seed", y = "Success %") +
  theme(axis.line = element_line(color = "black"))


# 23283 vs 70 vs 5
ggplot(df.result.seed.comb, aes(seed)) +
  #ggtitle("Accuracy over 24283 vs 70 vs 5 genes") +
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




