library(ggplot2)
library(tidyr)

palette <- c("#12a8e5", "#23366a", "#f44a26", "#e4c19f", "#747A74")

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
