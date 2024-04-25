# script to plot features from data to understand their distributions
library(tidyverse)
library(ggpubr)

data <- read.csv("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/feature_data.csv")
data <- data[, !names(data) %in% c("X", "read_length", "insert_size")]
data_summary <- data %>% group_by(true) %>% summarize(across(everything(), list(mean = mean, sd = sd)))



more_data <- read.csv("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/feature_data_old.csv")
more_data <- more_data[, !names(more_data) %in% c("X", "read_length", "insert_size")]
more_data_summary <- more_data %>% group_by(true) %>% summarize(across(everything(), list(mean = mean, sd = sd)))

summary <- rbind(data_summary, more_data_summary) %>% arrange(true)

# Assuming df1 and df2 are your two dataframes with identical data columns

# Create a list of variable names to compare
variables <- names(data)

# Loop through each variable and compare the distributions
for (variable in variables) {
  compare_means(as.data.frame(data[[variable]]), as.data.frame(more_data[[variable]]), method = "t.test", 
                method.args = list(var.equal = TRUE), 
                group.labels = c("Dataset 1", "Dataset 2"),
                title = paste("Comparison of", variable))
}

data$true <- as.factor(data$true)
more_data$true <- as.factor(more_data$true)

ggplot(data = data, aes(x = TE_specific_Orphan_Read_sd)) +
  geom_histogram(aes(fill = true), 
                 alpha = 0.5, 
                 position = "identity") +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_sd",
       y = "Frequency", title="TE_specific_Orphan_Read_sd by class")+
  xlim(0, 25) + ylim(0,1500)

ggplot(data = data, aes(x = TE_specific_Orphan_Read_mean)) +
  geom_histogram(aes(fill = true), 
                 alpha = 0.5, 
                 position = "identity") +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Orphan_Read_mean by class")+
  xlim(0, 25) + ylim(0,1500)

ggplot(data = more_data, aes(x = TE_specific_Orphan_Read_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Orphan_Read_mean by class")

ggplot(data = data, aes(x = TE_specific_Insert_Size_mean)) +
  geom_histogram(aes(fill = true), 
                 alpha = 0.5, 
                 position = "identity") +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Insert_Size_mean by class")+
  xlim(0, 2000) #+ ylim(0,1500)

ggplot(data = data, aes(x = TE_specific_Insert_Size_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Insert_Size_mean",
       y = "Frequency", title="TE_specific_Insert_Size_mean by class")

plot(data$TE_specific_Orphan_Read_mean, as.numeric(data$true)-1)

head(data)
data$orphan_read_fraction <- data$TE_specific_Orphan_Read_mean / data$Orphan_Read_mean
data$Cigar1_frac <- data$TE_specific_Paired_mean / data$Paired_mean

ggplot(data = data, aes(x = Cigar1_frac)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "orphan_read_fraction",
       y = "Frequency", title="orphan_read_fraction by class")


ggplot(data = data, aes(x = TE_specific_Orphan_Read_mean / Proper_Pair_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "orphan_read_fraction",
       y = "Frequency", title="orphan_read_fraction by class")+
  xlim(0, 2.5)

ggplot(data = data, aes(x = TE_specific_Paired_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "orphan_read_fraction",
       y = "Frequency", title="orphan_read_fraction by class")
data$orphan_read_fraction <- NULL



ggplot(data = more_data, aes(x = TE_specific_Quality_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Quality_mean by class")

ggplot(data = data, aes(x = TE_specific_Quality_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Quality_mean by class")

ggplot(data = more_data, aes(x = TE_specific_Parallel_Read_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Parallel_Read_mean by class")

ggplot(data = data, aes(x = TE_specific_Parallel_Read_mean)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Parallel_Read_mean by class")
