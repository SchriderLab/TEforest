# script to plot features from data to understand their distributions
library(tidyverse)
library(ggpubr)

data <- read.csv("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/feature_data.csv")
data <- read.csv("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/feature_data_shrink.csv")
data <- data[, !names(data) %in% c("X", "insert_size")]
data_summary <- data %>% group_by(true,read_length) %>% summarize(across(everything(), list(mean = mean, sd = sd)))

data$true <- data$true/2

more_data <- read.csv("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/feature_data_new.csv")
more_data <- more_data[, !names(more_data) %in% c("X", "insert_size")]
more_data_summary <- more_data %>% group_by(true,read_length) %>% summarize(across(everything(), list(mean = mean, sd = sd)))
more_data$true <- more_data$true/2

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
data$read_length <- as.factor(data$read_length)

more_data$true <- as.factor(more_data$true)
more_data$read_length <- as.factor(more_data$read_length)

ggplot(data = data, aes(x = TE_specific_Orphan_Read_mean)) +
  geom_boxplot(aes(fill = true), 
                 alpha = 0.5, 
                 position = "identity") +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "Frequency", title="TE_specific_Orphan_Read_mean by class")+
  #xlim(0, 25) + ylim(0,500) +
  facet_wrap("read_length", ncol=1)

ggplot(data = data, aes(x = TE_specific_Orphan_Read_sd)) +
  geom_histogram(aes(fill = true), 
                 alpha = 0.5, 
                 position = "identity") +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_sd",
       y = "Frequency", title="TE_specific_Orphan_Read_sd by class")+
  xlim(0, 25) + ylim(0,500) +
  facet_wrap("read_length", ncol=1)


ggplot(data = more_data, aes(x = TE_specific_Orphan_Read_IQR)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = T) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "zygosity", title="TE_specific_Orphan_Read_mean by class") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

ggplot(data = more_data, aes(x = TE_specific_Orphan_Read_IQR, y = true)) +
  geom_boxplot(aes(color = true), width = 0.7, show.legend = T) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "zygosity", title="TE_specific_Orphan_Read_IQR") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

ggplot(data = more_data, aes(x = Proper_Pair_IQR, y = true)) +
  geom_boxplot(aes(color = true), width = 0.3, show.legend = T) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Orphan_Read_mean",
       y = "zygosity", title="Proper_Pair_IQR") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

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


### PCA 
pca_data_subset <- more_data %>% ungroup() %>% subset(true!=0)
pca_data <- pca_data_subset %>% select_if(~var(.) != 0) 

# Select only the columns to be used for PCA (excluding `true` and `read_length`)
pca_data <- pca_data %>% 
  select(-true, -read_length) 

# Perform PCA
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# Create a data frame with the PCA results
pca_df <- as_tibble(pca_result$x)

# Add the grouping variables back to the PCA result
pca_df <- bind_cols(pca_data_subset %>% select(true, read_length), pca_df)

# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(true), shape = as.factor(true))) +
  geom_point(size = 3) +
  labs(title = "PCA of Grouped Data",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Read Length",
       shape = "Zygosity") +
  theme_minimal()

# Get the loadings (contributions of each variable to the principal components)
loadings <- as_tibble(pca_result$rotation, rownames = "variable")

# Identify the top contributing variables for PC1 and PC2
top_pc1 <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n = 10) %>%
  pull(variable)

top_pc2 <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n = 10) %>%
  pull(variable)

# Print the top contributing variables
cat("Top variables contributing to PC1:\n", top_pc1, "\n")
cat("Top variables contributing to PC2:\n", top_pc2, "\n")




# Define the initial feature names
feature_names <- c(
  "Cigar1",
  "Cigar2",
  "Cigar3",
  "Cigar4",
  "Cigar5",
  "Paired",
  "Proper_Pair",
  "Is_Read1_Unmapped",
  "Is_Read2_Unmapped",
  "Is_Read1_Rev_Comp",
  "Is_Read2_Rev_Comp",
  "Is_First_Read",
  "Is_Second_Read",
  "Split",
  "Long_Insert",
  "Short_Insert",
  "Parallel_Read",
  "Everted_Read",
  "Orphan_Read",
  "Insert_Size",
  "Quality"
)

# Create the extended feature list
feature_list <- c()
for (x in feature_names) {
  feature_list <- c(feature_list, paste0(x, "_mean"), paste0(x, "_median"), paste0(x, "_sd"), paste0(x, "_IQR"))
}

# Extend the feature list with TE specific features
feature_list_extended <- c(feature_list, paste0("TE_specific_", feature_list))
feature_list_extended_df <- feature_list_extended
feature_list_extended_df

data <- read.table("/nas/longleaf/home/adaigle/work/test_TEforest/test_basenorm_feats/3L3RX/feature_importances_all.txt", header=T)
data$Feature <- feature_list_extended[data$Feature + 1]  # +1 to account for 0-based indexing
data <- data %>% filter(Importance!=0)
data <- data[1:10,]
plot <- ggplot(data, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = Importance - Std, ymax = Importance + Std), width = 0.2) +
  coord_flip() +
  labs(title = "", x = "Feature", y = "Importance") +
  theme_minimal() +
  theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), axis.text.y = element_text(size=16), axis.text.x = element_text(size=16))

output_path <- "/nas/longleaf/home/adaigle/TEforest/plots/feat_importance_30X.png"
ggsave(output_path, plot = plot, width = 10, height = 5)
