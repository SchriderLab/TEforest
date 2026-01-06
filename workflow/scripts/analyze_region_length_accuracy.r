# Load required libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)

# Specify your file paths (adjust these paths as needed)
files <- c("/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/3L3RX_classifer/A2_A3.csv",
           "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/3L3RX_classifer/AKA-017_GIM-024.csv",
           "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/3L3RX_classifer/JUT-008_MUN-009.csv")

# Read and combine CSV files, adding a 'sample' column from the file name (without extension)
data_list <- lapply(files, function(f) {
  df <- read_csv(f, col_types = cols())
  # Extract sample name from the CSV file name (e.g., "A2_A3" from "A2_A3.csv")
  sample_name <- tools::file_path_sans_ext(basename(f))
  df <- df %>% mutate(sample = sample_name)
  return(df)
})
data <- bind_rows(data_list)


# The 'file' column contains strings like:
# "A2_A3-3L-10089483-10090059-2-Stalker2"
# Split by "-" and extract: sample, chrom, start, stop (assumed to be positions 1,2,3,4 respectively)
data <- data %>%
  mutate(parts = str_split(file, "-", simplify = TRUE),
         sample_extracted = parts[,1],  # may duplicate 'sample'
         chrom  = parts[,2],
         start  = as.numeric(parts[,3]),
         stop   = as.numeric(parts[,4]),
         region_length = stop - start)

# Classify each row: 
# - True Positive if (true, pred) are (2,2) or (1,1)
# - True Negative if (true, pred) are (0,0)
# - Otherwise, mark as Misclassified
data <- data %>%
  mutate(classification = case_when(
    (true == 2 & pred == 2) | (true == 1 & pred == 1) ~ "True Positive",
    (true == 0 & pred == 0) ~ "True Negative",
    TRUE ~ "Misclassified"
  ))

# First Plot: Distribution of region lengths by the true value (0, 1, 2)
p1 <- ggplot(data, aes(x = region_length, fill = factor(true), group = factor(true))) +
  geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))), bins = 50, position = "identity", alpha = 0.7) +
  labs(title = "Distribution of Region Lengths by True Value",
       x = "Region Length",
       y = "Proportion",
       fill = "True Value") +
  theme_minimal()
print(p1)

# Create an accuracy indicator: correct if classification is either True Positive or True Negative
data <- data %>%
  mutate(correct = ifelse(classification == "Misclassified", 0, 1))

# Option 1: Bin region lengths and compute mean accuracy per bin
# Here we use deciles as bins
data <- data %>%
  mutate(length_bin = cut(region_length,
                          breaks = quantile(region_length, probs = seq(0, 1, 0.1), na.rm = TRUE),
                          include.lowest = TRUE))

accuracy_by_bin <- data %>%
  group_by(length_bin) %>%
  summarise(mean_accuracy = mean(correct),
            mean_length = mean(region_length),
            count = n())

p2 <- ggplot(accuracy_by_bin, aes(x = mean_length, y = mean_accuracy)) +
  geom_point() +
  geom_line() +
  labs(title = "Accuracy vs. Region Length (Binned)",
       x = "Mean Region Length (per bin)",
       y = "Accuracy") +
  theme_minimal()
print(p2)

# Option 2: Plot accuracy as a continuous variable using a loess smoother
p3 <- ggplot(data, aes(x = region_length, y = correct)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess") +
  labs(title = "Accuracy vs. Region Length (Loess Smoothing)",
       x = "Region Length",
       y = "Accuracy") +
  theme_minimal()
print(p3)

# Filter data for rows with true value 1 or 2
data_filtered <- data %>% 
  filter(true %in% c(1, 2))

# Option 1: Binned Plot
# Create bins for region_length using quantiles (ignoring NA values)
data_filtered <- data_filtered %>%
  mutate(length_bin = cut(region_length, 
                          breaks = quantile(region_length, probs = seq(0, 1, 0.1), na.rm = TRUE), 
                          include.lowest = TRUE))

# Compute average accuracy (where correct==1 for True Positive) in each bin
accuracy_by_bin_filtered <- data_filtered %>%
  group_by(length_bin) %>%
  summarise(mean_accuracy = mean(correct),
            mean_length = mean(region_length),
            count = n())

p4 <- ggplot(accuracy_by_bin_filtered, aes(x = mean_length, y = mean_accuracy)) +
  geom_point() +
  geom_line() +
  labs(title = "Accuracy vs. Region Length (Binned) for True=1 or 2",
       x = "Mean Region Length (per bin)",
       y = "Accuracy") +
  theme_minimal()
print(p4)

# Option 2: Continuous Plot using Loess Smoothing
p5 <- ggplot(data_filtered, aes(x = region_length, y = correct)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess") +
  labs(title = "Accuracy vs. Region Length (Loess Smoothing) for True=1 or 2",
       x = "Region Length",
       y = "Accuracy") +
  theme_minimal()
print(p5)

# Filter for misclassified entries and split into two types
misclassified_data <- data %>% 
  filter(classification == "Misclassified") %>%
  mutate(misclass_type = case_when(
    (true == 1 & pred == 2) | (true == 2 & pred == 1) ~ "Misclassified Genotype",
    (true == 2 & pred == 0) | (true == 1 & pred == 0) ~ "False Negative",
    (true == 0 & pred == 2) | (true == 0 & pred == 1) ~ "False Positive",
    TRUE ~ "Other"  # Catch-all for any other misclassification (if any)
  ))


# Plot histogram of region lengths for misclassified TEs
p_misclass_hist <- ggplot(misclassified_data, aes(x = region_length)) +
  geom_histogram(bins = 30, fill = "tomato", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Misclassified TEs by Region Length",
       x = "Region Length",
       y = "Count") +
  theme_minimal()
print(p_misclass_hist)



# Plot histogram of region lengths for each misclassification type
p_misclass_hist_split <- ggplot(misclassified_data, aes(x = region_length)) +
  geom_histogram(bins = 30, fill = "tomato", color = "black", alpha = 0.7) +
  facet_wrap(~ misclass_type, scales = "free_y") +
  labs(title = "Histogram of Misclassified TEs by Region Length",
       x = "Region Length",
       y = "Count") +
  theme_minimal()

print(p_misclass_hist_split)


# Plot an overlaid histogram with proportions on the y-axis
p_misclass_overlay <- ggplot(misclassified_data, aes(x = region_length, fill = misclass_type)) +
  geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), 
                 bins = 10, position = "identity", alpha = 0.5) +
  labs(title = "Histogram of Misclassified TEs by Region Length",
       x = "Region Length",
       y = "Proportion",
       fill = "Misclassification Type") +
  theme_minimal()

print(p_misclass_overlay)

# Create or update the plot_group variable to classify each row
data <- data %>%
  mutate(plot_group = ifelse(classification == "Misclassified",
    case_when(
      (true == 1 & pred == 2) | (true == 2 & pred == 1) ~ "Misclassified Genotype",
      (true == 2 & pred == 0) | (true == 1 & pred == 0) ~ "False Negative",
      (true == 0 & (pred == 2 | pred == 1)) ~ "False Positive",
      TRUE ~ "Other"  # Catch-all for any unexpected misclassification types
    ),
    paste("True", true)
  ))

# Create the violin plot of region lengths for each group
p_violin <- ggplot(data, aes(x = plot_group, y = region_length, fill = plot_group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  labs(title = "Distribution of Region Lengths by Group",
       x = "Group",
       y = "Region Length") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tilt labels for readability

print(p_violin)

# Filter data for rows with true value 0 (true negatives and any misclassifications among true negatives)
tn_data <- data %>% 
  filter(true == 0)

# Create decile bins for region length (ignoring any NA values)
tn_data <- tn_data %>%
  mutate(length_bin = cut(region_length, 
                          breaks = quantile(region_length, probs = seq(0, 1, 0.1), na.rm = TRUE), 
                          include.lowest = TRUE))

# Calculate the accuracy for true negatives in each decile
# Here, accuracy is the proportion of rows correctly classified as true negatives (i.e., pred==0)
accuracy_by_decile <- tn_data %>%
  group_by(length_bin) %>%
  summarise(mean_accuracy = mean(correct),     # 'correct' should be 1 for true negatives (pred==0) and 0 otherwise
            mean_length = mean(region_length),
            count = n())

# Plot the accuracy vs. the mean region length of each decile
p_accuracy_tn <- ggplot(accuracy_by_decile, aes(x = mean_length, y = mean_accuracy)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Accuracy on True Negatives across Region Length Deciles",
       x = "Mean Region Length (per decile)",
       y = "Accuracy (Proportion of True Negatives)") +
  theme_minimal()

print(p_accuracy_tn)




# Bin the data into deciles based on region_length
data_binned <- data %>%
  mutate(length_bin = cut(region_length, 
                          breaks = quantile(region_length, probs = seq(0, 1, 0.1), na.rm = TRUE), 
                          include.lowest = TRUE)) %>%
  group_by(length_bin) %>%
  summarise(mean_length = mean(region_length, na.rm = TRUE),
            TP = sum((true %in% c(1, 2)) & (pred %in% c(1, 2))),
            FN = sum((true %in% c(1, 2)) & (pred == 0)),
            TN = sum((true == 0) & (pred == 0)),
            FP = sum((true == 0) & (pred %in% c(1, 2))),
            .groups = "drop") %>%
  mutate(sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA),
         specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA),
         balanced_accuracy = (sensitivity + specificity) / 2)

# Plot balanced accuracy vs. mean region length in each decile
p_balanced_acc <- ggplot(data_binned, aes(x = mean_length, y = balanced_accuracy)) +
  geom_point(size = 3) +
  geom_line() +
  labs(title = "Balanced Accuracy across Region Length Deciles",
       x = "Mean Region Length (per decile)",
       y = "Balanced Accuracy") +
  theme_minimal()

print(p_balanced_acc)

data %>% filter(true==0 & pred==2) 
