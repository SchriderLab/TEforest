# script to plot features from data to understand their distributions
library(tidyverse)
library(ggpubr)
library(kableExtra)

feature_data <- read.csv("/nas/longleaf/home/adaigle/TEforest/plots/predictions_with_features.csv")
feature_data <- feature_data[, !names(feature_data) %in% c("X", "insert_size")]

# Create 'read_length' column based on the 'file' prefix
feature_data <- feature_data %>%
  mutate(
    read_length = case_when(
      startsWith(file, "A2") ~ 54,
      startsWith(file, "AKA") ~ 125,
      startsWith(file, "JUT") ~ 154,
      TRUE ~ NA_real_  # Assign NA if no condition matches
    )
  )

feature_data <- feature_data %>%
  mutate(
    class = case_when(
      true == 0 & pred == 0 ~ "TN",
      true == 0 & pred == 1 ~ "FP_het",
      true == 0 & pred == 2 ~ "FP_homo",
      true == 1 & pred == 0 ~ "FN_het",
      true == 1 & pred == 1 ~ "TP_het",
      true == 1 & pred == 2 ~ "misclass_homo",
      true == 2 & pred == 0 ~ "FN_homo",
      true == 2 & pred == 1 ~ "misclass_het",
      true == 2 & pred == 2 ~ "TP_homo",
      TRUE ~ "Unknown"  # Catch-all for unexpected values
    )
  )

#feature_data_summary <- feature_data %>% group_by(true,read_length) %>% summarize(across(everything(), list(mean = mean, sd = sd)))
feature_data$true <- feature_data$true/2

feature_data$true <- as.factor(feature_data$true)
feature_data$read_length <- as.factor(feature_data$read_length)

feature_data$class <- factor(
  feature_data$class,
  levels = rev(c(
    "TP_homo",        # True positive - homozygous
    "misclass_homo",   # Misclassified as homozygous
    "FP_homo",        # False positive - homozygous
    "TP_het",         # True positive - heterozygous
    "misclass_het",  # Misclassified as heterozygous
    "FP_het",         # False positive - heterozygous
    "misclass_FN",    # Misclassified as false negative
    "FN_homo",
    "FN_het", 
    "TN"              # False negative
  )
))

TE_orphan_IQR_plt <- ggplot(data = feature_data, aes(x = TE_specific_Orphan_Read_IQR, y = true)) +
  geom_boxplot(aes(color = true), width = 0.7, show.legend = T) +
  labs(x = "TE_specific_Orphan_Read_IQR",
       y = "zygosity", title="") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

feature_data_filter <- feature_data %>%
  filter(TE_specific_Orphan_Read_IQR >= 0, TE_specific_Orphan_Read_IQR <= 1)

TE_orphan_IQR_plt_complex <- ggplot(data = feature_data_filter, aes(x = TE_specific_Orphan_Read_IQR, y = class)) +
  geom_boxplot(aes(color = class), width = 0.7, show.legend = TRUE) +
  labs(
    x = "TE_specific_Orphan_Read_IQR",
    y = "Class",
  ) +
  xlim(0, 1) +
  facet_wrap(~ read_length, ncol = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

feature_data_filter_150 <- feature_data_filter %>% filter(read_length==154)
feature_data_filter_150 <- feature_data_filter %>%
  filter(read_length == 154) %>%
  mutate(
    pred = case_when(
      pred == 0 ~ "Absent",
      pred == 1 ~ "Heterozygote",
      pred == 2 ~ "Homozygote",
      TRUE ~ as.character(pred)  # Retain original value if none of the conditions match
    )
  )

TE_orphan_IQR_plt_150 <- ggplot(data = feature_data_filter_150, aes(x = TE_specific_Orphan_Read_IQR, y = class)) +
  geom_boxplot(aes(color = factor(pred)), width = 0.7, show.legend = TRUE) +
  labs(
    x = "TE specific Discordant Read IQR",
    y = "Classification",
    color = "Predicted Zygosity"  # Updated legend title
  ) +
  xlim(0, 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/TE_orphan_IQR_plt_150.png", plot = TE_orphan_IQR_plt_150, width = 8.5, height = 8.5, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/TE_orphan_IQR_plt_150.svg", plot = TE_orphan_IQR_plt_150, width = 8.5, height = 8.5, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/TE_orphan_IQR_plt_150.tiff", plot = TE_orphan_IQR_plt_150, width = 8.5, height = 8.5, dpi = 300)



TE_orphan_IQR_plt <- ggplot(data = feature_data, aes(x = TE_specific_Orphan_Read_IQR, y = true)) +
  geom_boxplot(aes(color = true), width = 0.7, show.legend = T) +
  labs(x = "TE_specific_Orphan_Read_IQR",
       y = "zygosity", title="") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

feature_data_filter_paired <- feature_data %>%
  filter(Paired_mean >= 0.99, Paired_mean <= 1)

Paired_mean_plt <- ggplot(data = feature_data, aes(x = Paired_mean, y = true)) +
  geom_boxplot(aes(color = true), width = 0.7, show.legend = T) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "Paired_mean",
       y = "zygosity", title="") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

Paired_mean_plt_complex <- ggplot(data = feature_data, aes(x = Paired_sd, y = class)) +
  geom_boxplot(aes(color = class), width = 0.7, show.legend = TRUE) +
  labs(
    x = "Paired_mean",
    y = "Class",
  ) +
  xlim(0, .5) +
  facet_wrap(~ read_length, ncol = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

feature_data_filter_cigar4 <- feature_data %>%
  filter(Cigar4_sd >= 0, Cigar4_sd <= 1)

cigar4_plt <- ggplot(data = feature_data, aes(x = Cigar4_sd, y = true)) +
  geom_boxplot(aes(color = true), width = 0.7, show.legend = T) +
  scale_fill_manual(values = c("darkorange","purple","cyan4")) +
  labs(x = "TE_specific_Paired_mean",
       y = "zygosity", title="") +
  xlim(0, 1) +
  facet_wrap("read_length", ncol=1)

cigar4_plt_complex <- ggplot(data = feature_data_filter_paired, aes(x = Cigar4_sd, y = class)) +
  geom_boxplot(aes(color = class), width = 0.7, show.legend = TRUE) +
  labs(
    x = "Cigar4_sd",
    y = "Class",
  ) +
  xlim(0, 0.4) +
  facet_wrap(~ read_length, ncol = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

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

data <- read.table("/nas/longleaf/home/adaigle/TEforest/plots/feature_importances_all.txt", header=T, sep="\t")
data$Feature <- feature_list_extended[data$Feature + 1]  # +1 to account for 0-based indexing
#data <- data %>% filter(Built.in.Importance!=0)
# Select the top 25 features based on Permutation Importance
top_features <- data %>%
  arrange(desc(Permutation.Importance)) %>%
  head(10) %>%
  mutate(
    Permutation.Performance = sprintf("%.6f ± %.6f", Permutation.Importance, Permutation.Std)  # Format Importance ± Sd
  ) %>%
  select(
     Feature,
    `Permutation Performance` = Permutation.Performance,
    `Impurity-Based Importance` = Built.in.Importance
  )

# Create a nice table
kable(top_features, format = "markdown", caption = "Top 10 Features by Permutation Importance")# %>%
  #kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed")) %>%
  #column_spec(2, bold = TRUE)

# Normalize Importance values
data <- data %>%
  mutate(
    Built.in.Importance.Normalized = Built.in.Importance / sum(Built.in.Importance),
    Permutation.Importance.Normalized = Permutation.Importance / sum(Permutation.Importance)
  )

# Reshape data for side-by-side plotting (unnormalized)
data_long_unnormalized <- data %>%
  pivot_longer(cols = c("Built.in.Importance", "Permutation.Importance"), 
               names_to = "Importance.Type", 
               values_to = "Importance.Value") %>%
  mutate(
    ymin = ifelse(Importance.Type == "Permutation.Importance", 
                  Importance.Value - data$Permutation.Std, 
                  NA),
    ymax = ifelse(Importance.Type == "Permutation.Importance", 
                  Importance.Value + data$Permutation.Std, 
                  NA)
  )

# Reshape data for side-by-side plotting (normalized)
data_long_normalized <- data %>%
  pivot_longer(cols = c("Built.in.Importance.Normalized", "Permutation.Importance.Normalized"), 
               names_to = "Importance.Type", 
               values_to = "Importance.Value") %>%
  mutate(
    ymin = ifelse(Importance.Type == "Permutation.Importance.Normalized", 
                  Importance.Value - data$Permutation.Std / sum(data$Permutation.Importance), 
                  NA),
    ymax = ifelse(Importance.Type == "Permutation.Importance.Normalized", 
                  Importance.Value + data$Permutation.Std / sum(data$Permutation.Importance), 
                  NA)
  )

# Plot unnormalized values
plot_unnormalized <- ggplot(data_long_unnormalized, aes(x = reorder(Feature, Importance.Value), y = Importance.Value, fill = Importance.Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(data = data_long_unnormalized %>% filter(Importance.Type == "Permutation.Importance"),
                aes(ymin = ymin, ymax = ymax), 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  coord_flip() +
  labs(title = "Feature Importances (Unnormalized)", x = "Feature", y = "Importance") +
  scale_fill_manual(values = c("skyblue", "orange"), 
                    labels = c("Built-in Importance", "Permutation Importance")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

# Plot normalized values
plot_normalized <- ggplot(data_long_normalized, aes(x = reorder(Feature, Importance.Value), y = Importance.Value, fill = Importance.Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(data = data_long_normalized %>% filter(Importance.Type == "Permutation.Importance.Normalized"),
                aes(ymin = ymin, ymax = ymax), 
                position = position_dodge(width = 0.9), 
                width = 0.3) +
  coord_flip() +
  labs(title = "Feature Importances (Normalized)", x = "Feature", y = "Normalized Importance") +
  scale_fill_manual(values = c("skyblue", "orange"), 
                    labels = c("Built-in Importance", "Permutation Importance")) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

output_path <- "/nas/longleaf/home/adaigle/TEforest/plots/feat_importance_30X.png"
#ggsave(output_path, plot = plot, width = 10, height = 5)
