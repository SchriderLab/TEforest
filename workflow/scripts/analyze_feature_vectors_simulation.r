# Load necessary libraries
library(tidyverse)
library(stringr)
library(RcppCNPy)
library(randomForest)
library(corrplot)
library(GenomicRanges)
library(ggpubr)
# Define the base directory containing the feature vectors
base_dir <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/teforest_quantile_length2/feature_vectors"

# List all .npy files recursively
files <- list.files(base_dir, pattern = "\\.npy$", full.names = TRUE, recursive = TRUE)
files <- files[!grepl("_readcount", files)]
         
# Create a tibble and parse file names.
# The assumed file name pattern is:
#   <sample>-<chrom>-<start>-<end>-<feature_index>-<TE>.npy
# For example:
#   RL100IS200_rep0_rev-3L-4638-5361-0-roo.npy
df <- tibble(filepath = files) %>%
  mutate(
    filename = basename(filepath),
    # Extract components with regex:
    # Group 1: sample (e.g. RL100IS200_rep0_rev)
    # Group 2: chrom (e.g. 3L)
    # Group 3: start (e.g. 4638)
    # Group 4: end (e.g. 5361)
    # Group 5: feature index (ignored here)
    # Group 6: TE (e.g. roo)
    pattern = "^(.+)-([^-]+)-([0-9]+)-([0-9]+)-([0-9]+)-([^-]+)\\.npy$",
    sample_full = str_match(filename, pattern)[,2],
    chrom = str_match(filename, pattern)[,3],
    start = as.integer(str_match(filename, pattern)[,4]),
    end = as.integer(str_match(filename, pattern)[,5]),
    TE = str_match(filename, pattern)[,7],
    region_length = end - start
  ) %>%
  # Now extract RL, IS, replicate number, and orientation from the sample name.
  mutate(
    RL = str_extract(sample_full, "RL[0-9]+") %>% str_remove("RL") %>% as.integer(),
    IS = str_extract(sample_full, "IS[0-9]+") %>% str_remove("IS") %>% as.integer(),
    rep = str_extract(sample_full, "rep[0-9]+") %>% str_remove("rep") %>% as.integer(),
    orientation = str_extract(sample_full, "(fwd|rev)")
  ) %>%
  select(RL, IS, rep, orientation, chrom, start, end, region_length, filepath)

# Optionally, load the numpy feature vector for each file.
# The feature vector is stored in a list-column.
df <- df %>%
  rowwise() %>%
  mutate(features = list(npyLoad(filepath))) %>%
  ungroup()

# For clarity, group the tibble by RL, IS, and orientation.
df_grouped <- df %>%
  group_by(RL, IS, orientation)

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

# Create the extended feature list (each feature has mean, median, sd, IQR)
feature_list <- c()
for (x in feature_names) {
  feature_list <- c(feature_list, paste0(x, "_mean"), paste0(x, "_median"), paste0(x, "_sd"), paste0(x, "_IQR"))
}

# Extend the feature list with TE-specific features (adds another set of 84 features)
feature_list_extended <- c(feature_list, paste0("TE_specific_", feature_list))

# Now, assuming your df tibble has a column "features" where each element is a numeric vector of length 168,
# assign names to each vector and then unnest it into separate columns.
df_grouped <- df_grouped %>% 
  mutate(features = map(features, ~ set_names(.x, feature_list_extended))) %>% 
  unnest_wider(features)

# df now has 168 new columns with names from feature_list_extended.
ggplot(df_grouped, 
       aes(x = factor(IS), 
           y = TE_specific_Orphan_Read_IQR, 
           fill = factor(RL),
           group = interaction(IS, RL))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Insert Size",
       y = "TE Specific Orphan Read IQR",
       fill = "Read Length") +
  theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = TE_specific_Paired_median, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired Median",
           color = "Read length") +
      theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = Paired_mean, color = IS)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Paired_mean",
           color = "Insert Size") +
      theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = RL)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Orphan Read IQR",
           color = "Read length") +
      theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = TE_specific_Cigar4_mean, color = RL)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Cigar4_mean",
           color = "Read length") +
      theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = Quality_mean, color = RL)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Quality",
           color = "Read length") +
      theme_minimal()

ggplot(df_grouped, aes(x = IS, y = TE_specific_Proper_Pair_mean, color = RL)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Proper_Pair_mean",
           color = "Read length") +
      theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = Short_Insert_mean, color = IS)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Short_Insert_mean",
           color = "Insert size") +
      theme_minimal()

ggplot(df_grouped, aes(x = region_length, y = Insert_Size_mean, color = IS)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Short_Insert_mean",
           color = "Insert size") +
      theme_minimal()

# Make sure df_grouped is ungrouped for modeling
df_model <- df_grouped %>% ungroup()

# Build a random forest model for read length (RL)
rf_RL <- randomForest(RL ~ ., 
                      data = df_model %>% select(RL, one_of(feature_list_extended)),
                      importance = TRUE)
print(rf_RL)
varImpPlot(rf_RL, main = "Variable Importance for Read Length")

# Build a random forest model for insert size (IS)
rf_IS <- randomForest(IS ~ ., 
                      data = df_model %>% select(IS, one_of(feature_list_extended)),
                      importance = TRUE)
print(rf_IS)
varImpPlot(rf_IS, main = "Variable Importance for Insert Size")


df_model_filtered <- df_model
df_model_filtered$Quality_mean <- NULL
df_model_filtered$Quality_median <- NULL
df_model_filtered$Quality_sd <- NULL
df_model_filtered$Quality_IQR <- NULL
df_model_filtered$TE_specific_Quality_mean <- NULL
df_model_filtered$TE_specific_Quality_median <- NULL
df_model_filtered$TE_specific_Quality_sd <- NULL
df_model_filtered$TE_specific_Quality_IQR <- NULL
# Build a random forest model for region length
rf_region <- randomForest(region_length ~ ., 
                          data = df_model_filtered %>% select(region_length, one_of(feature_list_extended)),
                          importance = TRUE)
print(rf_region)
varImpPlot(rf_region, main = "Variable Importance for Region Length")

# Create a data frame with region_length and all feature columns,
# ensuring they are numeric.
# Select only the numeric columns (this will remove 'orientation')
# Select only the numeric columns (drop character columns like "orientation")

df_corr <- df_grouped %>% 
  select(region_length, all_of(feature_list_extended)) %>% 
  mutate(across(everything(), as.numeric))
df_corr_numeric <- df_corr 
df_corr_numeric$orientation <- NULL
# Check structure to confirm only numeric columns remain
# Compute correlation matrix (selecting only numeric columns)
corr_matrix <- cor(df_corr_numeric, use = "pairwise.complete.obs")

# Build a two-column tibble from the 'region_length' correlations:
region_cor <- tibble(
  Feature = names(corr_matrix["region_length",]),
  Correlation = as.numeric(corr_matrix["region_length",])
) %>%
  filter(Feature != "region_length") %>%  # Remove self-correlation
  drop_na(Correlation) %>%                # Remove any NA correlations
  arrange(desc(Correlation))              # Sort descending by Correlation

# Plot using ggplot2:
ggplot(region_cor, aes(x = reorder(Feature, Correlation), y = Correlation, fill = Correlation)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Correlation of Features with Region Length",
       x = "Feature",
       y = "Correlation Coefficient") +
  theme_minimal()


### Now we load in features from training data to compare to simulated data
# Define the directory with your files:
data_dir <- "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/3L3RX"
#data_dir <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/teforest_rescale_length2/feature_vectors"
#data_dir <- "/nas/longleaf/home/adaigle/work/test_TEforest/test_celegans_fullnorm/feature_vectors/RW7000"

# List files that have "-2-" in the name and end with .npy
files <- list.files(data_dir, pattern = "-2-.*\\.npy$", full.names = TRUE)
#files <- list.files(data_dir, pattern = "-0-.*\\.npy$", full.names = TRUE, recursive = TRUE)

# Define a key mapping sample prefixes to RL and IS values:
rl_is_key <- tibble(
  sample_prefix = c("A2_A3", "AKA-017_GIM-024", "JUT-008_MUN-009"),
  RL = c(54, 124, 151),
  IS = c(287, 208, 436)
)

# Build the data frame with metadata extracted from the filename.
# The new file names have this structure:
#   sample_prefix - chrom - start - end - feature_index - TEelement.npy
# For example: A2_A3-X-13276022-13276697-2-S_element.npy
df_empirical <- tibble(filepath = files) %>%
  mutate(
    filename = basename(filepath),
    # The regex below uses 6 capture groups:
    #   Group 1: sample_prefix (e.g. A2_A3)
    #   Group 2: chrom (e.g. X)
    #   Group 3: start (e.g. 13276022)
    #   Group 4: end (e.g. 13276697)
    #   Group 5: feature index (e.g. 2; here always "2")
    #   Group 6: TE type (e.g. S_element)
    pattern = "^(.+)-([^-]+)-([0-9]+)-([0-9]+)-([0-9]+)-([^-]+)\\.npy$",
    #pattern = "^A2_A3-([^-]+)-([0-9]+)-([0-9]+)-([0-9]+)-([^-]+)\\.npy$",

    #pattern = "^([^-]+\\.modref)-([^-]+)-([0-9]+)-([0-9]+)-([0-9]+)-([^-]+)\\.npy$",
    sample_prefix = str_match(filename, pattern)[,2],
    #sample_prefix = "A2_A3",
    chrom = str_match(filename, pattern)[,3],
    start = as.integer(str_match(filename, pattern)[,4]),
    end = as.integer(str_match(filename, pattern)[,5]),
    TE = str_match(filename, pattern)[,6],
    region_length = end - start
  ) %>%
  # Add the RL and IS columns by joining with the key table:
  left_join(rl_is_key, by = "sample_prefix") %>%
  select(RL, IS, sample_prefix, chrom, start, end, region_length, TE, filepath)

#df_empirical <- df_empirical %>%
#  mutate(feature_vector = map(filepath, ~ npyLoad(.x)))

df_empirical <- df_empirical %>%
  rowwise() %>%
  mutate(features = list(npyLoad(filepath))) %>%
  ungroup()

# For clarity, group the tibble by RL, IS, and orientation.
df_empirical_grouped <- df_empirical %>%
  group_by(RL, IS)

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

# Create the extended feature list (each feature has mean, median, sd, IQR)
feature_list <- c()
for (x in feature_names) {
  feature_list <- c(feature_list, paste0(x, "_mean"), paste0(x, "_median"), paste0(x, "_sd"), paste0(x, "_IQR"))
}

# Extend the feature list with TE-specific features (adds another set of 84 features)
feature_list_extended <- c(feature_list, paste0("TE_specific_", feature_list))

# Now, assuming your df tibble has a column "features" where each element is a numeric vector of length 168,
# assign names to each vector and then unnest it into separate columns.
df_empirical_grouped <- df_empirical_grouped %>% 
  mutate(features = map(features, ~ set_names(.x, feature_list_extended))) %>% 
  unnest_wider(features)

df_empirical_grouped <- df_empirical_grouped %>%
  mutate(IS_sd = case_when(
    sample_prefix == "A2_A3" ~ 50,
    sample_prefix == "JUT-008_MUN-009" ~ 113,
    sample_prefix == "AKA-017_GIM-024" ~ 125,
    TRUE ~ NA_real_  # Assigns NA to any other sample
  ))

df_empirical_grouped_roo <- df_empirical_grouped %>% filter(TE=="roo")


ggplot(df_empirical_grouped, 
       aes(x = factor(IS), 
           y = TE_specific_Orphan_Read_IQR, 
           fill = factor(RL),
           group = interaction(IS, RL))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Insert Size",
       y = "TE Specific Orphan Read IQR",
       fill = "Read Length") +
  theme_minimal()

ggplot(df_empirical_grouped, aes(x = region_length, y = TE_specific_Paired_median, color = IS)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired Median",
           color = "Insert Size") +
      theme_minimal()

ggplot(df_empirical_grouped, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = IS)) +
      geom_point(alpha = 0.8) +
      ylim(0,1) +
      labs(x = "Region Length",
           y = "TE_specific_Orphan_Read_IQR",
           color = "Insert Size") +
      theme_minimal()

ggplot(df_empirical_grouped, aes(x = region_length, y = TE_specific_Orphan_Read_median, color = IS)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Orphan_Read_median",
           color = "Insert Size") +
      theme_minimal()

ggplot(df_empirical_grouped, aes(x = region_length, y = TE_specific_Cigar4_mean, color = IS)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Cigar4_mean",
           color = "Insert Size") +
      theme_minimal() + ylim(0,1)

df_grouped$IS_sd <- 100
combined_df <- bind_rows(df_empirical_grouped, df_grouped)



ggplot(combined_df, aes(x = region_length, y = TE_specific_Orphan_Read_median, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Orphan_Read_median",
           color = "Read length") +
      theme_minimal() 

ggplot(combined_df, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Orphan Read IQR",
           color = "Read length") +
      theme_minimal() + ylim(0,1)

ggplot(combined_df, aes(x = region_length, y = TE_specific_Paired_median, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired Median",
           color = "Insert Size") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = TE_specific_Paired_IQR, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired IQR",
           color = "Read Length") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = Cigar4_mean, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Cigar4_mean",
           color = "Read Length") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = Cigar4_IQR, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Cigar4_IQR",
           color = "Read Length") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = Paired_mean, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      ylim(0.98,1) +
      labs(x = "Region Length",
           y = "Paired_mean",
           color = "Read Length") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = Paired_sd, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Paired_sd",
           color = "Read Length") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = Proper_Pair_mean, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Proper_Pair_mean",
           color = "Read Length") +
      theme_minimal()

ggplot(combined_df, aes(x = region_length, y = Proper_Pair_sd, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Proper_Pair_sd",
           color = "Read Length") +
      theme_minimal()


##plotting sd of insert size against various features
ggplot(combined_df, aes(x = IS_sd, y = TE_specific_Paired_median, color = RL)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired Median",
           color = "Insert Size") +
      theme_minimal()

ggplot(combined_df, aes(x = IS_sd, y = TE_specific_Paired_IQR, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired IQR",
           color = "Read Length") +
      theme_minimal()

feature_columns <- combined_df %>% select(10:ncol(.))
feature_columns$RL <- NULL
feature_columns$IS <- NULL
feature_columns$IS_sd <- NULL

# Remove non-numeric columns
combined_df_cleaned <- combined_df %>%
  select(where(is.numeric))  # Keep only numeric columns

# Ensure `IS_sd` is in the dataset (if it was removed by accident)
if (!"IS_sd" %in% colnames(combined_df_cleaned)) {
  combined_df_cleaned <- combined_df_cleaned %>%
    left_join(combined_df %>% select(sample, IS_sd), by = "sample")
}

# Select numeric feature columns (excluding RL, IS, and IS_sd, which are predictors)
feature_columns <- setdiff(colnames(combined_df_cleaned)[10:ncol(combined_df_cleaned)], c("RL", "IS", "IS_sd"))

# Fit linear models for each feature
lm_results <- lapply(feature_columns, function(feature) {
  # Create a subset of data removing NAs and Inf values
  df_filtered <- combined_df_cleaned %>%
    select(all_of(feature), RL, IS, IS_sd) %>%
    filter(!is.na(.data[[feature]]), !is.infinite(.data[[feature]]),
           !is.na(RL), !is.na(IS), !is.na(IS_sd))  # Remove NA/Inf values
  
  # Ensure there are enough data points
  if (nrow(df_filtered) < 2) return(NA)  # Skip feature if not enough data

  # Fit the model
  formula <- as.formula(paste0("`", feature, "` ~ IS_sd"))  # Handle special characters
  model <- lm(formula, data = df_filtered)

  # Extract R-squared value
  return(summary(model)$r.squared)
})

# Convert to dataframe
lm_results_df <- data.frame(Feature = feature_columns, R_squared = unlist(lm_results))

# Remove NA results
lm_results_df <- lm_results_df %>% filter(!is.na(R_squared))

# Sort results
lm_results_df <- lm_results_df %>% arrange(desc(R_squared))

# Print results
print(lm_results_df)


ggplot(lm_results_df, aes(x = reorder(Feature, -R_squared), y = R_squared)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "R² of Linear Models for Each Feature",
       x = "Feature",
       y = "R² Value")

# Filter only features with "IQR" or "sd" in the name
lm_results_filtered <- lm_results_df %>%
  filter(grepl("IQR|sd", Feature, ignore.case = TRUE))

# Plot the R² values for the selected features
ggplot(lm_results_filtered, aes(x = reorder(Feature, -R_squared), y = R_squared)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "R² of Linear Models for IQR and SD Features",
       x = "Feature",
       y = "R² Value")

selected_features <- c(
  "TE_specific_Orphan_Read_IQR",
  "Paired_mean",
  "Paired_sd",
  "TE_specific_Orphan_Read_mean",
  "Proper_Pair_mean",
  "Cigar4_sd",
  "Cigar1_sd",
  "Quality_IQR",
  "Proper_Pair_sd",
  "Orphan_Read_median"
)

# Filter the dataset for only the selected features
lm_results_selected <- lm_results_df %>%
  filter(Feature %in% selected_features)

# Create the plot
ggplot(lm_results_selected, aes(x = reorder(Feature, -R_squared), y = R_squared)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "R² of Linear Models for Selected Features",
       x = "Feature",
       y = "R² Value") +
  theme(axis.text.y = element_text(size = 10))


# Set the R² threshold
R2_cutoff <- 0.3

# Extract feature names with R² < 0.3
low_correlation_features <- lm_results_df %>%
  filter(R_squared < R2_cutoff) %>%
  pull(Feature)  # Extract feature names as a character vector
low_correlation_features <- low_correlation_features[!grepl("rep", low_correlation_features, ignore.case = TRUE)]

# Print the character vector
print(low_correlation_features)
## train RF to predict genotype, to quickly test if removing features is viable
# List files that have "-2-" in the name and end with .npy
files_train <- list.files(data_dir, pattern = "*\\.npy$", full.names = TRUE)

# Build the data frame with metadata extracted from the filename.
# The new file names have this structure:
#   sample_prefix - chrom - start - end - feature_index - TEelement.npy
# For example: A2_A3-X-13276022-13276697-2-S_element.npy
df_empirical_train <- tibble(filepath = files_train) %>%
  mutate(
    filename = basename(filepath),
    # The regex below uses 6 capture groups:
    #   Group 1: sample_prefix (e.g. A2_A3)
    #   Group 2: chrom (e.g. X)
    #   Group 3: start (e.g. 13276022)
    #   Group 4: end (e.g. 13276697)
    #   Group 5: feature index (e.g. 2; here always "2")
    #   Group 6: TE type (e.g. S_element)
    pattern = "^(.+)-([^-]+)-([0-9]+)-([0-9]+)-([0-9]+)-([^-]+)\\.npy$",
    sample_prefix = str_match(filename, pattern)[,2],
    chrom = str_match(filename, pattern)[,3],
    start = as.integer(str_match(filename, pattern)[,4]),
    end = as.integer(str_match(filename, pattern)[,5]),
    genotype = str_match(filename, pattern)[,6],
    TE = str_match(filename, pattern)[,7],
    region_length = end - start
  ) %>%
  # Add the RL and IS columns by joining with the key table:
  left_join(rl_is_key, by = "sample_prefix") %>%
  select(RL, IS, sample_prefix, chrom, start, end, region_length, TE, filepath,genotype)


df_empirical_train <- df_empirical_train %>%
  rowwise() %>%
  mutate(features = list(npyLoad(filepath))) %>%
  ungroup()

df_empirical_train <- df_empirical_train %>% 
  mutate(features = map(features, ~ set_names(.x, feature_list_extended))) %>% 
  unnest_wider(features)


# Subset the training data to only include 'genotype' and the low correlation feature columns
#train_data <- df_empirical_train[, c("genotype", low_correlation_features)]
train_data <- df_empirical_train %>% 
  mutate(IS_sd = case_when(
    sample_prefix == "A2_A3" ~ 50,
    sample_prefix == "JUT-008_MUN-009" ~ 113,
    sample_prefix == "AKA-017_GIM-024" ~ 125,
    TRUE ~ NA_real_  # Assigns NA to any other sample
  )) %>% select(-sample_prefix, -chrom, -start, -end, -TE, -filepath)

# Ensure the response variable 'genotype' is a factor
train_data$genotype <- as.factor(train_data$genotype)

# Train the random forest model; importance=TRUE to extract variable importance measures
rf_model <- randomForest(genotype ~ ., data = train_data, importance = TRUE)

# Print the model summary
print(rf_model)

# Manually create confusion matrix from provided data
conf_matrix <- matrix(c(18657, 57, 28,
                        169,    608, 16,
                        45,    9, 526),
                      nrow=3, byrow=TRUE)

rownames(conf_matrix) <- colnames(conf_matrix) <- c("0", "1", "2")

# Precision, Recall, F1 calculations
calc_metrics <- function(conf_matrix, class_label){
  TP <- conf_matrix[class_label, class_label]
  FP <- sum(conf_matrix[-class_label, class_label])
  FN <- sum(conf_matrix[class_label, -class_label])

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * (precision * recall) / (precision + recall)

  return(c(Precision=precision, Recall=recall, F1=f1))
}

metrics <- sapply(1:3, function(x) calc_metrics(conf_matrix, x))
colnames(metrics) <- c("0","1","2")

# Display results
print(metrics)
# View variable importance
print(importance(rf_model))

# Convert the importance matrix to a data frame
importance_df <- as.data.frame(importance(rf_model))

# Arrange rows by MeanDecreaseAccuracy (descending)
importance_df_sorted <- importance_df %>%
  arrange(desc(MeanDecreaseGini))


#df_grouped <- df_grouped %>% ungroup() %>%
#  mutate(
#    RL = as.numeric(RL),
#    IS = as.numeric(IS),
#    IS_sd = as.numeric(IS_sd)
#  )

df_sim_test <- df_grouped %>% ungroup() %>% select(-rep, -orientation, -chrom, -start, -end, -filepath)
# test on simulations
df_sim_test <- df_sim_test %>%
  mutate(predicted_genotype = predict(rf_model, newdata = .))

# Calculate overall accuracy assuming the true genotype should be "2"
overall_accuracy <- mean(df_sim_test$predicted_genotype == "2")
cat("Overall accuracy (proportion predicted as '2'): ", overall_accuracy, "\n")

# Now, calculate accuracy across RL, IS, and IS_sd groups
accuracy_by_group <- df_sim_test %>%
  group_by(RL, IS, IS_sd) %>%
  summarise(n = n(),
            accuracy = mean(predicted_genotype == "2"),
            .groups = "drop")

print(accuracy_by_group)

# if on sep session run this based on curr accuracies of rf model

# Create the lookup table from your provided data
accuracy_lookup <- read.table(text = "
RL   IS  IS_sd  n  accuracy
50  350   10   20 1.00
75  300   10   20 1.00
75  350   10   20 1.00
75  400   10   17 1.00
75  450   10   20 1.00
75  500   10    7 1.00
100 300   10   20 1.00
100 350   10    1 1.00
100 400   10   20 1.00
100 450   10   20 1.00
100 500   10   20 1.00
125 400   10   20 1.00
125 450   10   20 1.00
125 500   10   20 1.00
150 350   10   20 1.00
150 400   10   18 1.00
150 450   10   20 1.00
150 500   10   20 1.00
50  400   10   20 0.95
125 300   10   20 0.95
50  300   10   20 0.90
50  450   10   20 0.90
50  500   10   20 0.90
75  250   10   20 0.80
100 250   10   20 0.60
150 300   10   20 0.60
50  250   10   20 0.55
50  200   10   20 0.30
75  200   10   20 0.30
125 250   10   20 0.25
100 200   10   20 0.00
125 200   10   20 0.00
150 200   10    5 0.00
150 250   10   20 0.00
", header = TRUE)


# Now join the lookup table to your combined_df based on RL, IS, and IS_sd
combined_df_accuracy <- combined_df %>%
  left_join(accuracy_lookup %>% select(RL, IS, IS_sd, accuracy),
            by = c("RL", "IS", "IS_sd"))


ggplot(combined_df_accuracy, aes(x = region_length, y = TE_specific_Orphan_Read_median, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Orphan Read median",
           color = "Accuracy") +
      scale_color_gradient(low = "red", high = "green") +

      theme_minimal()

ggplot(combined_df_accuracy, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Orphan Read IQR",
           color = "Accuracy") +
      theme_minimal()

ggplot(combined_df_accuracy, aes(x = region_length, y = accuracy, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Accuracy",
           color = "Accuracy")


## which features are correlated with accuracy? 
df_corr <- combined_df_accuracy %>% 
  select(accuracy, all_of(feature_list_extended)) %>% 
  mutate(across(everything(), as.numeric))
df_corr_numeric <- df_corr 
df_corr_numeric$orientation <- NULL
# Check structure to confirm only numeric columns remain
# Compute correlation matrix (selecting only numeric columns)
corr_matrix <- cor(df_corr_numeric, use = "pairwise.complete.obs")

# Build a two-column tibble from the 'region_length' correlations:
region_cor <- tibble(
  Feature = names(corr_matrix["accuracy",]),
  Correlation = as.numeric(corr_matrix["accuracy",])
) %>%
  filter(Feature != "accuracy") %>%  # Remove self-correlation
  drop_na(Correlation) %>%                # Remove any NA correlations
  arrange(desc(Correlation))              # Sort descending by Correlation

# Plot using ggplot2:
ggplot(region_cor, aes(x = reorder(Feature, Correlation), y = Correlation, fill = Correlation)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Correlation of Features with Region Length",
       x = "Feature",
       y = "Correlation Coefficient") +
  theme_minimal() + theme(axis.text.y=element_text(size=5))


ggplot(combined_df_accuracy, aes(x = region_length, y = Orphan_Read_median, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Orphan Read median",
           color = "Accuracy")

ggplot(combined_df_accuracy, aes(x = region_length, y = Short_Insert_IQR, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Orphan Read median",
           color = "Accuracy")

filter_values <- data.frame(
  RL = c(50, 50, 50, 75, 100, 150, 50, 50, 75, 125),
  IS = c(300, 450, 500, 250, 250, 300, 250, 200, 200, 250)
)

# Filter the combined_df_accuracy to only include rows with matching RL and IS
combined_df_accuracy_intermediate <- merge(combined_df_accuracy, filter_values, by = c("RL", "IS"))

## which features are correlated with accuracy? 
df_corr <- combined_df_accuracy_intermediate %>% 
  select(accuracy, all_of(feature_list_extended)) %>% 
  mutate(across(everything(), as.numeric))
df_corr_numeric <- df_corr 
df_corr_numeric$orientation <- NULL
# Check structure to confirm only numeric columns remain
# Compute correlation matrix (selecting only numeric columns)
corr_matrix <- cor(df_corr_numeric, use = "pairwise.complete.obs")

# Build a two-column tibble from the 'region_length' correlations:
region_cor <- tibble(
  Feature = names(corr_matrix["accuracy",]),
  Correlation = as.numeric(corr_matrix["accuracy",])
) %>%
  filter(Feature != "accuracy") %>%  # Remove self-correlation
  drop_na(Correlation) %>%                # Remove any NA correlations
  arrange(desc(Correlation))              # Sort descending by Correlation

ggplot(region_cor, aes(x = reorder(Feature, Correlation), y = Correlation, fill = Correlation)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Correlation of Features with Region Length",
       x = "Feature",
       y = "Correlation Coefficient") +
  theme_minimal() + theme(axis.text.y=element_text(size=5))

ggplot(combined_df_accuracy, aes(x = region_length, y = Cigar4_mean, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Cigar4_mean",
           color = "Accuracy")

ggplot(combined_df_accuracy_intermediate, aes(x = region_length, y = Quality_mean, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Quality_mean",
           color = "Accuracy")

# attempt to normalize by power law... 
ggplot(combined_df_accuracy, aes(x = region_length, y = Cigar4_mean*region_length, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Cigar4_mean X region length",
           color = "Accuracy")

ggplot(combined_df_accuracy, aes(x = region_length, y = Orphan_Read_median*region_length, color = accuracy)) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Orphan Read median X insert size",
           color = "Accuracy")


### laod in worms
data_dir <- "/nas/longleaf/home/adaigle/work/test_TEforest/test_celegans_fullnorm/feature_vectors/RW7000"

data_dir <- "/nas/longleaf/home/adaigle/TEforest/workflow/scripts/featvec_csvs/RW7000"

data_dir <- "/nas/longleaf/home/adaigle/users/rob_flies/teforest_runs/fullMA_fullnorm/feature_vectors/MA19F"

# Load files for each genotype and add the genotype column
# 1. Read and parse the predictions file
##predictions <- read_csv("/nas/longleaf/home/adaigle/work/test_TEforest/test_celegans_fullnorm/output/RW7000/predictions.csv") %>%
##  # Split the file column into components:
##  separate(file, into = c("sample", "chrom", "start", "stop", "nothing", "TE"),
##           sep = "-", convert = TRUE)
predictions <- read_csv("/nas/longleaf/home/adaigle/users/rob_flies/teforest_runs/fullMA_fullnorm/output/MA19F/predictions.csv") %>%
  # Split the file column into components:
  separate(file, into = c("sample", "chrom", "start", "stop", "nothing", "TE"),
           sep = "-", convert = TRUE)
# 2. Load your files (all files have "-0-" in the name)
files <- list.files(data_dir, pattern = "-0-.*\\.npy$", full.names = TRUE, recursive = TRUE)
df_all <- tibble(filepath = files)

# 3. Define the key mapping for RL and IS values
rl_is_key <- tibble(
  sample_prefix = c("A2_A3", "AKA-017_GIM-024", "JUT-008_MUN-009"),
  RL = c(54, 124, 151),
  IS = c(287, 208, 436)
)

# 4. Build the empirical data frame by extracting metadata from each filename.
#    The filename structure is: sample_prefix-chrom-start-end-feature_index-TEelement.npy
df_elegans <- df_all %>%
  mutate(
    filename = basename(filepath),
    pattern = "^(.+)-([^-]+)-([0-9]+)-([0-9]+)-([0-9]+)-([^-]+)\\.npy$",
    sample_prefix = str_match(filename, pattern)[,2],
    chrom = str_match(filename, pattern)[,3],
    start = as.integer(str_match(filename, pattern)[,4]),
    end = as.integer(str_match(filename, pattern)[,5]),
    TE = str_match(filename, pattern)[,7],
    region_length = end - start
  ) %>%
  left_join(rl_is_key, by = "sample_prefix") %>%
  select(sample_prefix, chrom, start, end, region_length, TE, filepath, RL, IS)

# 5. Load the feature vectors from each file
df_elegans <- df_elegans %>%
  rowwise() %>%
  mutate(features = list(npyLoad(filepath))) %>%
  ungroup()

# Group by RL and IS (for clarity in subsequent operations)
df_elegans_grouped <- df_elegans# %>% group_by(RL, IS)

# 6. Define the feature names and create the extended feature list.
feature_names <- c(
  "Cigar1", "Cigar2", "Cigar3", "Cigar4", "Cigar5",
  "Paired", "Proper_Pair", "Is_Read1_Unmapped", "Is_Read2_Unmapped",
  "Is_Read1_Rev_Comp", "Is_Read2_Rev_Comp", "Is_First_Read", "Is_Second_Read",
  "Split", "Long_Insert", "Short_Insert", "Parallel_Read", "Everted_Read",
  "Orphan_Read", "Insert_Size", "Quality"
)

feature_list <- c()
for (x in feature_names) {
  feature_list <- c(feature_list, paste0(x, "_mean"), paste0(x, "_median"),
                    paste0(x, "_sd"), paste0(x, "_IQR"))
}
feature_list_extended <- c(feature_list, paste0("TE_specific_", feature_list))

# Unnest the features into separate columns
df_elegans_grouped <- df_elegans_grouped %>% 
  mutate(features = map(features, ~ set_names(.x, feature_list_extended))) %>% 
  unnest_wider(features)

# Optionally, add an IS_sd column based on sample_prefix
df_elegans_grouped <- df_elegans_grouped %>%
  mutate(IS_sd = case_when(
    sample_prefix == "A2_A3" ~ 50,
    sample_prefix == "JUT-008_MUN-009" ~ 113,
    sample_prefix == "AKA-017_GIM-024" ~ 125,
    TRUE ~ NA_real_
  ))

# 7. Join the predictions with the empirical data.
#    We match on: chrom, start, TE, and end (which corresponds to 'stop' in predictions)
df_elegans_grouped <- df_elegans_grouped %>%
  left_join(predictions %>% select(chrom, start, stop, TE, pred),
            by = c("chrom", "start", "TE", "end" = "stop"))

df_elegans_grouped <- df_elegans_grouped %>% ungroup()
# Rename the 'pred' column to 'genotype'
df_elegans_grouped <- df_elegans_grouped %>%
  rename_with(~ "genotype", .cols = "pred")



#rf_geno <- randomForest(genotype ~ ., 
#                      data = df_elegans_grouped %>% select(genotype, one_of(feature_list_extended)),
#                      importance = TRUE)
#print(rf_geno)
#varImpPlot(rf_geno, main = "Variable Importance for genotype")


# 8. Plot using genotype to color the points.
ggplot(df_elegans_grouped, aes(x = region_length, y = TE_specific_Orphan_Read_IQR,
                                 color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  ylim(0, 1) +
  xlim(0,3000)
  labs(x = "Region Length",
       y = "TE_specific_Orphan_Read_IQR",
       color = "Genotype") +
  theme_minimal()

# Combine into two panels (side by side)
ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))

ggplot(df_elegans_grouped, aes(x = region_length, y = Proper_Pair_sd,
                                 color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  ylim(0, .3) +
  xlim(0,3000) +
  labs(x = "Region Length",
       y = "Proper_Pair_sd",
       color = "Genotype") +
  theme_minimal()



# Path to your BED file
bed_path <- "/nas/longleaf/home/adaigle/users/elegans_mcclintock/new_consensus/RW7000_1/results/temp2/RW7000_1_temp2_nonredundant.bed"
bed_path <- "/nas/longleaf/home/adaigle/users/rob_flies/trimmed_reads_test/MA19F_50_1/results/temp2/MA19F_50_1_temp2_nonredundant.bed"

# Read the BED file, skipping the first header line
bed_data <- read_tsv(bed_path, col_names = FALSE, skip = 1)

# Parse and separate fields
parsed_bed <- bed_data %>%
  select(chrom = X1, start = X2, stop = X3, info = X4, strand = X6) %>%
  separate(info, into = c("TE", "type", "frequency", "sample", "caller", "rp", "num"), sep = "\\|") %>%
  mutate(
    start = as.integer(start),
    stop = as.integer(stop),
    frequency = as.numeric(frequency),
    num = as.integer(num)
  ) %>% filter(type=="non-reference")

# 1. Convert parsed_bed to GRanges
gr_bed <- GRanges(
  seqnames = parsed_bed$chrom,
  ranges = IRanges(start = parsed_bed$start, end = parsed_bed$stop)
)

# 2. Convert df_elegans_grouped to GRanges
gr_elegans <- GRanges(
  seqnames = df_elegans_grouped$chrom,
  ranges = IRanges(start = df_elegans_grouped$start, end = df_elegans_grouped$end)
)

# 3. Find overlaps between gr_elegans and gr_bed
hits <- findOverlaps(gr_elegans, gr_bed)

# 4. Create an indicator vector for overlaps with matching TE
overlap_flag <- rep(0, length(gr_elegans))
hits_df <- data.frame(query = queryHits(hits), subject = subjectHits(hits))

# Loop over each unique element in gr_elegans that has an overlap
for (i in unique(hits_df$query)) {
  # Get all corresponding subject indices from parsed_bed that overlap with gr_elegans[i]
  subject_indices <- hits_df$subject[hits_df$query == i]
  
  # Check if any overlapping element has a matching TE
  # (Assuming that the TE column exists in both data frames)
  if (any(df_elegans_grouped$TE[i] == parsed_bed$TE[subject_indices])) {
    overlap_flag[i] <- 1
  }
}

# 5. Add the new overlap column to the original data frame
df_elegans_grouped <- df_elegans_grouped %>%
  mutate(overlaps_bed = overlap_flag)

df_elegans_grouped %>% filter(overlaps_bed==1) %>% filter(genotype==0)


# Define a subset for the special red-highlighted points
highlight_df <- df_elegans_grouped %>%
  filter(overlaps_bed == 1 & genotype == 0)

# Base plot style to reuse
base_theme <- theme_minimal()
xlim_range <- c(0, 3000)

ggplot(df_elegans_grouped, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = TE_specific_Orphan_Read_IQR),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "TE_specific_Insert_Size_mean", color = "Genotype") +
  base_theme

# PLOT 1
p1 <- ggplot(df_elegans_grouped, aes(x = region_length, y = TE_specific_Insert_Size_mean, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = TE_specific_Insert_Size_mean),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 100) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "TE_specific_Insert_Size_mean", color = "Genotype") +
  base_theme

# PLOT 2
p2 <- ggplot(df_elegans_grouped, aes(x = region_length, y = TE_specific_Insert_Size_sd, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = TE_specific_Insert_Size_sd),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 100) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "TE_specific_Insert_Size_sd", color = "Genotype") +
  base_theme

# PLOT 3
p3 <- ggplot(df_elegans_grouped, aes(x = region_length, y = TE_specific_Cigar1_mean, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = TE_specific_Cigar1_mean),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "TE_specific_Cigar1_mean", color = "Genotype") +
  base_theme

# PLOT 4
p4 <- ggplot(df_elegans_grouped, aes(x = region_length, y = Cigar5_median, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = Cigar5_median),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "Cigar5_median", color = "Genotype") +
  base_theme

# Combine all plots into 2x2 grid
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

ggplot(df_elegans_grouped, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = Cigar5_median),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "Cigar5_median", color = "Genotype") +
  base_theme


# Make newgenotype a factor (classifier labels)
df_elegans_grouped <- df_elegans_grouped %>%
  mutate(newgenotype = factor(if_else(overlaps_bed == 1 & genotype == 0, 3, genotype)))

df_elegans_grouped_filter <- df_elegans_grouped %>% filter(newgenotype!=0)


df_elegans_grouped_filter <- df_elegans_grouped_filter %>%
  mutate(newgenotype = factor(if_else(overlaps_bed == 1 & genotype == 0, 1, 0)))

# 2. Check class distribution
table(df_elegans_grouped_filter$newgenotype)


# 3. Determine balanced sampling sizes:
class_counts <- table(df_elegans_grouped_filter$newgenotype)
min_class <- min(class_counts)
sampsize <- c(`0` = min_class, `1` = min_class)
# This ensures that for each tree, both classes are represented equally.

# 4. Train the random forest classifier with balanced sampling.
rf_newgenotype <- randomForest(
  newgenotype ~ ., 
  data = df_elegans_grouped_filter %>% select(newgenotype, one_of(feature_list_extended)),
  importance = TRUE,
  sampsize = sampsize
)

# Print the model summary:
print(rf_newgenotype)

# 5. Plot variable importance:
varImpPlot(rf_newgenotype, main = "Variable Importance: Group 3 vs Others")


ggplot(df_elegans_grouped, aes(x = region_length, y = Insert_Size_sd, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = Insert_Size_sd),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "Cigar5_median", color = "Genotype") +
  base_theme

ggplot(df_elegans_grouped, aes(x = region_length, y = Insert_Size_IQR, color = factor(genotype))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = Insert_Size_IQR),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "Insert_Size_IQR", color = "Genotype") +
  base_theme

ggplot(df_elegans_grouped, aes(x = region_length, y = Cigar4_sd, color = factor(region_length))) +
  geom_point(alpha = 0.8) +
  geom_point(data = highlight_df, aes(x = region_length, y = Cigar4_sd),
             color = "red", size = 3, alpha = 1, shape = 21, stroke = 1.5) +
  #ylim(0, 1) +
  xlim(xlim_range) +
  labs(x = "Region Length", y = "Cigar4_sd", color = "Genotype") +
  base_theme


df_elegans_grouped$RL <- 125
df_elegans_grouped$IS <- 361
elegans_training_df <- bind_rows(df_empirical_grouped, df_elegans_grouped)


ggplot(elegans_training_df, aes(x = region_length, y = TE_specific_Orphan_Read_median, color = factor(genotype))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Orphan_Read_median",
           color = "Genotype") +
      theme_minimal() +ylim(0,1)

ggplot(elegans_training_df, aes(x = factor(genotype), y = TE_specific_Orphan_Read_median, fill = factor(genotype))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(x = "Genotype",
       y = "TE Specific Orphan Read (Median)",
       fill = "Genotype") +
  theme_minimal() +
  ylim(0, 1)

ggplot(elegans_training_df, aes(x = region_length, y = TE_specific_Orphan_Read_IQR, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Orphan Read IQR",
           color = "Read length") +
      theme_minimal() + ylim(0,1)

ggplot(elegans_training_df, aes(x = factor(genotype), y = TE_specific_Orphan_Read_IQR, fill = factor(genotype))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(x = "Genotype",
       y = "TE Specific Orphan Read (IQR)",
       fill = "Genotype") +
  theme_minimal() +
  ylim(0, 1)

ggplot(elegans_training_df, aes(x = region_length, y = TE_specific_Paired_median, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired Median",
           color = "Insert Size") +
      theme_minimal()

ggplot(elegans_training_df, aes(x = region_length, y = TE_specific_Paired_IQR, color = factor(RL))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE Specific Paired IQR",
           color = "Read Length") +
      theme_minimal()

ggplot(elegans_training_df, aes(x = factor(genotype), y = TE_specific_Paired_IQR, fill = factor(genotype))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(x = "Genotype",
       y = "TE_specific_Paired_IQR",
       fill = "Genotype") +
  theme_minimal() +
  ylim(0, 1)

ggplot(elegans_training_df, aes(x = region_length, y = Proper_Pair_sd, color = factor(genotype))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "Insert_Size_sd",
           color = "Insert Size") +
      theme_minimal()

ggplot(elegans_training_df, aes(x = factor(genotype), y = Proper_Pair_sd, fill = factor(genotype))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(x = "Genotype",
       y = "Proper_Pair_sd",
       fill = "Genotype") +
  theme_minimal()

ggplot(elegans_training_df, aes(x = region_length, y = TE_specific_Insert_Size_mean, color = factor(genotype))) +
      geom_point(alpha = 0.8) +
      labs(x = "Region Length",
           y = "TE_specific_Insert_Size_mean",
           color = "Insert Size") +
      theme_minimal() + ylim(0,100)

ggplot(elegans_training_df, aes(x = factor(genotype), y = Cigar4_sd, fill = factor(genotype))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(x = "Genotype",
       y = "TE_specific_Insert_Size_mean",
       fill = "Genotype") +
  theme_minimal() + ylim(0,100)
