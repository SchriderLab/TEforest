library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
library(viridis)
library(yardstick)
library(UpSetR) 

# load benchmarking data
#load("/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/2L_2R_plots/JUT-008_MUN-009.RData")
#load("/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/2L_2R_plots/A2_A3.RData")

#load("/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/2L_2R_plots/AKA-017_GIM-024.RData")

load("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/benchmark_mapping_results_150_50X.RData")


classification_df <- read.csv("/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/3L3RX_classifer/JUT-008_MUN-009.csv")
# Step 1: Split the strings, making sure "JUT-008_MUN-009" stays intact
split_files <- strsplit(sub("^JUT-008_MUN-009-", "", classification_df$file), "-")

# Step 2: Extract the parts and add the fixed part "JUT-008_MUN-009" to the front
df_expanded <- do.call(rbind, lapply(split_files, function(x) {
  c("JUT-008_MUN-009", x[1], x[2], x[3], x[4], paste(x[5:length(x)], collapse = "-"))
}))

# Step 3: Convert the result into a data frame and rename the columns
df_expanded <- as.data.frame(df_expanded, stringsAsFactors = FALSE)
colnames(df_expanded) <- c("JUT_MUN", "seqnames", "start", "end", "Strand", "TE")

# Step 4: Combine with the original data frame (excluding the old file column)
classification_df_final <- cbind(df_expanded, classification_df[ , -1])
false_negatives_classifier <- classification_df_final %>% filter(true %in% c(1, 2) & pred == 0)
false_negatives_classifier <- false_negatives_classifier #%>%
  #rename(seqnames = Chromosome, start = Start, end = End)
#GenomicRanges::Grange(classification_df_final)

false_negatives_classifier_Pelement <- false_negatives_classifier %>% filter(TE=="P_element")

benchmark_mapping_results_TEforest <- benchmark_mapping_results %>% 
    filter(caller=="TEforest_classifier_bps")


tefam_f1score <- benchmark_mapping_results_TEforest %>%
  select(TE, f1_score, precision, recall, nonref_calls_true_positives_length, nonref_calls_false_positives_length, nonref_calls_false_negatives_length) %>%        # Select the desired columns
  #drop_na %>%                # Remove columns with all NA values
  arrange(desc(nonref_calls_false_negatives_length)) %>% 
  # Filter out rows where all three specified columns are zero
  filter(!(nonref_calls_true_positives_length == 0 & 
           nonref_calls_false_positives_length == 0 & 
           nonref_calls_false_negatives_length == 0)) %>%
  # Replace all remaining NA values with 0
  replace_na(list(
    f1_score = 0,
    precision = 0,
    recall = 0,
    nonref_calls_true_positives_length = 0,
    nonref_calls_false_positives_length = 0,
    nonref_calls_false_negatives_length = 0
  ))

false_positives <- do.call(
    c, c(benchmark_mapping_results_TEforest$nonref_false_positives)
    )

false_negatives <- do.call(
    c, c(benchmark_mapping_results_TEforest$nonref_false_negatives)
    )

true_positives_inrefcoords <- do.call(
    c, c(benchmark_mapping_results_TEforest$nonref_true_positives_truth_dataset_coords)
    )

# Convert each GRanges object to a data frame
df_list <- lapply(false_negatives, function(gr) {
  as.data.frame(gr)
})

# Combine all data frames into one
final_df <- do.call(rbind, df_list)
final_df$start <- final_df$start + 500
final_df$end <- final_df$end - 500

# Convert each GRanges object to a data frame
df_list_tp <- lapply(true_positives_inrefcoords, function(gr) {
  as.data.frame(gr)
})

# Combine all data frames into one
final_df_tp_inrefcoords <- do.call(rbind, df_list_tp)


# Read the CSV, skipping the first line (which is the unwanted header)
te_data <- read.csv("/nas/longleaf/home/adaigle/Rech_TE_info.csv",
    skip = 1, header = TRUE) %>%
    select(Chr,Start,End,Class,Order,SuperFamily,Family,Genomes,
        TE.length,TE.ratio.vs.family.consensus.sequence..in.strains.)
te_data <- te_data %>%
  filter(!grepl("ISO-1", Genomes))
te_data %>% filter(Family=="P_element")
#te_data_test <- te_data[-420,]
te_data$Genomes <- gsub(";\\.;", ";100;", te_data$Genomes)
te_data$TE.length <- gsub(";\\.;", ";100;", te_data$TE.length)
te_data$TE.ratio.vs.family.consensus.sequence..in.strains. <- gsub(";\\.;", ";NA;", te_data$TE.ratio.vs.family.consensus.sequence..in.strains.)
te_data$Family <- gsub("-", "_", te_data$Family)



te_data_long <- te_data %>%
  separate_rows(
    Genomes, TE.length, TE.ratio.vs.family.consensus.sequence..in.strains., 
    sep = ";"
    )

te_data_long_mygenomes <- 
    te_data_long %>% filter(Genomes %in% c("JUT-008", "MUN-009"))
#te_data_long_mygenomes <- 
#    te_data_long %>% filter(Genomes %in% c("JUT-008"))

te_data_long_mygenomes <- te_data_long_mygenomes %>%
  #rename(
  #  seqnames = "Chr",
  #  start = "Start",
  #  end = "End",
  #  TE = "Family" 
  #) 
  rename(
    Chr = "seqnames",
    Start = "start",
    End = "end",
    Family = "TE" 
  ) 

te_data_long_mygenomes <- te_data_long_mygenomes %>% filter(seqnames %in% c("2L", "2R"))

#match my false negatives data to the rech supplemental information
false_negatives <- inner_join(te_data_long_mygenomes, final_df, by = c("seqnames", "start", "end", "TE")) 
true_positives_inrefcoords <- inner_join(te_data_long_mygenomes, final_df_tp_inrefcoords, by = c("seqnames", "start", "end", "TE")) 




#these are assumed to be true positives since they do not show up in the false negatives set
#true_positives <- anti_join(te_data_long_mygenomes, final_df, by = c("seqnames", "start", "end", "TE")) %>%
#  filter(seqnames %in% c("2L", "2R"))

true_positives <- inner_join(te_data_long_mygenomes, final_df_tp_inrefcoords, by = c("seqnames", "start", "end", "TE")) 

false_negatives$Genomes <- NULL
false_negatives <- false_negatives %>% distinct()
false_negatives$TE.length <- as.numeric(false_negatives$TE.length)

gr_false_negatives <- GRanges(
  seqnames = false_negatives$seqnames,
  ranges = IRanges(start = false_negatives$start, end = false_negatives$end),
  TE = false_negatives$TE
)

# Create GRanges object for false_negatives_classifier
gr_false_negatives_classifier <- GRanges(
  seqnames = false_negatives_classifier$seqnames,
  ranges = IRanges(start = as.numeric(false_negatives_classifier$start), end = as.numeric(false_negatives_classifier$end)),
  TE = false_negatives_classifier$TE
)
overlaps <- findOverlaps(gr_false_negatives, gr_false_negatives_classifier)
# Convert overlaps to a data frame for easier manipulation
overlap_df <- as.data.frame(overlaps)

# Add the TE information to the overlap data frame
overlap_df$TE_query <- mcols(gr_false_negatives)$TE[overlap_df$queryHits]
overlap_df$TE_subject <- mcols(gr_false_negatives_classifier)$TE[overlap_df$subjectHits]

# Filter overlaps where TE values are equal
overlap_df_TE_equal <- overlap_df[overlap_df$TE_query == overlap_df$TE_subject, ]

overlaps_with_equal_TE <- rep(FALSE, length(gr_false_negatives))

# Mark TRUE for queries that have overlaps with equal TE
overlaps_with_equal_TE[unique(overlap_df_TE_equal$queryHits)] <- TRUE

# Add the overlap information to the false_negatives data frame
false_negatives$overlaps_with_classifier_TE_equal <- overlaps_with_equal_TE



false_negatives_LTR <- false_negatives %>% filter(Order=="LTR")

df_summary <- false_negatives_LTR %>%
  group_by(TE) %>%
  summarize(
    mean_length = mean(TE.length, na.rm = TRUE),
    sd_length = sd(TE.length, na.rm = TRUE),
    n = n(),
    se_length = sd_length / sqrt(n)
  )

# Plot lengths of different TE families
te_length_plot <- ggplot(df_summary, aes(x = TE, y = mean_length, fill=TE)) +
  # Bar plot for the mean values
  geom_bar(stat = "identity",width = 0.7) +
  # Error bars for the standard error
  geom_errorbar(aes(ymin = mean_length - se_length, ymax = mean_length + se_length), 
                width = 0.2) +
  # Add jittered points for individual observations from the original df
  geom_jitter(data = false_negatives_LTR, aes(x = TE, y = TE.length), width = 0.2, size = 2, alpha = 0.6) +
  # Add the mean point on top of the bars
  geom_point(aes(y = mean_length), size = 3) +
  labs(title = "Length Grouped by TE Class with Error Bars and Individual Points", 
       x = "TE Class", y = "Length") +
  
  theme_minimal()


# Make sure TE.length is numeric in both dataframes
false_negatives$TE.length <- as.numeric(false_negatives$TE.length)
true_positives$TE.length <- as.numeric(true_positives$TE.length)
true_positives$TE.ratio.vs.family.consensus.sequence..in.strains. <- as.numeric(true_positives$TE.ratio.vs.family.consensus.sequence..in.strains.)

# Step 1: Combine the two data frames into one, adding a column to indicate their origin
false_negatives$Source <- "False negatives"
true_positives$Source <- "True positives"


# Make sure TE.ratio.vs.family.consensus.sequence..in.strains. is numeric in both dataframes
false_negatives$TE.ratio.vs.family.consensus.sequence..in.strains. <- as.numeric(false_negatives$TE.ratio.vs.family.consensus.sequence..in.strains.)
true_positives$TE.ratio.vs.family.consensus.sequence..in.strains. <- as.numeric(true_positives$TE.ratio.vs.family.consensus.sequence..in.strains.)

combined_df <- bind_rows(false_negatives, true_positives)


# Step 1: Combine the two data frames into one, adding a column to indicate their origin
false_negatives$Source <- "False negatives"
true_positives$Source <- "True positives"

plot_TE_comparison <- function(TE_list, false_negatives, true_positives) {
  
    # Step 1: Filter both dataframes to only include rows where TE matches the current TE_name
    false_negatives_filtered <- false_negatives %>% filter(TE %in% TE_list)
    true_positives_filtered <- true_positives %>% filter(TE %in% TE_list)
    print(nrow(false_negatives_filtered))
    # Step 2: Combine the filtered data frames into one, adding a column to indicate their origin
    false_negatives_filtered$Source <- "False negatives"
    true_positives_filtered$Source <- "True positives"
    
    combined_df <- bind_rows(false_negatives_filtered, true_positives_filtered)
    
    # Check if there is sufficient data to perform the test and generate plots
    if (nrow(combined_df) > 0) {
      
      # Step 3: Perform Wilcoxon rank-sum test between the two groups
      wilcox_test <- wilcox.test(
        TE.length ~ Source, 
        data = combined_df
      )
      
      # Extract the p-value from the test
      p_value <- wilcox_test$p.value
      print(p_value)
      # Step 4: Create the plot comparing the TE.ratio.vs.family.consensus.sequence..in.strains. 
      # between "False negatives" and "True positives", and add the p-value to the plot title
      plot <- ggplot(combined_df, aes(x = Source, y = TE.length, color = Source)) +
        # Add jittered points for individual observations from both data sets
        geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
        
        # Add boxplot to show the distribution of TE.ratio.vs.family.consensus.sequence..in.strains. in both sets
        geom_boxplot(outlier.shape = NA, alpha = 0.2, width = 0.4) +
        
        # Add p-value to the plot title
        labs(title = paste("Comparison of TEs Between False Negatives and True Positives", 
                           "\nWilcoxon p-value:", round(p_value, 5)), 
             x = "Source", 
             y = "TE Length vs Family Consensus Sequence (in strains)") +
        
        theme_minimal()
      
      # Print the plot for each TE in the list
      print(plot)
      
    } else {
      # In case there are no rows for a given TE, print a message
      message(paste("No data available for TE:", TE_name))
    }
}

plot_TE_comparison_ratio <- function(TE_list, false_negatives, true_positives) {
  
    # Step 1: Filter both dataframes to only include rows where TE matches the current TE_name
    false_negatives_filtered <- false_negatives %>% filter(TE %in% TE_list)
    true_positives_filtered <- true_positives %>% filter(TE %in% TE_list)
    print(nrow(false_negatives_filtered))
    # Step 2: Combine the filtered data frames into one, adding a column to indicate their origin
    false_negatives_filtered$Source <- "False negatives"
    true_positives_filtered$Source <- "True positives"
    
    combined_df <- bind_rows(false_negatives_filtered, true_positives_filtered)
    
    # Check if there is sufficient data to perform the test and generate plots
    if (nrow(combined_df) > 0) {
      
      # Step 3: Perform Wilcoxon rank-sum test between the two groups
      wilcox_test <- wilcox.test(
        TE.ratio.vs.family.consensus.sequence..in.strains. ~ Source, 
        data = combined_df
      )
      
      # Extract the p-value from the test
      p_value <- wilcox_test$p.value
      print(p_value)
      # Step 4: Create the plot comparing the TE.ratio.vs.family.consensus.sequence..in.strains. 
      # between "False negatives" and "True positives", and add the p-value to the plot title
      plot <- ggplot(combined_df, aes(x = Source, y = TE.ratio.vs.family.consensus.sequence..in.strains., color = Source)) +
        # Add jittered points for individual observations from both data sets
        geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
        
        # Add boxplot to show the distribution of TE.ratio.vs.family.consensus.sequence..in.strains. in both sets
        geom_boxplot(outlier.shape = NA, alpha = 0.2, width = 0.4) +
        
        # Add p-value to the plot title
        labs(title = paste("Comparison of TEs Between False Negatives and True Positives", 
                           "\nWilcoxon p-value:", round(p_value, 5)), 
             x = "Source", 
             y = "TE Length vs Family Consensus Sequence (in strains)") +
        
        theme_minimal()
      
      # Print the plot for each TE in the list
      print(plot)
      
    } else {
      # In case there are no rows for a given TE, print a message
      message(paste("No data available for TE:", TE_name))
    }
}
#comparing some random te lengths
#plot_TE_comparison("roo", false_negatives, true_positives)
#plot_TE_comparison("Copia", false_negatives, true_positives)
#plot_TE_comparison("jockey", false_negatives, true_positives)
#plot_TE_comparison("FB4", false_negatives, true_positives)
#plot_TE_comparison("Micropia", false_negatives, true_positives)
#plot_TE_comparison("pogo", false_negatives, true_positives)
#plot_TE_comparison("INE_1", false_negatives, true_positives)
#plot_TE_comparison(c("roo","Copia","pogo"), false_negatives, true_positives)

# compare all tes 
alltes <- unique(false_negatives$TE, true_positives$TE)
plot_TE_comparison(alltes, false_negatives, true_positives)
plot_TE_comparison_ratio(alltes, false_negatives, true_positives)

# compare nested/nonnested TEs--nested fn are longer
false_negatives_nested <- false_negatives %>% filter(overlap_category!="no_overlap")
false_negatives_nonnested <- false_negatives %>% filter(overlap_category=="no_overlap")

plot_TE_comparison(alltes, false_negatives_nested, true_positives)
plot_TE_comparison(alltes, false_negatives_nonnested, true_positives)
plot_TE_comparison_ratio(alltes, false_negatives_nested, true_positives)

plot_TE_comparison_ratio(alltes, false_negatives_nonnested, true_positives)
plot_TE_comparison_ratio(alltes, false_negatives_nonnested, false_negatives_nested)

# compare fn from candidate region stage to classifier stage
false_negatives_classifier <- false_negatives %>% filter(overlaps_with_classifier_TE_equal==TRUE)
false_negatives_candidateregions <- false_negatives %>% filter(overlaps_with_classifier_TE_equal==FALSE)


plot_TE_comparison(alltes, false_negatives_classifier, false_negatives_candidateregions)

plot_TE_comparison(alltes, false_negatives_classifier, true_positives)
plot_TE_comparison(alltes, false_negatives_candidateregions, true_positives)
plot_TE_comparison_ratio(alltes, false_negatives_classifier, true_positives)
#plot_TE_comparison(alltes, false_negatives_classifier_nested, true_positives)

# side point--nested vs nonnested in classifier not sig different
false_negatives_classifier_nested <- false_negatives_classifier %>% filter(overlap_category!="no_overlap")
false_negatives_classifier_nonnested <- false_negatives_classifier %>% filter(overlap_category=="no_overlap")


plot_TE_comparison_ratio(alltes, false_negatives_classifier_nested, false_negatives_classifier_nonnested)
plot_TE_comparison(alltes, false_negatives_classifier_nested, false_negatives_classifier_nonnested)
plot_TE_comparison_ratio(alltes, false_negatives_classifier_nonnested, true_positives)

#now put both together -- remove nested tes and only look at candidate regions
false_negatives_nonested_candidateregions <- false_negatives %>% filter(overlaps_with_classifier_TE_equal==FALSE) %>% filter(overlap_category=="no_overlap")
false_negatives_nested_candidateregions <- false_negatives %>% filter(overlaps_with_classifier_TE_equal==FALSE) %>% filter(overlap_category!="no_overlap")

plot_TE_comparison_ratio(alltes, false_negatives_classifier_nested, false_negatives_nested_candidateregions)

plot_TE_comparison(alltes, false_negatives_nonested_candidateregions, true_positives)
plot_TE_comparison_ratio(alltes, false_negatives_nonested_candidateregions, true_positives)

plot_TE_comparison(alltes, false_negatives_nested_candidateregions, false_negatives_nonested_candidateregions)
plot_TE_comparison_ratio(alltes, false_negatives_nested_candidateregions, false_negatives_nonested_candidateregions)

high_ratio_falsenegatives <- false_negatives %>% filter(TE.ratio.vs.family.consensus.sequence..in.strains. > .8)

#do not double count heterozygotes or nested TEs of the same type here.
true_positives_unique <- true_positives %>% select(seqnames, start, end, TE,overlap_category)
true_positives_unique <- unique(true_positives_unique)

# Step 1: Get counts of TEs in each dataset
counts_merged <- table(false_negatives$TE)
counts_non_matching <- table(true_positives$TE)
counts_candidate_region_fn <- table(false_negatives_candidateregions$TE)

# Convert these tables to data frames
df_merged <- as.data.frame(counts_merged)
df_non_matching <- as.data.frame(counts_non_matching)
df_counts_candidate_region_fn <- as.data.frame(counts_candidate_region_fn)

# Rename columns for clarity
colnames(df_merged) <- c("TE", "Total_False_Negatives")
colnames(df_non_matching) <- c("TE", "Total_True_positives")
colnames(df_counts_candidate_region_fn) <- c("TE", "Total_Candidate_Region_False_Negatives")

# Step 2: Calculate proportions for each TE
total_merged <- sum(df_merged$Total_False_Negatives)
total_non_matching <- sum(df_non_matching$Total_True_positives)
total_candidate_region_fn <- sum(df_non_matching$Total_Candidate_Region_False_Negatives)

df_merged <- df_merged %>% mutate(Proportion_False_Negatives = Total_False_Negatives / total_merged)
df_non_matching <- df_non_matching %>% mutate(Proportion_True_Positives = Total_True_positives / total_non_matching)
df_counts_candidate_region_fn <- df_counts_candidate_region_fn %>% mutate(Proportion_Candidate_Region_False_Negatives = Total_Candidate_Region_False_Negatives / total_non_matching)

# Step 3: Merge the two data frames by TE, so we can plot them side by side
combined_df <- full_join(df_merged, df_non_matching, by = "TE")
combined_df <- full_join(combined_df, df_counts_candidate_region_fn, by = "TE")

# Replace NAs with 0 for cases where a TE is present in one dataset but not the other
combined_df[is.na(combined_df)] <- 0

# Step 4: Calculate SDs for each TE between the datasets
combined_df <- combined_df %>%
  mutate(SD = sqrt((Proportion_False_Negatives^2 + Proportion_True_Positives^2) / 2))

TEfam_performance_df <- combined_df %>% arrange(-Total_False_Negatives)
TEfam_performance_df$SD <- NULL

write.csv(TEfam_performance_df, "/nas/longleaf/home/adaigle/TEforest/plots/TEfam_perf_df.csv", row.names = FALSE)


# Step 5: Reshape data for plotting (gather proportions into long format for ggplot)
combined_long <- combined_df %>%
  select(TE, Proportion_False_Negatives, Proportion_True_Positives) %>%
  gather(key = "Dataset", value = "Proportion", -TE)

# Step 6: Plot the proportions side by side in a bar plot with SDs
ggplot(combined_long, aes(x = TE, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Proportion of TEs in Merged and Non-Matching Rows",
       x = "TE", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_TE_comparison_four <- function(TE_list, true_positives, candidate_false_negatives, classifier_false_negatives, nested_false_negatives) {
  
  # Step 1: Filter each dataframe to include only rows where TE matches the current TE_list
  true_positives_filtered <- true_positives %>% filter(TE %in% TE_list)
  candidate_fn_filtered <- candidate_false_negatives %>% filter(TE %in% TE_list)
  classifier_fn_filtered <- classifier_false_negatives %>% filter(TE %in% TE_list)
  nested_fn_filtered <- nested_false_negatives %>% filter(TE %in% TE_list)
  
  # Add a source label to each dataframe
  true_positives_filtered$Source <- "True Positives"
  candidate_fn_filtered$Source <- "Candidate FN"
  classifier_fn_filtered$Source <- "Classifier FN"
  nested_fn_filtered$Source <- "Nested FN"
  
  # Combine all dataframes into one for plotting
  combined_df <- bind_rows(true_positives_filtered, candidate_fn_filtered, classifier_fn_filtered, nested_fn_filtered)
  
  # Check if there is sufficient data to perform the test and generate plots
  if (nrow(combined_df) > 0) {
    
    # Step 2: Perform Wilcoxon tests between each "false negatives" group and "true positives"
    wilcox_candidate <- wilcox.test(candidate_fn_filtered$TE.length, true_positives_filtered$TE.length)$p.value
    wilcox_classifier <- wilcox.test(classifier_fn_filtered$TE.length, true_positives_filtered$TE.length)$p.value
    wilcox_nested <- wilcox.test(nested_fn_filtered$TE.length, true_positives_filtered$TE.length)$p.value
    
    # Step 3: Create significance labels based on p-values
    significance_labels <- function(p) {
      if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else if (p < 0.05) return("*")
      else return("ns")
    }
    
    candidate_label <- significance_labels(wilcox_candidate)
    classifier_label <- significance_labels(wilcox_classifier)
    nested_label <- significance_labels(wilcox_nested)
    
    # Step 4: Create the plot with significance bars
    plot <- ggplot(combined_df, aes(x = Source, y = TE.length)) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
      geom_boxplot(outlier.shape = NA, alpha = 0.2, width = 0.4) +
      labs(
        x = "Classification",
        y = "TE length (bp)"
      ) + 
      theme_minimal() + theme(legend.title=element_blank())+
      
      # Step 5: Add significance bars with labels
      stat_pvalue_manual(
        data = data.frame(
          group1 = c("True Positives", "True Positives", "True Positives"),
          group2 = c("Candidate FN", "Classifier FN", "Nested FN"),
          p.adj = c(wilcox_candidate, wilcox_classifier, wilcox_nested),
          y.position = c(
            max(combined_df$TE.length, na.rm = TRUE) * 1.1,
            max(combined_df$TE.length, na.rm = TRUE) * 1.2,
            max(combined_df$TE.length, na.rm = TRUE) * 1.3
          ),
          label = c(candidate_label, classifier_label, nested_label)
        ),
        label = "label"
      )
    
    # Print the plot
    return(plot)
    
  } else {
    message("No data available for the specified TEs.")
  }
}

plot_TE_comparison_ratio_four <- function(TE_list, true_positives, candidate_false_negatives, classifier_false_negatives, nested_false_negatives) {
  
  # Step 1: Filter each dataframe to include only rows where TE matches the current TE_list
  true_positives_filtered <- true_positives %>% filter(TE %in% TE_list)
  candidate_fn_filtered <- candidate_false_negatives %>% filter(TE %in% TE_list)
  classifier_fn_filtered <- classifier_false_negatives %>% filter(TE %in% TE_list)
  nested_fn_filtered <- nested_false_negatives %>% filter(TE %in% TE_list)
  
  # Add a source label to each dataframe
  true_positives_filtered$Source <- "True Positives"
  candidate_fn_filtered$Source <- "Candidate FN"
  classifier_fn_filtered$Source <- "Classifier FN"
  nested_fn_filtered$Source <- "Nested FN"
  
  # Combine all dataframes into one for plotting
  combined_df <- bind_rows(true_positives_filtered, candidate_fn_filtered, classifier_fn_filtered, nested_fn_filtered)
  
  # Check if there is sufficient data to perform the test and generate plots
  if (nrow(combined_df) > 0) {
    
    # Step 2: Perform Wilcoxon tests between each "false negatives" group and "true positives"
    wilcox_candidate <- wilcox.test(candidate_fn_filtered$TE.ratio.vs.family.consensus.sequence..in.strains., true_positives_filtered$TE.ratio.vs.family.consensus.sequence..in.strains.)$p.value
    wilcox_classifier <- wilcox.test(classifier_fn_filtered$TE.ratio.vs.family.consensus.sequence..in.strains., true_positives_filtered$TE.ratio.vs.family.consensus.sequence..in.strains.)$p.value
    wilcox_nested <- wilcox.test(nested_fn_filtered$TE.ratio.vs.family.consensus.sequence..in.strains., true_positives_filtered$TE.ratio.vs.family.consensus.sequence..in.strains.)$p.value
    
    # Step 3: Create significance labels based on p-values
    significance_labels <- function(p) {
      if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else if (p < 0.05) return("*")
      else return("ns")
    }
    
    candidate_label <- significance_labels(wilcox_candidate)
    classifier_label <- significance_labels(wilcox_classifier)
    nested_label <- significance_labels(wilcox_nested)
    
    # Step 4: Create the plot with significance bars
    plot <- ggplot(combined_df, aes(x = Source, y = TE.ratio.vs.family.consensus.sequence..in.strains.)) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
      geom_boxplot(outlier.shape = NA, alpha = 0.2, width = 0.4) +
      labs(
        x = "Classification",
        y = "TE length relative to family consensus sequence length"
      ) +
      theme_minimal() + theme(legend.title=element_blank()) +
      
      # Step 5: Add significance bars with labels
      stat_pvalue_manual(
        data = data.frame(
          group1 = c("True Positives", "True Positives", "True Positives"),
          group2 = c("Candidate FN", "Classifier FN", "Nested FN"),
          p.adj = c(wilcox_candidate, wilcox_classifier, wilcox_nested),
          y.position = c(
            max(combined_df$TE.ratio.vs.family.consensus.sequence..in.strains., na.rm = TRUE) * 1.1,
            max(combined_df$TE.ratio.vs.family.consensus.sequence..in.strains., na.rm = TRUE) * 1.2,
            max(combined_df$TE.ratio.vs.family.consensus.sequence..in.strains., na.rm = TRUE) * 1.3
          ),
          label = c(candidate_label, classifier_label, nested_label)
        ),
        label = "label"
      )
    
    # Print the plot
    return(plot)
    
  } else {
    message("No data available for the specified TEs.")
  }
}

TElength_comparision <- plot_TE_comparison_four(alltes, true_positives, false_negatives_nonested_candidateregions, false_negatives_classifier_nonnested, false_negatives_nested)
TEratio_comparision <- plot_TE_comparison_ratio_four(alltes, true_positives, false_negatives_nonested_candidateregions, false_negatives_classifier_nonnested, false_negatives_nested)


### analyze unique true positives

benchmark_mapping_results_TEforest <- benchmark_mapping_results %>% 
        filter(caller=="TEforest_classifier_bps")

benchmark_mapping_results_not_TEforest <- benchmark_mapping_results %>% 
        filter(caller %in% c("popoolationte", "retroseq", "teflon", "temp", "temp2"))

benchmark_mapping_results_retroseq <- benchmark_mapping_results %>% 
        filter(caller %in% c("retroseq"))

benchmark_mapping_results_all <- benchmark_mapping_results %>% 
        filter(caller %in% c("TEforest_classifier_bps", "popoolationte", "retroseq", "teflon", "temp", "temp2"))

benchmark_mapping_results_all <- benchmark_mapping_results_all %>%
  mutate(nonref_true_positives_truth_dataset_coords = map2(nonref_true_positives_truth_dataset_coords, caller, ~ {
    # If the GRanges object is not empty, assign the caller value
    if (length(.x) > 0) {
      mcols(.x)$caller <- .y
    } else {
      # For empty GRanges, we can leave it as is or explicitly set an empty column
      # mcols(.x)$caller <- character(0) # optional if you want an empty column present
    }
    .x
  }))

false_positives_all <- do.call(
    c, c(benchmark_mapping_results_all$nonref_false_positives)
    )

# Step 1: Group GRanges by TE family names
# Extract unique TE names
unique_te_names <- unique(names(false_positives_all))

# Step 1: Combine GRanges by TE family name
combine_granges_by_name <- function(gr_list) {
  # Get unique TE family names
  unique_te_names <- unique(names(gr_list))
  
  # Combine all GRanges with the same TE name
  combined <- lapply(unique_te_names, function(te_name) {
    granges_with_name <- gr_list[names(gr_list) == te_name]
    do.call(c, granges_with_name) # Concatenate all GRanges with the same name
  })
  
  # Name the list based on TE families
  #names(combined) <- unique_te_names
  combined
}

# Step 2: Flatten nested lists
flatten_combined_granges <- function(combined_granges) {
  lapply(combined_granges, function(gr_list) {
    do.call(c, unname(unlist(gr_list))) # Flatten and unname
  })
}

combined_granges <- combine_granges_by_name(false_positives_all)
flattened_granges <- flatten_combined_granges(combined_granges)

# Step 2: Merge overlapping ranges for each TE family
merged_granges <- lapply(flattened_granges, GenomicRanges::reduce)
names(merged_granges) <- unique_te_names

benchmark_mapping_results_all <- benchmark_mapping_results_all %>%
  mutate(
    nonref_falsepositives_unifiedcoords = map2(
      nonref_false_positives, TE,
      ~ subsetByOverlaps(merged_granges[[.y]], .x)
    )
  )

benchmark_mapping_results_all <- benchmark_mapping_results_all %>%
  mutate(nonref_falsepositives_unifiedcoords = map2(nonref_falsepositives_unifiedcoords, caller, ~ {
    # If the GRanges object is not empty, assign the caller value
    if (length(.x) > 0) {
      mcols(.x)$caller <- .y
    } else {
      # For empty GRanges, we can leave it as is or explicitly set an empty column
      # mcols(.x)$caller <- character(0) # optional if you want an empty column present
    }
    .x
  }))

granges_list_to_df <- function(df_column) {
calls <- do.call(c, c(df_column))

gr_list_df <- lapply(calls, function(gr) {
  df <- as.data.frame(gr)
  return(df)
})

# Combine all data frames into a single data frame
final_df <- bind_rows(gr_list_df, .id = "source")  # Adds a column for the list element ID
return(final_df)
}

not_TEforest_df <- granges_list_to_df(benchmark_mapping_results_not_TEforest$nonref_true_positives_truth_dataset_coords)
allcallers_tp_df <- granges_list_to_df(benchmark_mapping_results_all$nonref_true_positives_truth_dataset_coords)
allcallers_fp_df <- granges_list_to_df(benchmark_mapping_results_all$nonref_falsepositives_unifiedcoords)



#TEforest_df <- granges_list_to_df(benchmark_mapping_results_TEforest$nonref_true_positives_truth_dataset_coords)
#
#retroseq_df <- granges_list_to_df(benchmark_mapping_results_retroseq$nonref_true_positives_truth_dataset_coords)
#
#unique_tp_TEforest <- anti_join(TEforest_df, not_TEforest_df)
#
#unique_tp_notTEforest <- anti_join(not_TEforest_df, TEforest_df)
#
#unique_tp_TEforest_vsRetroseq <- anti_join(TEforest_df, retroseq_df)
#
#unique_tp_TEforest_vsRetroseq_inverse <- anti_join(retroseq_df, TEforest_df)
#
#
#unique_tp_TEforest_fulldata <- inner_join(te_data_long_mygenomes, unique_tp_TEforest, by = c("seqnames", "start", "end", "TE")) 
#unique_tp_TEforest_fulldata$TE.length <- as.numeric(unique_tp_TEforest_fulldata$TE.length)
#unique_tp_TEforest_fulldata$TE.ratio.vs.family.consensus.sequence..in.strains. <- as.numeric(unique_tp_TEforest_fulldata$TE.ratio.vs.family.consensus.sequence..in.strains.)
#unique_tp_TEforest_fulldata$Source <- "Unique true positives"
#
#unique_tp_TEforest_vsRetroseq_fulldata <- unique(inner_join(te_data_long_mygenomes, unique_tp_TEforest_vsRetroseq, by = c("seqnames", "start", "end", "TE")))
#unique_tp_TEforest_vsRetroseq_fulldata$TE.length <- as.numeric(unique_tp_TEforest_vsRetroseq_fulldata$TE.length)
#unique_tp_TEforest_vsRetroseq_fulldata$TE.ratio.vs.family.consensus.sequence..in.strains. <- as.numeric(unique_tp_TEforest_vsRetroseq_fulldata$TE.ratio.vs.family.consensus.sequence..in.strains.)
#unique_tp_TEforest_vsRetroseq_fulldata$Source <- "Unique true positives"
#
#plot_TE_comparison(alltes, unique_tp_TEforest_fulldata, true_positives)
#plot_TE_comparison_ratio(alltes, unique_tp_TEforest_fulldata, true_positives)
#plot_TE_comparison(alltes, unique_tp_TEforest_vsRetroseq_fulldata, true_positives)
#plot_TE_comparison_ratio(alltes, unique_tp_TEforest_vsRetroseq_fulldata, true_positives)

### I am upset
allcallers_tp_df <- allcallers_tp_df %>%
        mutate(
            caller = case_when(
                caller == "TEforest_classifier_bps" ~ "TEforest",
                caller == "temp" ~ "TEMP",
                caller == "temp2" ~ "TEMP2",
                caller == "teflon" ~ "TEFLoN",
                caller == "retroseq" ~ "RetroSeq",
                caller == "popoolationte" ~ "PopoolationTE",
                caller == "popoolationte2" ~ "PopoolationTE2",
                TRUE ~ caller # Keep the original value for other cases
            )
        )

# Identify all callers present in the data
unique_callers <- unique(allcallers_tp_df$caller)


# For each unique combination of non-caller columns, gather the callers
allcallers_tp_df$ID <- NULL

# Columns other than 'caller' that define uniqueness
non_caller_cols <- setdiff(names(allcallers_tp_df), "caller")

df_for_upset <- unique(allcallers_tp_df) %>% #Do not double count nested tes of the same type in tp
  group_by(across(all_of(non_caller_cols))) %>%
  summarize(callers = list(unique(caller)), .groups = "drop")

#I need to filter out entries within 500 basepairs of each other, as I would be double counting these otherwise
nearby_entries <- df_for_upset %>%
  # Group by seqnames and TE to ensure comparisons are only within these groups
  group_by(seqnames, TE) %>%
  # Only proceed if there's more than one row per group
  filter(n() > 1) %>%
  # Self-join the group to compare each row with every other row in the same group
  inner_join(df_for_upset, by = c("seqnames", "TE"), suffix = c("", ".other")) %>%
  # Exclude identical pairs (where start and end match exactly the same row)
  filter(!(start == start.other & end == end.other)) %>%
  # Now apply the ±500 bp condition
  filter(
    # Check if other start or end falls into the current row's expanded range
    (start.other >= (start - 500) & start.other <= (end + 500)) |
    (end.other   >= (start - 500) & end.other   <= (end + 500)) |
    # Also check if the current row falls into the other's expanded range
    (start       >= (start.other - 500) & start       <= (end.other + 500)) |
    (end         >= (start.other - 500) & end         <= (end.other + 500))
  ) %>%
  # Keep only unique rows from the original set
  distinct(source,seqnames, TE, start, end,width,strand,length,TE,heterozygosity, overlap_category, callers) %>%
  ungroup()

#everything is in pairs, so just keeping one of the pairs
odd_rows_nearby_entries <- nearby_entries[seq(1, nrow(nearby_entries), by = 2), ]
df_for_upset <- anti_join(df_for_upset,odd_rows_nearby_entries)

presence_absence_mat <- sapply(unique_callers, function(cl) {
  sapply(df_for_upset$callers, function(x) cl %in% x)
})

# Convert the matrix to a data frame
presence_absence_df <- as.data.frame(presence_absence_mat)
colnames(presence_absence_df) <- unique_callers

# Bind the new presence/absence columns onto df_for_upset
df_for_upset <- bind_cols(df_for_upset, presence_absence_df)

df_for_upset_sets_only_tp <- df_for_upset[, unique_callers]

tp_df <- as.data.frame(lapply(df_for_upset_sets_only_tp, as.integer))

allcallers_fp_df <- allcallers_fp_df %>%
        mutate(
            caller = case_when(
                caller == "TEforest_classifier_bps" ~ "TEforest",
                caller == "temp" ~ "TEMP",
                caller == "temp2" ~ "TEMP2",
                caller == "teflon" ~ "TEFLoN",
                caller == "retroseq" ~ "RetroSeq",
                caller == "popoolationte" ~ "PopoolationTE",
                caller == "popoolationte2" ~ "PopoolationTE2",
                TRUE ~ caller # Keep the original value for other cases
            )
        )

# Identify all callers present in the data
unique_callers <- unique(allcallers_fp_df$caller)

# Columns other than 'caller' that define uniqueness
non_caller_cols <- setdiff(names(allcallers_fp_df), "caller")

# For each unique combination of non-caller columns, gather the callers
df_for_upset <- allcallers_fp_df %>%
  group_by(across(all_of(non_caller_cols))) %>%
  summarize(callers = list(unique(caller)), .groups = "drop")

presence_absence_mat <- sapply(unique_callers, function(cl) {
  sapply(df_for_upset$callers, function(x) cl %in% x)
})

# Convert the matrix to a data frame
presence_absence_df <- as.data.frame(presence_absence_mat)
colnames(presence_absence_df) <- unique_callers

# Bind the new presence/absence columns onto df_for_upset
df_for_upset <- bind_cols(df_for_upset, presence_absence_df)

df_for_upset_sets_only_fp <- df_for_upset[, unique_callers]

fp_df <- as.data.frame(lapply(df_for_upset_sets_only_fp, as.integer))



# get false negatives
full_te_list <- unique(te_data$Family)

caller_te_list <- benchmark_mapping_results_all %>%
  group_by(caller) %>%
  summarize(te_list = list(unique(TE)), .groups = "drop")

# 2) Find missing TEs for each caller
caller_te_list <- caller_te_list %>%
  mutate(missing_tes = map(te_list, ~ setdiff(full_te_list, .x)))


truth_dataset <- te_data_long_mygenomes
truth_dataset <- unique(truth_dataset)
truth_dataset$Class <- NULL
truth_dataset$Order <- NULL
truth_dataset$SuperFamily <- NULL
truth_dataset$TE.length <- NULL
truth_dataset$TE.ratio.vs.family.consensus.sequence..in.strains. <- NULL

benchmark_mapping_results_all <- benchmark_mapping_results_all %>%
  mutate(nonref_false_negatives = map2(nonref_false_negatives, caller, ~ {
    # If the GRanges object is not empty, assign the caller value
    if (length(.x) > 0) {
      mcols(.x)$caller <- .y
    } else {
      # For empty GRanges, we can leave it as is or explicitly set an empty column
      # mcols(.x)$caller <- character(0) # optional if you want an empty column present
    }
    .x
  }))
allcallers_fn_df <- granges_list_to_df(benchmark_mapping_results_all$nonref_false_negatives)

missing_entries <- caller_te_list %>%
  # Each row now contains a caller and a vector of missing TEs. Unnest them:
  unnest(missing_tes) 

names(missing_entries)[names(missing_entries) == 'missing_tes'] <- 'TE'

missing_entries_df <- missing_entries %>%
  mutate(missing_te_coords = map(TE, ~ {truth_dataset[truth_dataset$TE==.x,]}))

missing_entries_df <- missing_entries_df %>%
  mutate(missing_te_coords = map2(missing_te_coords, caller, ~ {
    # If the GRanges object is not empty, assign the caller value
    if (nrow(.x) > 0) {
      .x$caller <- .y
    } else {
      # For empty GRanges, we can leave it as is or explicitly set an empty column
      # mcols(.x)$caller <- character(0) # optional if you want an empty column present
    }
    .x
  }))

missing_fn <- do.call(rbind, c(missing_entries_df$missing_te_coords))

allcallers_fn_df <- allcallers_fn_df %>% select(seqnames, start, end, TE, caller)

all_fn <- rbind(missing_fn, allcallers_fn_df)

all_fn <- all_fn %>%
        mutate(
            caller = case_when(
                caller == "TEforest_classifier_bps" ~ "TEforest",
                caller == "temp" ~ "TEMP",
                caller == "temp2" ~ "TEMP2",
                caller == "teflon" ~ "TEFLoN",
                caller == "retroseq" ~ "RetroSeq",
                caller == "popoolationte" ~ "PopoolationTE",
                caller == "popoolationte2" ~ "PopoolationTE2",
                TRUE ~ caller # Keep the original value for other cases
            )
        )

# Identify all callers present in the data
unique_callers <- unique(all_fn$caller)

# Columns other than 'caller' that define uniqueness
non_caller_cols <- setdiff(names(all_fn), "caller")

# For each unique combination of non-caller columns, gather the callers
df_for_upset <- all_fn %>%
  group_by(across(all_of(non_caller_cols))) %>%
  summarize(callers = list(unique(caller)), .groups = "drop")

presence_absence_mat <- sapply(unique_callers, function(cl) {
  sapply(df_for_upset$callers, function(x) cl %in% x)
})

# Convert the matrix to a data frame
presence_absence_df <- as.data.frame(presence_absence_mat)
colnames(presence_absence_df) <- unique_callers

# Bind the new presence/absence columns onto df_for_upset
df_for_upset <- bind_cols(df_for_upset, presence_absence_df)

df_for_upset_sets_only <- df_for_upset[, unique_callers]

fn_df <- as.data.frame(lapply(df_for_upset_sets_only, as.integer))

callers <- c("TEforest", "PopoolationTE", "RetroSeq", "TEFLoN", "TEMP", "TEMP2")

# Function to create a combination ID
get_combination_id <- function(df) {
  apply(df[, callers], 1, function(row) {
    active_callers <- callers[row == 1]
    if (length(active_callers) == 0) {
      return("None")
    } else {
      return(paste(sort(active_callers), collapse = "+"))
    }
  })
}
# Generate all combinations of callers (from single sets to all)
all_combinations <- unlist(lapply(1:length(callers), function(k) {
  apply(combn(callers, k), 2, list)
}), recursive = FALSE)

# Updated total true positives
totaltp <- 401

# Calculate performance metrics for each combination of callers
performance_list <- lapply(all_combinations, function(subset_callers) {
  subset_callers <- unlist(subset_callers)
  
  # Calculate total true positives (TP) and false positives (FP) for the given subset
  total_TP <- sum(rowSums(tp_df[, subset_callers, drop = FALSE]) > 0)
  total_FP <- sum(rowSums(fp_df[, subset_callers, drop = FALSE]) > 0)
  
  # Calculate total false negatives (FN) based on total_tp
  total_FN <- totaltp - total_TP
  
  # Compute precision, recall, and F1 score
  precision <- if ((total_TP + total_FP) > 0) total_TP / (total_TP + total_FP) else NA
  recall <- if ((total_TP + total_FN) > 0) total_TP / (total_TP + total_FN) else NA
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall > 0)) {
    2 * precision * recall / (precision + recall)
  } else {
    NA
  }

  # Return a data frame with the metrics
  data.frame(
    combination = paste(subset_callers, collapse = "+"),
    TP = total_TP,
    FP = total_FP,
    FN = total_FN,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
})

perf_df <- do.call(rbind, performance_list)

expand_combination <- function(combination, callers) {
  out <- setNames(as.data.frame(matrix(0, ncol = length(callers), nrow = 1)), callers)
  if (combination != "None") {
    active <- unlist(strsplit(combination, "\\+"))
    out[active] <- 1
  }
  return(out)
}


# Create a data frame suitable for ComplexUpset
upset_data <- do.call(rbind, lapply(perf_df$combination, expand_combination, callers=callers))
upset_data <- cbind(perf_df, upset_data)
upset_data <- na.omit(upset_data)

# Now upset_data has columns for each caller (binary), plus TP, FN, FP, Recall, Precision, F1.
# We can use these as annotations in ComplexUpset.
library(ComplexUpset)

upset_data <- upset_data %>%
  mutate(combination = rownames(upset_data)) %>%  # Add combination names if needed
  arrange(desc(F1))  # Sort by F1 Score in descending order

# Define the sets of interest
sets <- c("TEforest", "PopoolationTE", "RetroSeq", "TEFLoN", "TEMP", "TEMP2")
# Compute intersections and add F1, Recall, and Precision metrics
upset_data_prepared <- upset_data %>%
  mutate(
    combination = apply(select(., all_of(sets)), 1, function(x) paste(names(x)[x == 1], collapse = ", "))
  )

# Aggregate metrics by intersection
intersection_summary <- upset_data_prepared %>%
  group_by(combination) %>%
  summarize(
    Recall = mean(Recall, na.rm = TRUE),
    Precision = mean(Precision, na.rm = TRUE),
    F1 = mean(F1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(F1))

# Generate the ordered list of intersections based on F1
ordered_intersections <- intersection_summary %>%
  pull(combination) %>%
  strsplit(", ")  # Convert back to list of sets for ComplexUpset

# Prepare data and filter zero intersections
#df_for_upset_sets_only_tp <- df_for_upset_sets_only_tp %>%
#  filter(rowSums(select(., starts_with("TE"))) > 0)  # Adjust as needed

unique_tp_upset <- upset(
  df_for_upset_sets_only_tp,
  intersect = c("TEforest", "PopoolationTE", "RetroSeq", "TEFLoN", "TEMP", "TEMP2"),
  set_sizes = upset_set_size() + ggtitle("Total true positives"),
  themes = upset_modify_themes(
    list(
      'intersect' = theme(axis.title.x = element_text(size = 14, hjust = 0.5))
    )
  )
)+
  theme(
    axis.text.x = element_blank(), axis.title.x = element_blank()  # Ensure x-axis text is removed
  )

# Prepare data and filter zero intersections
#df_for_upset_sets_only_fp <- df_for_upset_sets_only_fp %>%
#  filter(rowSums(select(., starts_with("TE"))) > 0)  # Adjust as needed

unique_fp_upset <- upset(
  df_for_upset_sets_only_fp,
  intersect = c("TEforest", "PopoolationTE", "RetroSeq", "TEFLoN", "TEMP", "TEMP2"),
  set_sizes = upset_set_size() + ggtitle("Total false positives"),
  themes = upset_modify_themes(
    list(
      'intersect' = theme(axis.title.x = element_text(size = 14, hjust = 0.5))
    )
  )
) +
  theme(
    axis.text.x = element_blank(), axis.title.x = element_blank()  # Ensure x-axis text is removed
  )

# Create the UpSet plot
iamupset <- upset(
  upset_data,
  intersect = sets,
  intersections = ordered_intersections,
  sort_intersections = FALSE,  # Use the custom order provided
  base_annotations = list(
    # Replace the default intersection size with F1 Score
    'F1 Score' = (
      ggplot(mapping = aes(x = intersection, y = F1)) +
      geom_bar(stat = 'identity', fill = '#F8766D') +  # Set F1 Score color to red
      theme_minimal() +
      labs(y = 'F1 Score')  +
      scale_y_continuous(limits = c(0.5, 1), oob = scales::squish) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())  # Remove x-axis text
    ),
    'Precision' = (
      ggplot(mapping = aes(x = intersection, y = Precision)) +
      geom_bar(stat = 'identity', fill = '#00BA38') +  # Set Precision color to green
      theme_minimal() +
      labs(y = 'Precision') +
      scale_y_continuous(limits = c(0.8, 1), oob = scales::squish) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())  # Remove x-axis text
    ),
    'Recall' = (
      ggplot(mapping = aes(x = intersection, y = Recall)) +
      geom_bar(stat = 'identity', fill = '#619CFF') +  # Set Recall color to blue
      theme_minimal() +
      labs(y = 'Recall') +
      scale_y_continuous(limits = c(0.4, 1), oob = scales::squish) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())  # Remove x-axis text
    )
  ),
  set_sizes = FALSE  # Remove the set size bars
) +
  theme(
    axis.text.x = element_blank(), axis.title.x = element_blank()  # Ensure x-axis text is removed
  )


length_upset_fig <- ggarrange(TElength_comparision, iamupset, 
                    labels = c("A", "B"),
                    font.label = list(size = 22, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    heights = c(0.333, 0.666))

length_upset_fig <- ggarrange(unique_tp_upset, unique_fp_upset, 
                    labels = c("A", "B"),
                    font.label = list(size = 22, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2)


ggsave("/nas/longleaf/home/adaigle/TEforest/plots/upset.svg", plot = length_upset_fig, width = 8.5, height = 10, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/upset.jpg", plot = length_upset_fig, width = 8.5, height = 10, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/upset.tiff", plot = length_upset_fig, width = 8.5, height = 10, dpi = 300)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/TEratio.svg", plot = TEratio_comparision, width = 8.5, height = 5, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/TEratio.jpg", plot = TEratio_comparision, width = 8.5, height = 5, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/TEratio.tiff", plot = TEratio_comparision, width = 8.5, height = 5, dpi = 300)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/upset_tpfn.svg", plot = length_upset_fig, width = 8.5, height = 10, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/upset_tpfn.jpg", plot = length_upset_fig, width = 8.5, height = 10, dpi = 300)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/upset_tpfn.tiff", plot = length_upset_fig, width = 8.5, height = 10, dpi = 300)


te_data_long_mygenomes


#compare benchmarking data
load("/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/2L_2R_plots/JUT-008_MUN-009.RData")
