library(tidyverse)
library(GenomicRanges)
library(viridis)

colors <- c(
  "truth" = "#000000",      # Black for 'truth'
  "TEforest" = "#F8766D", 
  "temp2" = "#CD9600", 
  "retroseq" = "#7CAE00", 
  "temp" = "#00BE67", 
  "popoolationte" = "#00BFC4", 
  "teflon" = "#00A9FF", 
  "popoolationte2" = "#C77CFF", 
  "tepid" = "#FF61CC"
)



ISO1 <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/ISO-1_Ref_Coord.bed") #%>%

dir_path <- "/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates"

# List all files in the directory
files <- list.files(path = dir_path, pattern = "_Ref_Coord.bed$", full.names = TRUE)

dspr <- c(paste0("A", 1:7), paste0("B", 1:4), "B6", "B8")
#dspr <- "RAL"

files <- unlist(lapply(files, function(file) grep(paste(dspr, collapse = "|"), file, value = TRUE)))

# Function to read a file and return a list with file name and data frame
read_bed_file <- function(file_path) {
  # Extract the base name without the suffix
  base_name <- sub("_Ref_Coord.bed$", "", basename(file_path))
  # Read the BED file into a data frame
  bed_df <- read.table(file_path, header = FALSE)
  # Return a list with the name and the data frame
  list(name = base_name, file_path = file_path, df = bed_df)
}

read_mcclintock_het_format <- function(path) {
  #this function reads in the alternative format of mcclintock files. 
  #it needs to be a separate funciton so it outputs a df to a tibble (readLines is tricky)
  #has been modified in this script to only keep calls on chrs. 2L and 2R
  lines <- readLines(path)
  data <- strsplit(lines, "\t|\\|")
  df <- data.frame(do.call(rbind, data))
  df <- df[-1,]
  #df <- df %>% filter(df$X1%in%c("2L","2R"))
  #df <- df %>% filter(df$X1%in%c("3L","3R", "X"))
  return(df)
}

extend <- function(x, upstream=0, downstream=0) {
    #if (any(strand(x) == "*"))
    #    warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

fold_vector <- function(vec) {
  n <- length(vec)
  half_n <- ceiling(n / 2)
  folded_vec <- vec[1:half_n] + rev(vec)[1:half_n]
  return(folded_vec)
}

# Apply the function to all files and combine the results into a data frame
dfs_list <- lapply(files, read_bed_file)

# Create a data frame of data frames
dfs_df <- tibble(
  name = sapply(dfs_list, function(x) x$name),
  file_path = sapply(dfs_list, function(x) x$file_path),
  df = I(sapply(dfs_list, function(x) x$df %>% filter(V7=="roo"), simplify = FALSE)),
)

ISO1_ids <- ISO1$V4

dfs_df <- dfs_df %>%
  filter(name!="ISO-1") %>%
  mutate(IDS = map(df, ~ pluck(.x, "V4") %>% as.character),
    reference_IDS = map(IDS, ~ .x[.x %in% ISO1_ids]),
    nonreference_IDS = map(IDS, ~ .x[!(.x %in% ISO1_ids)])
  )

reference_ids_vector <- unlist(dfs_df$reference_IDS)
# Count the frequency of each ID
ref_id_counts <- tibble(reference_ids = reference_ids_vector) %>%
  count(reference_ids, name = "frequency")

ggplot(ref_id_counts, aes(x = frequency)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(title = "Histogram of Reference ID Overlaps",
       x = "Number of Overlaps",
       y = "Frequency")

nonreference_ids_vector <- unlist(dfs_df$nonreference_IDS)
# Count the frequency of each ID
nonref_id_counts <- tibble(nonreference_IDS = nonreference_ids_vector) %>%
  count(nonreference_IDS, name = "frequency")

ggplot(nonref_id_counts, aes(x = frequency)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(title = "Histogram of Noneference ID Overlaps",
       x = "Number of Overlaps",
       y = "Frequency")

#normalized version
ggplot(nonref_id_counts, aes(x = frequency)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, color = "black") +
  labs(title = "Histogram of Nonreference ID Overlaps",
       x = "Number of Overlaps",
       y = "Density")

nonref_sfs <- as.numeric(table(nonref_id_counts$frequency))
nonref_sfs_normalized <- nonref_sfs / sum(nonref_sfs)
nonref_sfs_normalized_df <- as.data.frame(nonref_sfs_normalized)
nonref_sfs_normalized_df$Row <- as.numeric(rownames(nonref_sfs_normalized_df))

df_long <- pivot_longer(nonref_sfs_normalized_df, 
  cols = -Row, 
  names_to = "Column", 
  values_to = "Value")
df_long$Column <- factor(df_long$Column)


ggplot(df_long, aes(x = Row, y = Value, color = Column, fill = Column, group = Column)) +
  #geom_bar(stat = "identity", colour = "black", position = "dodge") +
  geom_line() +
  scale_y_log10() +
  scale_fill_viridis_d(option = "viridis") +
  scale_color_manual(values = rep("black", length(unique(df_long$Column)))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", title = "log scale SFS") +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25), 
        strip.text = element_text(size = 15), 
        plot.title = element_text(size = 25), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15))

nonref_sfs_normalized_df_folded <- fold_vector(nonref_sfs_normalized)
nonref_sfs_normalized_df_folded <- as.data.frame(nonref_sfs_normalized_df_folded)
nonref_sfs_normalized_df_folded$Row <- as.numeric(rownames(nonref_sfs_normalized_df_folded))

df_long <- pivot_longer(nonref_sfs_normalized_df_folded, 
  cols = -Row, 
  names_to = "Column", 
  values_to = "Value")
df_long$Column <- factor(df_long$Column)


#ISO1_gr <- GRanges(
#    seqnames = ISO1$V1,
#    ranges = IRanges(start = ISO1$V2, end = ISO1$V3),
#    mcols=ISO1$V7
#    )
#
## Find overlaps
#overlaps <- findOverlaps(ISO1_gr, ISO1_gr)
#
## Count overlaps for each range
#overlap_counts <- countOverlaps(ISO1_gr, ISO1_gr)
#
## Create a data frame of counts
#df <- data.frame(counts = overlap_counts)
#
## Summarize the counts
#sfs <- df %>%
#  group_by(counts) %>%
#  summarize(frequency = n()) %>%
#  arrange(counts)
#
#
##appyling to a list of genomes rather than one
#dfs_df <- dfs_df %>%
#  filter(name!="ISO-1") %>%
#  mutate(IDS = map(df, ~ pluck(.x, "V4") %>% as.character),
#    reference_IDS = map(IDS, ~ .x[.x %in% ISO1_ids]),
#    nonreference_IDS = map(IDS, ~ .x[!(.x %in% ISO1_ids)])
#  )
#
#dfs_df <- dfs_df %>% mutate(
#    results = map(df, # reduce to get rid of nested TEs of same type
#        ~ GRanges(
#            seqnames = .x$V1,
#            ranges = IRanges(start = .x$V2, end = .x$V3),
#            mcols = .x$V7
#)),
#    reference = map(results,
#        ~ as.data.frame(subsetByOverlaps(.x, ISO1_gr, type=c("equal"), ignore.strand=TRUE))),
#    nonreference = map(results, #extended bc a tp is within 500 bp of a tp
#        ~ as.data.frame(subsetByOverlaps(.x, ISO1_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE)))
#)
#
#
#split_dfs <- dfs_df %>%
#  unnest(cols = c(nonreference)) %>%
#  group_by(mcols) %>%
#  group_split()
#
## Convert each data frame in split_dfs to GRanges and combine
#gr_list <- map(split_dfs, ~ GRanges(
#  seqnames = .x$seqnames,
#  ranges = IRanges(start = .x$start, end = .x$end),
#  mcols = .x$mcols
#))
#
## Combine all GRanges into one
#combined_gr <- do.call(c, gr_list)
## Find overlaps
#overlaps <- findOverlaps(combined_gr, combined_gr)
#
## Count overlaps for each range
#overlap_counts <- countOverlaps(combined_gr, combined_gr)
#
## Create a data frame of counts
#df <- data.frame(counts = overlap_counts)
#
## Summarize the counts to create the site frequency spectrum
#sfs <- df %>%
#  group_by(counts) %>%
#  summarize(frequency = n()) %>%
#  arrange(counts)
#
## Print the site frequency spectrum
#print(sfs)
#


#genome1 <- "ISO-1"
#genome2 <- "JUT-008"
euchromatin_coordinates_path <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/euchromatin.txt"
#caller_name <- "TEforest_regressor"
#caller_name2 <- "TEforest_classifier"
#plt_dir <- "/nas/longleaf/home/adaigle/work/test_TEforest/Iso1_reference_training/2L_2R_plots"
#basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/Iso1_reference_training"
euchromatin_coordinates <- makeGRangesFromDataFrame(read.table(euchromatin_coordinates_path), seqnames.field="V1", start.field="V2", end.field="V3")

#dspr <- c("A1", "A2", "A3")
# Process each genome and combine results
mcclintock_results <- map_dfr(dspr, function(genome) {
  result_path <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/from_the_ashes/", genome, "_1/results/")
  mcclintock_results <- tibble(
    caller = list.files(result_path)[!grepl("summary", list.files(result_path))],
    result_file = unlist(map(paste0(result_path, caller), ~list.files(.x, pattern = "_nonredundant.bed"))),
    fullpath = paste0(result_path, caller, "/", result_file),
    data = lapply(fullpath, read_mcclintock_het_format),
    df_list = map(data, ~split(.x, .x[["X4"]])), # split dataframes based on individual TEs
    genome = str_extract(fullpath, "(?<=/)[A-Z0-9-]+(?=_1/)")
  ) %>%
  unnest(cols = df_list) %>%
  mutate(
    TE = unlist(map(df_list, ~.x$X4[1])),
    length = unlist(map(df_list, ~length(.x$X1)))
  )
  
  return(mcclintock_results)
})

mcclintock_results <- subset(mcclintock_results, length !=0) %>% # might bias results upwards
  mutate(
  nonreference = map(df_list, ~ filter(.x, X5 == "non-reference")), 
  reference = map(df_list, ~ filter(.x, X5 == "reference")),
  #0 to 1 based
  granges = map(df_list, ~ GRanges(
    seqnames = .x$X1,
    ranges = IRanges(start = as.numeric(.x$X2)+1, end = as.numeric(.x$X3)))),
  nonref_gr = map(nonreference, ~ GRanges(
    seqnames = .x$X1,
    ranges = IRanges(start = as.numeric(.x$X2)+1, end = as.numeric(.x$X3)),
    TE = .x$X4, heterozygosity= .x$X6)))

mcclintock_results <- mcclintock_results %>% mutate(
  nonref_gr_filter = map(nonref_gr, ~ subsetByOverlaps(GenomicRanges::reduce(extend(.x, 500,500)), euchromatin_coordinates, ignore.strand = TRUE))
)

mcclintock_results_notruth <- mcclintock_results %>% filter(caller != "truth") 
truth_calls <- mcclintock_results %>% filter(caller == "truth") %>% 
  select(genome, TE, truth_calls = nonref_gr_filter)

#mcclintock_results_truthcalls <- left_join(mcclintock_results_notruth, truth_calls,by = join_by(genome, TE)) %>%
#  mutate(
#
#  )



nested_results <- mcclintock_results %>% select(caller,TE,nonref_gr_filter) %>% 
  group_by(caller,TE) %>%
  nest() %>%
  mutate(
    combined_gr = map(data, ~ do.call(c, unname(.x$nonref_gr_filter)))
  )

nested_results_sfs <- nested_results %>% mutate(
    combined_gr_reduced = map(combined_gr, GenomicRanges::reduce),
    overlap_counts = map2(combined_gr_reduced, combined_gr, ~ countOverlaps(.x, .y)),
    sfs = map(overlap_counts, ~ {
      df <- data.frame(counts = .x)
      df %>%
        group_by(counts) %>%
        summarize(frequency = n()) %>%
        arrange(counts)
    })
  )

nested_results_sfs_new <- nested_results %>%
  mutate(
    combined_gr_reduced = map(combined_gr, GenomicRanges::reduce),
    overlap_counts = map2(combined_gr_reduced, combined_gr, ~ countOverlaps(.x, .y)),
    combined_gr_reduced = map2(combined_gr_reduced, overlap_counts, ~ {
      mcols(.x)$overlap_counts <- .y
      .x
    }),
    sfs = map(overlap_counts, ~ {
      df <- data.frame(counts = .x)
      
      # Create a data frame with all counts from 1 to 13
      all_counts <- data.frame(counts = 1:13)
      
      # Group by counts and summarize frequencies
      freq_df <- df %>%
        group_by(counts) %>%
        summarize(frequency = n()) %>%
        arrange(counts)
      
      # Join with all_counts to ensure all counts from 1 to 13 are present
      complete_freq_df <- all_counts %>%
        left_join(freq_df, by = "counts") %>%
        replace_na(list(frequency = 0))
      
      complete_freq_df
    })
  )

nested_results_sfs_new_notruth <- nested_results_sfs_new %>% filter(caller != "truth") 

truth_calls_counts <- nested_results_sfs_new %>% ungroup() %>% filter(caller == "truth") %>% 
  select(TE, truth_calls = combined_gr_reduced) %>% group_by(TE)

nested_results_truthjoin <- left_join(nested_results_sfs_new_notruth, truth_calls_counts, join_by(TE))




# Function to compare two GRanges objects
compare_granges <- function(grange1, grange2) {
  results <- data.frame(
    grange1_idx = integer(),
    grange2_idx = integer(),
    overlap_difference = numeric()
  )
  
  # Iterate over each range in grange1
  for (i in seq_along(grange1)) {
    overlaps <- findOverlaps(grange1[i], grange2)
    
    # If there are overlapping ranges
    if (length(overlaps) > 0) {
      for (j in subjectHits(overlaps)) {
        diff <- mcols(grange1)$overlap_counts[i] - mcols(grange2)$overlap_counts[j]
        results <- rbind(results, data.frame(
          grange1_idx = i,
          grange2_idx = j,
          overlap_difference = diff
        ))
      }
    }
  }
  
  return(results$overlap_difference)
}

compare_granges_filtered <- function(grange1, grange2, allowed_intervals) {
  results <- data.frame(
    grange1_idx = integer(),
    grange2_idx = integer(),
    overlap_difference = numeric()
  )
  
  # Helper function to check if a value falls into any of the allowed intervals
  in_allowed_intervals <- function(value, intervals) {
    any(sapply(intervals, function(interval) {
      value >= interval[1] && value <= interval[2]
    }))
  }
  
  # Iterate over each range in grange1
  for (i in seq_along(grange1)) {
    overlaps <- findOverlaps(grange1[i], grange2)
    
    # If there are overlapping ranges
    if (length(overlaps) > 0) {
      for (j in subjectHits(overlaps)) {
        gr2_count <- mcols(grange2)$overlap_counts[j]
        
        # Check if gr2_count is in allowed intervals
        if (in_allowed_intervals(gr2_count, allowed_intervals)) {
          diff <- mcols(grange1)$overlap_counts[i] - gr2_count
          results <- rbind(results, data.frame(
            grange1_idx = i,
            grange2_idx = j,
            overlap_difference = diff
          ))
        }
      }
    }
  }
  
  return(results$overlap_difference)
}

compare_granges_filtered <- function(grange1, grange2, allowed_intervals = list(c(1,3), c(5,12))) {
  results <- data.frame(
    grange1_idx = integer(),
    grange2_idx = integer(),
    overlap_difference = numeric()
  )
  
  # Helper function to check if a value falls into any of the allowed intervals
  in_allowed_intervals <- function(value, intervals) {
    any(sapply(intervals, function(interval) {
      value >= interval[1] && value <= interval[2]
    }))
  }
  
  # Iterate over each range in grange1
  for (i in seq_along(grange1)) {
    overlaps <- findOverlaps(grange1[i], grange2)
    
    # If there are overlapping ranges
    if (length(overlaps) > 0) {
      has_valid_overlap <- FALSE
      for (j in subjectHits(overlaps)) {
        gr2_count <- mcols(grange2)$overlap_counts[j]
        
        # Check if gr2_count is in allowed intervals
        if (in_allowed_intervals(gr2_count, allowed_intervals)) {
          diff <- mcols(grange1)$overlap_counts[i] - gr2_count
          results <- rbind(results, data.frame(
            grange1_idx = i,
            grange2_idx = j,
            overlap_difference = diff
          ))
          has_valid_overlap <- TRUE
        }
      }
    
    } else {
      # No overlaps found: record the grange1 overlap_counts value
      results <- rbind(results, data.frame(
        grange1_idx = i,
        grange2_idx = NA,
        overlap_difference = mcols(grange1)$overlap_counts[i]
      ))
    }
  }
  
  return(results$overlap_difference)
}

test1 <- nested_results_truthjoin$combined_gr_reduced[[951]]
test2 <- nested_results_truthjoin$truth_calls[[951]]
comparison_results <- compare_granges(test1, test2)



nested_results_truthjoin <- nested_results_truthjoin %>% 
  filter(!caller %in% c("allseeingeye", "allseeingeye2", "allseeingeye3", "allseeingeye4", "allseeingeye5", "2L_test", "2L_test_bps", "tepid", "te-locate", "ngs_te_mapper2")) %>%
  filter(!map_lgl(truth_calls, is.null))

nested_results_truthjoin <- nested_results_truthjoin %>%
  mutate(
    # predicted freq -true freq
    overlap_differences = map2(combined_gr_reduced, truth_calls, compare_granges)
  )

nested_results_truthjoin_quicktest <- nested_results_truthjoin %>% filter(TE=="roo")

nested_results_truthjoin_common <- nested_results_truthjoin_quicktest %>%
  mutate(
    overlap_differences = map2(
      combined_gr_reduced,
      truth_calls,
      compare_granges_filtered,
      allowed_intervals = list(c(2,20))
    )
  )

nested_results_truthjoin_rare <- nested_results_truthjoin %>%
  mutate(
    overlap_differences = map2(
      combined_gr_reduced,
      truth_calls,
      compare_granges_filtered,
      allowed_intervals = list(c(0,1))
    )
  )
saveRDS(nested_results_truthjoin, "/nas/longleaf/home/adaigle/TEforest/plots/nested_results_truthjoin_signed.rds")

nested_results_truthjoin <- readRDS("/nas/longleaf/home/adaigle/TEforest/plots/nested_results_truthjoin_signed.rds")
# Summarize the overlap differences for each caller
freq_comparision_summary <- nested_results_truthjoin %>%
  #filter(TE == "roo") %>%
  unnest(overlap_differences) %>%
  group_by(caller) %>%
  summarize(overlap_differences_combined = list(overlap_differences))

freq_comparision_summary_wdif <- freq_comparision_summary %>% mutate(
  avg_difference = unlist(map(overlap_differences_combined, ~ mean(.x))),
  sd_difference = unlist(map(overlap_differences_combined, ~ sd(.x) / sqrt(length(.x))))
) %>%
  mutate(caller = factor(caller, levels = c("temp", "TEforest", "temp2", "retroseq", "teflon", "popoolationte", "popoolationte2", "tepid"))) %>%
  group_by(caller)

freqdif <- ggplot(freq_comparision_summary_wdif, aes(x=caller, y=avg_difference/13, fill=caller))+ 
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin = avg_difference/13 - sd_difference/13, ymax = avg_difference/13 + sd_difference/13), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        legend.position = "right",
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size = 16)) +
  labs(x = "Caller",
       y = "Predicted Frequency - True frequency")

freq_comparision_summary_unnested <- nested_results_truthjoin %>%
  unnest(overlap_differences) %>%
  # Make sure caller is a factor with the desired order
  mutate(caller = factor(caller, 
                         levels = c("temp", "TEforest", "temp2", 
                                    "retroseq", "teflon", 
                                    "popoolationte", "popoolationte2", 
                                    "tepid")))

# Plot the distribution with a violin (or box) plot
ggplot(freq_comparision_summary_unnested, 
       aes(x = caller, y = overlap_differences / 13, fill = caller)) +
  geom_boxplot() +
  # geom_boxplot(width = 0.1, alpha=0.5)  # add or replace violin with boxplot
  scale_fill_manual(values = colors) + # use the same colors as in your original code
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "right",
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14), 
    strip.text = element_text(size = 16)
  ) +
  labs(
    x = "Caller",
    y = "Predicted Frequency - True frequency"
  ) + scale_y_log10() 

ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/plots/freqdif.png"), freqdif, width=8, height = 5)
ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/plots/freqdif.svg"), freqdif, width=8, height = 5)
ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/plots/freqdif.tiff"), freqdif, width=8, height = 5)

nested_results_truthjoin <- readRDS("/nas/longleaf/home/adaigle/TEforest/plots/nested_results_truthjoin_signed.rds")
# Summarize the overlap differences for each caller
freq_comparision_summary_common <- nested_results_truthjoin_common %>%
  #filter(TE == "roo") %>%
  unnest(overlap_differences) %>%
  group_by(caller) %>%
  summarize(overlap_differences_combined = list(overlap_differences))

freq_comparision_summary_wdif_common <- freq_comparision_summary_common %>% mutate(
  avg_difference = unlist(map(overlap_differences_combined, ~ mean(.x))),
  sd_difference = unlist(map(overlap_differences_combined, ~ sd(.x) / sqrt(length(.x))))
) %>%
  mutate(caller = factor(caller, levels = c("temp", "TEforest", "temp2", "retroseq", "teflon", "popoolationte", "popoolationte2", "tepid"))) %>%
  group_by(caller)

freqdif_common <- ggplot(freq_comparision_summary_wdif_common, aes(x=caller, y=avg_difference/13, fill=caller))+ 
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin = avg_difference/13 - sd_difference/13, ymax = avg_difference/13 + sd_difference/13), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        legend.position = "right",
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size = 16)) +
  labs(x = "Caller",
       y = "True frequency - Predicted Frequency")

# Summarize the overlap differences for each caller
freq_comparision_summary_rare <- nested_results_truthjoin_rare %>%
  #filter(TE == "roo") %>%
  unnest(overlap_differences) %>%
  group_by(caller) %>%
  summarize(overlap_differences_combined = list(overlap_differences))

freq_comparision_summary_wdif_rare <- freq_comparision_summary_rare %>% mutate(
  avg_difference = unlist(map(overlap_differences_combined, ~ mean(.x))),
  sd_difference = unlist(map(overlap_differences_combined, ~ sd(.x) / sqrt(length(.x))))
) %>%
  mutate(caller = factor(caller, levels = c("temp", "TEforest", "temp2", "retroseq", "teflon", "popoolationte", "popoolationte2", "tepid"))) %>%
  group_by(caller)

freqdif_rare <- ggplot(freq_comparision_summary_wdif_rare, aes(x=caller, y=avg_difference/13, fill=caller))+ 
  geom_bar(stat = "identity", position="dodge") +
  geom_errorbar(aes(ymin = avg_difference/13 - sd_difference/13, ymax = avg_difference/13 + sd_difference/13), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        legend.position = "right",
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size = 16)) +
  labs(x = "Caller",
       y = "True frequency - Predicted Frequency")




###code for SFS plotting and misc
# Unnest the overlap_differences_combined column
unnested_results <- freq_comparision_summary %>%
  unnest(overlap_differences_combined)
# Create the violin plot
ggplot(unnested_results, aes(x = caller, y = overlap_differences_combined, color = caller)) +
  #geom_violin() +
  geom_boxplot() +
  #geom_histogram() +
  #geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "Violin Plot of Overlap Differences by Caller",
       x = "Caller",
       y = "Overlap Differences") +
  theme_minimal() 
#roo <- nested_results_sfs %>% filter(TE=='roo')

nested_results_sfs_alltes <-  nested_results_sfs_new %>%
  group_by(caller) %>%
  summarize(
    sfs = list(purrr::reduce(sfs, ~ full_join(.x, .y, by = "counts") %>%
                        mutate(frequency = rowSums(across(starts_with("frequency"), ~ replace_na(.x, 0)))) %>%
                        select(counts, frequency))),
    .groups = 'drop'
  )

filter_and_transform_bycaller <- function(df, caller_choice) {
  # Filter by the chosen caller
  filtered_df <- df %>% filter(caller == caller_choice)
  
  # Extract the sfs data and bind it with TE names
  sfs_data <- filtered_df %>%
    select(TE, sfs) %>%
    unnest(sfs) %>%
    group_by(TE)
  
  return(sfs_data)
}

temp2 <- filter_and_transform_bycaller(nested_results_sfs, "temp2" )

normalize_sfs <- function(df) {
  df$frequency <- df$frequency / sum(df$frequency)
  return(df)
}

nested_results_sfs_normalized <- nested_results_sfs %>% mutate(
  sfs = map(sfs, normalize_sfs)
)
nested_results_sfs_alltes_normalized <- nested_results_sfs_alltes %>% mutate(
  sfs = map(sfs, normalize_sfs)
)

temp2_normalized <- filter_and_transform_bycaller(nested_results_sfs_normalized, "temp2" )
truth_normalized <- filter_and_transform_bycaller(nested_results_sfs_normalized, "truth" )

ggplot(truth_normalized, aes(x = factor(counts), y = frequency, fill = TE)) +
  geom_bar(stat = "identity", colour = "black", position = "dodge") +
  scale_fill_viridis_d(option="viridis") +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", title = paste0("folded SFS")) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.position = "none")


filter_and_transform_byTE <- function(df, TE_choice) {
  # Filter by the chosen caller
  filtered_df <- df %>% filter(TE == TE_choice)
  
  # Extract the sfs data and bind it with TE names
  sfs_data <- filtered_df %>%
    select(caller, sfs) %>%
    unnest(sfs) %>%
    group_by(caller)
  
  return(sfs_data)
}

#roo <- filter_and_transform_byTE(nested_results_sfs, "roo" )
plot_te_allcallers <- function(TE) {
  roo_normalized <- filter_and_transform_byTE(nested_results_sfs_normalized, TE )
  roo_normalized <- roo_normalized %>% filter(!caller %in% c("allseeingeye", "allseeingeye2", "allseeingeye3", "allseeingeye4", "allseeingeye5", "2L_test", "2L_test_bps", "tepid", "te-locate", "ngs_te_mapper2"))
  # Define the color assignments, including 'truth'
  colors <- c(
    "truth" = "#000000",      # Black for 'truth'
    "TEforest" = "#F8766D", 
    "temp2" = "#CD9600", 
    "retroseq" = "#7CAE00", 
    "temp" = "#00BE67", 
    "popoolationte" = "#00BFC4", 
    "teflon" = "#00A9FF", 
    "popoolationte2" = "#C77CFF", 
    "tepid" = "#FF61CC"
  )

  # Reorder the levels of the caller factor to match the order in colors, including 'truth'
  roo_normalized$caller <- factor(roo_normalized$caller, levels = names(colors))

  final_plot <- ggplot(roo_normalized, aes(x = factor(counts), y = frequency, fill = caller)) +
    geom_bar(stat = "identity", colour = "black", position = "dodge") +
    scale_fill_manual(values = colors) +  # Use scale_fill_manual to apply the custom colors
    labs(x = "Derived allele count", y = "Proportion of polymorphisms", title = "") +
    theme(
      axis.text.x = element_text(size = 15), 
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25), 
      strip.text = element_text(size = 15), 
      plot.title = element_blank(), 
      legend.position = "bottom"
    ) + theme_minimal()
  return(final_plot)
}

copia <- plot_te_allcallers("Copia")
f <- plot_te_allcallers("F_element")
roo <- plot_te_allcallers("roo")
jockey <- plot_te_allcallers("jockey")

ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/roo_SFS.pdf"), roo, width=14)
ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/f_SFS.pdf"), f, width=14)
ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/copia_SFS.pdf"), copia, width=14)

# Arrange the plots with one on top and two on the bottom
combined_plot <- ggarrange(
  freqdif,                          # Top plot
  ggarrange(roo, jockey,            # Two bottom plots
            labels = c("B", "C"),   # Labels for the bottom plots
            ncol = 2, nrow = 1,     # Two plots in one row
            common.legend = TRUE,   # Share a common legend
            legend = "bottom"),
  labels = c("A", ""),              # Label for the top plot
  ncol = 1, nrow = 2,               # One plot per row
  common.legend = TRUE,             # Share a common legend
  legend = "bottom",                # Place legend at the bottom
  font.label = list(size = 22, color = "black", face = "bold")
)

ggsave(paste0("/nas/longleaf/home/adaigle/TEforest/plots/sfs_crap.jpg"), combined_plot, width=10.5, height=10, dpi=300)


sfs_data <- nested_results_sfs_alltes_normalized %>%
  select(caller, sfs) %>%
  unnest(sfs) %>%
  group_by(caller) %>%
  filter(!caller %in% c("allseeingeye", "allseeingeye2", "allseeingeye3", "allseeingeye4", "allseeingeye5", "2L_test", "2L_test_bps"))

ggplot(sfs_data, aes(x = factor(counts), y = frequency, fill = caller)) +
  geom_bar(stat = "identity", colour = "black", position = "dodge") +
  scale_fill_viridis_d(option="viridis") +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", title = paste0("unfolded SFS")) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.position = "bottom")


save.image(file = "/nas/longleaf/home/adaigle/TEforest/workflow/scripts/sfs_analysis_slow_done.RData")
load("/nas/longleaf/home/adaigle/TEforest/workflow/scripts/sfs_analysis_slow_done.RData")

#code for moving teforest and truth files to mcclintock output dirs. 
# Base directory paths
basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/inference_test_SFS_allgenomes"
destination_base_folder <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/from_the_ashes"

# Loop through each genome
for (genome in dspr) {
  # Construct the source file path
  source_file <- paste0(basedir_outputs_path, "/output/", genome, "_TEforest_bps_nonredundant.bed")
  # Specify the destination directory
  destination_folder <- paste0(destination_base_folder, "/", genome, "_1/results/", "TEforest")
  # Create the destination directory if it does not exist
  if (!dir.exists(destination_folder)) {
    dir.create(destination_folder)
  }
  # Construct the destination file path
  destination_file <- paste0(destination_folder, "/", genome, "_1_", "TEforest", "_nonredundant.bed")
  # Use file.copy() to copy the file
  file.copy(from = source_file, to = destination_file, overwrite = TRUE)
  # Print status message
  cat("Copied file for genome:", genome, "\n")
}

dfs_list <- lapply(files, read_bed_file)

# Create a data frame of data frames
dfs_df <- tibble(
  name = sapply(dfs_list, function(x) x$name),
  file_path = sapply(dfs_list, function(x) x$file_path),
  df = I(sapply(dfs_list, function(x) x$df, simplify = FALSE)),
)

dfs_df <- dfs_df %>%
  filter(name!="ISO-1") %>%
  mutate(IDS = map(df, ~ pluck(.x, "V4") %>% as.character),
    reference_IDS = map(IDS, ~ .x[.x %in% ISO1_ids]),
    nonreference_IDS = map(IDS, ~ .x[!(.x %in% ISO1_ids)])
  )

add_reference_column <- function(df, reference_ids) {
  df %>%
    mutate(V8 = if_else(V4 %in% reference_ids, "reference", "non-reference"),
           V7 = str_replace_all(V7, "-", "_"))
}

dfs_df <- dfs_df %>%
  mutate(df = map2(df, reference_IDS, ~ add_reference_column(.x,.y)))

dfs_df<- dfs_df %>% mutate(
  mcclintock_format = map2(df, name, ~ data.frame(
    seqnames = .x$V1,
    start = .x$V2 - 1,
    end = .x$V3, 
    TE_string = paste(.x$V7, .x$V8, "1", .y, "truth", "rp", "1", sep="|"), 
    score = 0, 
    strand = "."
  ))
)

for (genome in dspr) {
  mcclintock_format_df <- dfs_df %>% filter(name==genome) %>% pull(mcclintock_format)
  destination_folder <- paste0(destination_base_folder, "/", genome, "_1/results/", "truth")
  if (!dir.exists(destination_folder)) {
    dir.create(destination_folder)
  }
  
  write.table(mcclintock_format_df, 
      file = paste0(destination_folder, "/", genome, "_1_", "truth", "_nonredundant.bed"),
      quote = F, sep = "\t", row.names = F
  )
}