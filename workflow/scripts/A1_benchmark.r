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

#dspr <- c(paste0("A", 1:7), paste0("B", 1:4), "B6", "B8")
dspr <- "A1"

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


#genome1 <- "ISO-1"
#genome2 <- "JUT-008"
euchromatin_coordinates_path <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/euchromatin.txt"
#caller_name <- "TEforest_regressor"
#caller_name2 <- "TEforest_classifier"
#plt_dir <- "/nas/longleaf/home/adaigle/work/test_TEforest/Iso1_reference_training/2L_2R_plots"
#basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/Iso1_reference_training"
euchromatin_coordinates <- makeGRangesFromDataFrame(read.table(euchromatin_coordinates_path), seqnames.field="V1", start.field="V2", end.field="V3")

dspr <- c("A1")
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

A1_truth <- read.table(paste0("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/", "A1", "_Ref_Coord.bed"))
A1_truth$heterozygosity <- "1"
colnames(A1_truth) <- c("seqnames", "start", "end", "ID", "length", "strand", "TE", "heterozygosity")
A1_truth$TE <- gsub("-", "_", A1_truth$TE) #critical
TE_list <- unique(sort(A1_truth$TE))
A1_truth$start <- A1_truth$start #do not need to make 1 based (they already are!)
A1_truth <- A1_truth 
truth_gr <- GRanges(
    seqnames = A1_truth$seqnames,
    ranges = IRanges(start = A1_truth$start, end = A1_truth$end)
    )

benchmark_mapping_results <- mcclintock_results %>% mutate(
    A1_truth_forTE = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x), keep.extra.columns = T )),
    A1_truth_ref = map(A1_truth_forTE,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref = map(A1_truth_forTE, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    A1_truth_nonref_noextend = map(A1_truth_forTE, #extended bc a tp is within 500 bp of a tp
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE)),
    #nonref_gr_filter = map2(nonref_gr_filter, TE, ~ GRanges(seqnames = seqnames(.x), 
    #    ranges = ranges(.x),
    #    TE = .y, heterozygosity=nonref_gr_filter$heterozygosity)),
    nonref_gr_filter = map(nonref_gr, ~ subsetByOverlaps(.x, euchromatin_coordinates, ignore.strand = TRUE)),
    nonreference_freq_filter = map(nonref_gr_filter, ~ .x %>% as.data.frame() %>% filter(heterozygosity >= 0.10) %>% makeGRangesFromDataFrame(keep.extra.columns=T)),
    f1_score = unlist(map2(nonreference_freq_filter, A1_truth_nonref, ~ f1_score(.x, .y))),
    nonref_false_negatives = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE, invert=TRUE)),
    nonref_true_positives = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)),
    nonref_true_positives_truthcoords = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE)),
    nonref_false_positives = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, invert = TRUE)),
    nonref_true_positives_length = unlist(map(A1_truth_nonref,~ length(.x))),
    nonref_calls_true_positives_length = unlist(map(nonref_true_positives,~ length(.x))),
    nonref_calls_false_positives_length = unlist(map(nonref_false_positives,~ length(.x))),
    nonref_calls_false_negatives_length = unlist(map(nonref_false_negatives,~ length(.x))),
    nonref_calls_length= unlist(map(nonreference_freq_filter,~ length(.x))),
    #stats = map2(nonref_gr_filter, A1_truth_nonref, ~ tp_fp_fn(.x, .y))
    precision = nonref_calls_true_positives_length / (nonref_calls_true_positives_length + nonref_calls_false_positives_length),
    recall = nonref_calls_true_positives_length / (nonref_calls_true_positives_length + nonref_calls_false_negatives_length),
    f1_score = 2 * (precision * recall) / (precision + recall)
)
benchmark_mapping_results_nodata <- benchmark_mapping_results_het[c(1,6,7,20,24:30,33, 45:48, 52:53, 55, 58:59, 61:62, 65)]

ISO1_og <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/ISO-1_Ref_Coord.bed") #%>%
ISO1_og_gr <- GRanges(
    seqnames = ISO1_og$V1,
    ranges = IRanges(start = ISO1_og$V2, end = ISO1_og$V3),
    mcols=ISO1_og$V7
    )
truth_nonreference_extend_gr <- extend(subsetByOverlaps(truth_gr, ISO1_og_gr, type=c("equal"), invert=TRUE, ignore.strand=TRUE), 500, 500)
nonref_te_number <- length(unique(truth_nonreference_extend_gr)) #SWAP

benchmark_mapping_results_summary <- benchmark_mapping_results %>%
  group_by(caller) %>%
  summarize(
    #per_fam_mean_f1_score = mean(f1_score, na.rm = TRUE),
    #per_fam_mean_precision = mean(precision, na.rm = TRUE),
    #per_fam_mean_recall = mean(recall, na.rm = TRUE),
    sum_nonref_true_positives_length = sum(nonref_true_positives_length),
    sum_true_positives = sum(nonref_calls_true_positives_length),
    sum_false_positives = sum(nonref_calls_false_positives_length),
    sum_false_negatives = sum(nonref_calls_false_negatives_length),
    recalc_false_negatives = nonref_te_number - sum_true_positives,
    #true_negatives_length = sum(true_negatives_length),
  )

benchmark_mapping_results_summary <- benchmark_mapping_results_summary %>% mutate(
    precision = sum_true_positives / (sum_true_positives + sum_false_positives),
    recall = sum_true_positives / (sum_true_positives + recalc_false_negatives),
    f1_score = 2 * (precision * recall) / (precision + recall)
)

benchmark_mapping_results %>% filter(caller=="TEforest")

tp_df <- benchmark_mapping_results %>%
  filter(caller == "TEforest") %>%
  transmute(TE_row = TE, gr = nonref_true_positives_truthcoords) %>%
  mutate(df = map2(gr, TE_row, ~{
    g <- .x
    if (length(g) == 0) {
      return(tibble(
        seqnames = character(),
        start = integer(),
        stop = integer(),
        strand = character(),
        TE = .y,
        heterozygosity = character()
      ))
    }
    d <- as.data.frame(g)  # has seqnames, start, end, strand + metadata
    tibble(
      seqnames = as.character(d$seqnames),
      start = as.integer(d$start)+500,
      stop  = as.integer(d$end)-500,
      strand = as.character(d$strand),    # <- no Biostrings:: here
      TE = if ("TE" %in% names(d)) as.character(d$TE) else .y,
      heterozygosity = if ("heterozygosity" %in% names(d))
                         as.character(d$heterozygosity) else NA_character_
    )
  })) %>%
  select(df) %>%
  unnest(df)


fn_df <- benchmark_mapping_results %>%
  filter(caller == "TEforest") %>%
  transmute(TE_row = TE, gr = nonref_false_negatives) %>%
  mutate(df = map2(gr, TE_row, ~{
    g <- .x
    if (length(g) == 0) {
      return(tibble(
        seqnames = character(),
        start = integer(),
        stop = integer(),
        strand = character(),
        TE = .y,
        heterozygosity = character()
      ))
    }
    d <- as.data.frame(g)  # has seqnames, start, end, strand + metadata
    tibble(
      seqnames = as.character(d$seqnames),
      start = as.integer(d$start) + 500,
      stop  = as.integer(d$end) - 500,
      strand = as.character(d$strand),
      TE = if ("TE" %in% names(d)) as.character(d$TE) else .y,
      heterozygosity = if ("heterozygosity" %in% names(d))
                         as.character(d$heterozygosity) else NA_character_
    )
  })) %>%
  select(df) %>%
  unnest(df)

library(fs)

# Where the source images live
image_dir <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/testing_truth_dataset/alignment_images"

# Where you want copies to go
tp_out_dir <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/testing_truth_dataset/alignment_images_TP"
fn_out_dir <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/testing_truth_dataset/alignment_images_FN"

# Helper to build filenames and copy
make_image_paths <- function(df, image_dir, sanitize_te = FALSE) {
  df %>%
    mutate(
      start_i = as.integer(start),
      stop_i  = as.integer(stop),
      TE_str  = as.character(TE),
      TE_str  = if (sanitize_te) gsub("[^A-Za-z0-9._-]", "_", TE_str) else TE_str,
      filename = paste0(seqnames, "_", start_i, "_", stop_i, "_", TE_str, ".stacked.png"),
      src = file.path(image_dir, filename)
    )
}

copy_image_set <- function(df, image_dir, out_dir, sanitize_te = FALSE, overwrite = FALSE) {
  dir_create(out_dir)
  df2 <- make_image_paths(df, image_dir, sanitize_te) %>%
    distinct(src, .keep_all = TRUE) %>%        # avoid duplicate copies
    mutate(exists = file_exists(src))

  ok  <- df2 %>% filter(exists)
  miss <- df2 %>% filter(!exists) %>%
    select(seqnames, start, stop, TE, filename, src)

  if (nrow(ok) > 0) {
    # build destination paths explicitly (fast + clear)
    dests <- file.path(out_dir, basename(ok$src))
    file_copy(ok$src, dests, overwrite = overwrite)
  }

  list(
    copied  = ok %>% transmute(src, dest = file.path(out_dir, basename(src))),
    missing = miss
  )
}

# --- Run for your two data frames ---
# assuming you already have tp_df and fn_df from earlier steps
tp_res <- copy_image_set(tp_df, image_dir, tp_out_dir, sanitize_te = FALSE, overwrite = FALSE)
fn_res <- copy_image_set(fn_df, image_dir, fn_out_dir, sanitize_te = FALSE, overwrite = FALSE)

# Quick summaries
message(sprintf("TP: %d copied, %d missing", nrow(tp_res$copied), nrow(tp_res$missing)))
message(sprintf("FN: %d copied, %d missing", nrow(fn_res$copied), nrow(fn_res$missing)))