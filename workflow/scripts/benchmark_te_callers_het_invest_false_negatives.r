#!/usr/bin/env Rscript
# the purpose of this script is to benchmark all TE callers and my caller 
# in an identical fashion. 
# I will benchmark each TE family individually to achive maximum accuracy.
# And calculate the avg distance of the calls from breakpoints. 
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
library(viridis)
library(yardstick)

args <- commandArgs(TRUE)
print(args)
genome1 <- args[1]
genome2 <- args[2]
euchromatin_coordinates_path <- args[3]
caller_name <- args[4]
caller_name2 <- args[5]
plt_dir <- args[6]
coverage <- args[7]
basedir_outputs_path <- getwd()
plt_dir <- paste0(basedir_outputs_path, "/", plt_dir)
if (!dir.exists(plt_dir)) {
  # Create the directory
  dir.create(plt_dir)
}

genome1 <- "JUT-008"
genome2 <- "MUN-009"
euchromatin_coordinates_path <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/euchromatin.txt"
caller_name <- "TEforest_regressor"
caller_name2 <- "TEforest_classifier"
plt_dir <- "/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/2L_2R_plots"
basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X"
coverage <- "50"

training_csv <- read.csv(paste0(basedir_outputs_path, "/3L3RX_classifer/", genome1, "_", genome2, ".csv"))

read_mcclintock_het_format <- function(path) {
  #this function reads in the alternative format of mcclintock files. 
  #it needs to be a separate funciton so it outputs a df to a tibble (readLines is tricky)
  #has been modified in this script to only keep calls on chrs. 2L and 2R
  lines <- readLines(path)
  data <- strsplit(lines, "\t|\\|")
  df <- data.frame(do.call(rbind, data))
  df <- df[-1,]
  df <- df %>% filter(df$X1%in%c("2L","2R"))
  #df <- df %>% filter(df$X1%in%c("3L","3R", "X"))
  return(df)
}

shrink <- function(x, upstream = 0, downstream = 0) {
    if (any(strand(x) == "*")) {
        warning("'*' ranges were treated as '+'")
    }
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) + ifelse(on_plus, upstream, downstream)
    new_end <- end(x) - ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

center <- function(gr) {
  #find center of your grange, return it as a new grange
  center <- mid(gr)
  df <- data.frame(seqnames=as.character(seqnames(gr)),
    ranges=center,
    strand=strand(gr))
  new_gr <- makeGRangesFromDataFrame(df,start.field="ranges", end.field="ranges")
  return(new_gr)
}

confusion_matrix_plus  <- function(calls_Granges_file, truth_Granges_file) {
    # Plots a confusion matrix, creates the matches file for you. 
    # matches_Granges_file <- subsetByOverlaps(calls_Granges_file, truth_Granges_file, ignore.strand = TRUE)
    matches_Granges_file <- subsetByOverlaps(truth_Granges_file, calls_Granges_file, ignore.strand = TRUE)
    true_positives <- length(matches_Granges_file) 
    true_negatives <- length
    false_positives <- length(calls_Granges_file)-length(matches_Granges_file)
    false_negatives <- length(truth_Granges_file)-length(matches_Granges_file)
    target <- c(0,1,0,1)
    prediction <- c(0,0,1,1)
    # n <- c(0,false_positives,false_negatives,true_positives)
    n <- c(0,false_negatives,false_positives,true_positives)
    plot <- plot_confusion_matrix(data.frame(target,prediction, n), 
                        target_col = "target", 
                        prediction_col = "prediction",
                        counts_col = "n")
    return(plot)}

f1_score <- function(calls_Granges_file, truth_Granges_file) {
    # Plots a confusion matrix, creates the matches file for you. 
    # matches_Granges_file <- subsetByOverlaps(calls_Granges_file, truth_Granges_file, ignore.strand = TRUE)
    matches_Granges_file <- subsetByOverlaps(truth_Granges_file, calls_Granges_file, ignore.strand = TRUE)
    true_positives <- length(matches_Granges_file) 
    true_negatives <- length
    false_positives <- length(calls_Granges_file)-length(matches_Granges_file)
    false_negatives <- length(truth_Granges_file)-length(matches_Granges_file)
    precision <- true_positives / (true_positives + false_positives)
    recall <- true_positives / (true_positives + false_negatives)
    F1 <- 2 * (precision * recall) / (precision + recall)
    return(F1)
}

tp_fp_fn <- function(calls_Granges_file, truth_Granges_file) {
    #list with num tp, fp, fn, precision, recall
    matches_Granges_file <- subsetByOverlaps(truth_Granges_file, calls_Granges_file, ignore.strand = TRUE)
    true_positives <- length(matches_Granges_file) 
    false_positives <- length(calls_Granges_file)-length(matches_Granges_file)
    false_negatives <- length(truth_Granges_file)-length(matches_Granges_file)
    precision <- true_positives / (true_positives + false_positives)
    recall <- true_positives / (true_positives + false_negatives)
    return(list(
        true_positives = true_positives, 
        false_positives = false_positives,
        false_negatives = false_negatives,
        precision = precision,
        recall = recall))

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

read_bed_file <- function(file_path) {
  # Extract the base name without the suffix
  base_name <- sub("_Ref_Coord.bed$", "", basename(file_path))
  # Read the BED file into a data frame
  bed_df <- read.table(file_path, header = FALSE)
  # Return a list with the name and the data frame
  list(name = base_name, file_path = file_path, df = bed_df)
}

ISO1_og <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/ISO-1_Ref_Coord.bed") #%>%
    #filter(V1 %in% c("2L", "2R"))

euchromatin_coordinates <- makeGRangesFromDataFrame(read.table(euchromatin_coordinates_path), seqnames.field="V1", start.field="V2", end.field="V3")

#ISO1_filtered <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/DeNovoCoordinates/ISO1.bed")
ISO1_og_gr <- GRanges(
    seqnames = ISO1_og$V1,
    ranges = IRanges(start = ISO1_og$V2, end = ISO1_og$V3),
    mcols=ISO1_og$V7
    )

#specify genome we are analyzing. 
#genome <- "COR-023"

#benchmark <- function(genome1,genome2) {
genome <- paste0(genome1, "_", genome2)
#result_path <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage,  "X/", genome, "_1/results/")
#result_path <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage, "X/", genome, "_1/results/")
result_path <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage, "X/", genome, "_1/results/")


nonreference_genome1_truth_path <- paste0("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/", genome1, "_Ref_Coord.bed")
nonreference_genome2_truth_path <- paste0("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/", genome2, "_Ref_Coord.bed")



#moving files
source_file <- paste0(basedir_outputs_path, "/3L3RX_model_validation_output/", genome1, "_", genome2, "_TEforest_nonredundant.bed")
# Specify the destination directory
destination_folder <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage,  "X/", genome, "_1/results/", caller_name)
if (!dir.exists(destination_folder)) {
  # Create the directory
  dir.create(destination_folder)
}

destination_file <- paste0(destination_folder, "/", genome, "_1_", caller_name, "_nonredundant.bed")
# Use file.copy() to copy the file
file.copy(from = source_file, to = destination_file, overwrite=T)

source_file <- paste0(basedir_outputs_path, "/3L3RX_model_validation_output/", genome1, "_", genome2, "_TEforest_bps_nonredundant.bed")
# Specify the destination directory
destination_folder <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage,  "X/", genome, "_1/results/", caller_name, "_bps")
if (!dir.exists(destination_folder)) {
  # Create the directory
  dir.create(destination_folder)
}

destination_file <- paste0(destination_folder, "/", genome, "_1_", caller_name, "_bps_nonredundant.bed")
# Use file.copy() to copy the file
file.copy(from = source_file, to = destination_file, overwrite=T)

# now repeat for classifier stuff
source_file <- paste0(basedir_outputs_path, "/3L3RX_model_validation_output_classifier/", genome1, "_", genome2, "_TEforest_nonredundant.bed")
# Specify the destination directory
destination_folder <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage,  "X/", genome, "_1/results/", caller_name2)
if (!dir.exists(destination_folder)) {
  # Create the directory
  dir.create(destination_folder)
}

destination_file <- paste0(destination_folder, "/", genome, "_1_", caller_name, "_nonredundant.bed")
# Use file.copy() to copy the file
file.copy(from = source_file, to = destination_file, overwrite=T)

source_file <- paste0(basedir_outputs_path, "/3L3RX_model_validation_output_classifier/", genome1, "_", genome2, "_TEforest_bps_nonredundant.bed")
# Specify the destination directory
destination_folder <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/coverage_tests/", coverage,  "X/", genome, "_1/results/", caller_name2, "_bps")
if (!dir.exists(destination_folder)) {
  # Create the directory
  dir.create(destination_folder)
}

destination_file <- paste0(destination_folder, "/", genome, "_1_", caller_name, "_bps_nonredundant.bed")
# Use file.copy() to copy the file
file.copy(from = source_file, to = destination_file, overwrite=T)



#creating heterozygote truth
#complex sequence to turn tp's on chroms 2r and 3r to hets, and rest of the genome just genome1 calls
nonreference_genome1_truth <- read.table(nonreference_genome1_truth_path) 
colnames(nonreference_genome1_truth) <- c("seqnames", "start", "end", "ID", "length", "strand", "TE")
nonreference_genome1_truth$TE <- gsub("-", "_", nonreference_genome1_truth$TE)

nonreference_genome2_truth <- read.table(paste0("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/", genome2, "_Ref_Coord.bed"))
colnames(nonreference_genome2_truth) <- c("seqnames", "start", "end", "ID", "length", "strand", "TE")
nonreference_genome2_truth$TE <- gsub("-", "_", nonreference_genome2_truth$TE)
#keep only chromosomes 2R and 3R
nonreference_genome2_truth <- nonreference_genome2_truth  %>% filter(seqnames %in% c("2R","3R"))
nonreference_genome1_truth_r <- nonreference_genome1_truth  %>% filter(seqnames %in% c("2R","3R"))
nonreference_genome1_truth_notr <- nonreference_genome1_truth %>% filter(seqnames %in% c("2L","3L", "X"))

grange1 <- makeGRangesFromDataFrame(nonreference_genome1_truth, keep.extra.columns = T)
grange2 <- makeGRangesFromDataFrame(nonreference_genome2_truth, keep.extra.columns = T)
grange1_r <- makeGRangesFromDataFrame(nonreference_genome1_truth_r, keep.extra.columns = T)
grange1_notr <- makeGRangesFromDataFrame(nonreference_genome1_truth_notr, keep.extra.columns = T)

combined_2R3R <- unique(as.data.frame(c(grange1_r, grange2)))
TEs <- unique(sort(c(ISO1_og$V7, grange1$TE, grange2$TE)))
het_df <- tibble(TE=TEs)

# Function to check for exact match and add 0.5 to a new column
addMatchColumn <- function(query, reference) {
  is_exact_match <- query %in% reference
  ref_values <- ifelse(is_exact_match, 0.5, 0)
  return(ref_values)
}
check_heterozygotes <- function(grange) {
    # Check for exact matches with grange1 and grange2 for query_gr
    grange$grange1 <- addMatchColumn(grange, grange1_r)
    grange$grange2 <- addMatchColumn(grange, grange2)
    grange$heterozygosity <- grange$grange1 + grange$grange2
    return(grange)
}

heterozygote_truth_creation <- het_df %>% mutate(
    nonreference_2R3R_truth_forTE = map(TE, 
        ~ makeGRangesFromDataFrame(combined_2R3R %>% filter(TE==.x), keep.extra.columns = T )),
    add_heterozygote_cols = map(nonreference_2R3R_truth_forTE, 
        ~ check_heterozygotes(.x))
)

heterozygosity_2R3R <- as.data.frame(do.call(c, c(heterozygote_truth_creation$add_heterozygote_cols)))
grange1_notr_heterozygosity <- nonreference_genome1_truth_notr
grange1_notr_heterozygosity$grange1 <- 1.0
grange1_notr_heterozygosity$grange2 <- 0.0
grange1_notr_heterozygosity$heterozygosity <- 1.0
heterozygosity_2R3R <- heterozygosity_2R3R %>% select(seqnames, start, end, ID, length, strand, TE, grange1, grange2, heterozygosity)
truth <- rbind(heterozygosity_2R3R, grange1_notr_heterozygosity)


ISO1_ids <- ISO1_og$V4

# Move truth dataframes to mcclintock output directories
# Create a data frame of data frames
genome_rename <- paste0(genome, "_1")
mcclintock_format_truth <- truth %>% mutate(
    start = start - 1, 
    ref = if_else(ID %in% ISO1_ids, "reference", "non-reference"),
    TE_string = paste(TE, ref, heterozygosity, genome_rename, "truth", "rp", "1", sep="|"),
    score = 0, 
    strand = "."
) %>% 
  select(seqnames, start, end, TE_string, score, strand) %>%
  arrange(as.character(seqnames),start)


destination_folder <- paste0(result_path, "truth")
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

write.table(mcclintock_format_truth, 
    file = paste0(destination_folder, "/", genome, "_1_", "truth", "_nonredundant.bed"),
    quote = F, sep = "\t", row.names = F
)



#filtering to just validation chromosomes
truth <- truth %>% filter(seqnames %in% c("2L","2R"))
#truth <- truth %>% filter(seqnames %in% c("3L","3R", "X"))

#splitting truth based on het or homozygote to get total counts
truth_het <- truth %>% filter(heterozygosity==0.5)
truth_homo <- truth %>% filter(heterozygosity==1)
truth_het_gr <- GRanges(
    seqnames = truth_het$seqnames,
    ranges = IRanges(start = truth_het$start, end = truth_het$end)
    )
truth_het_reference_gr <- subsetByOverlaps(truth_het_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)

truth_homo_gr <- GRanges(
    seqnames = truth_homo$seqnames,
    ranges = IRanges(start = truth_homo$start, end = truth_homo$end)
    )
truth_homo_reference_gr <- subsetByOverlaps(truth_homo_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)
#defining a TP as within 500 bp of true site
#right now i am not checking for nests/overlaps or accurate calls of the wrong TE class
truth_gr <- GRanges(
    seqnames = truth$seqnames,
    ranges = IRanges(start = truth$start, end = truth$end)
    )
truth_reference_gr <- subsetByOverlaps(truth_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)

truth_nonreference_extend_gr <- extend(subsetByOverlaps(truth_gr, ISO1_og_gr, type=c("equal"), invert=TRUE, ignore.strand=TRUE), 500, 500)
num_true_ref <- length(truth_reference_gr)
num_true_nonref <- length(truth_nonreference_extend_gr)

# find a random sample of true negatives to add to confusion matrix
# equal to the number of tp nonreference calls
nrow(training_csv)

training_csv_negative <- training_csv %>% filter(true==0)  
genome_name_length <- length(strsplit(genome, "-")[[1]])

if (genome_name_length == 1) {
    names <- c("genome", "seqnames", "start", "end", "class", "TE")
    training_csv_negative_downsampled <- training_csv_negative %>% 
    separate_wider_delim(file,names = names, delim = "-") %>% 
  sample_n(min(nrow(training_csv_negative),num_true_nonref))
} else if (genome_name_length == 2) {
    names <- c("genome1", "genome2", "seqnames", "start", "end", "class", "TE")
    training_csv_negative_downsampled <- training_csv_negative %>% 
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    sample_n(min(nrow(training_csv_negative),num_true_nonref))
} else if (genome_name_length == 3) {
    names <- c("genome1", "genome2", "genome3", "seqnames", "start", "end", "class", "TE")
    training_csv_negative_downsampled <- training_csv_negative %>% 
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  genome3, sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    sample_n(min(nrow(training_csv_negative),num_true_nonref))
} else {
  names <- "Error"
}

mcclintock_results <- tibble(
  caller = list.files(result_path)[!grepl("summary", list.files(result_path))],
  result_file = unlist(map(paste0(result_path, caller), ~list.files(.x, pattern = "_nonredundant.bed"))),
  fullpath = paste0(result_path, caller, "/", result_file),
  data = lapply(fullpath, read_mcclintock_het_format),
  df_list = map(data, ~split(.x, .x[["X4"]])), #split dataframes based on individual TEs
  genome = str_extract(fullpath, "(?<=/)[A-Z0-9-]+_[A-Z0-9-]+(?=_1/)"))%>%
  unnest(cols=df_list) %>% 
  mutate(TE = unlist(map(df_list, ~.x$X4[1])),
    length = unlist(map(df_list, ~length(.x$X1))))
 #use df_list to find TE name

#mcclintock_results <- mcclintock_results %>% filter(TE=="jockey")

zero_calls_df <- subset(mcclintock_results, length==0)
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
    TE = .x$X4, heterozygosity= .x$X6)), 
  nonref_gr_filter = map(nonref_gr, ~ subsetByOverlaps(.x, euchromatin_coordinates, ignore.strand = TRUE)),
  # I gave mcclintock a 1-based bed file for refs, so don't need to correct this
  # This means the programs were checking with a start +1 of where it needed to be
  reference_freq_filter = map(reference, ~ .x %>% filter(X6 >= 0.10)),
  ref_gr = map(reference_freq_filter, ~ GRanges(
    seqnames = .x$X1,
    #teflon ref coords are off by 1 for unknown reasons (mcclintock error?). Adding 1 to correct this. 
    ranges = IRanges(start = ifelse(.x$X8 == "teflon", as.numeric(.x$X2)+1, as.numeric(.x$X2)), end = as.numeric(.x$X3)),
    TE = .x$X4, heterozygosity = ifelse(.x$X6 == "NA", '1', .x$X6))),
  nonreference_freq_filter = map(nonref_gr_filter, ~ .x %>% as.data.frame() %>% filter(heterozygosity >= 0.10) %>% makeGRangesFromDataFrame(keep.extra.columns=T)),
)

A1_truth <- truth[c(1:7,10)]
colnames(A1_truth) <- c("seqnames", "start", "end", "ID", "length", "strand", "TE", "heterozygosity")
A1_truth$TE <- gsub("-", "_", A1_truth$TE) #critical
TE_list <- unique(sort(A1_truth$TE))
A1_truth$start <- A1_truth$start #do not need to make 1 based (they already are!)
A1_truth <- A1_truth

gr <- GRanges(seqnames = A1_truth$seqnames,
              ranges = IRanges(start = A1_truth$start, end = A1_truth$end),
              TE = A1_truth$TE,
              ID = A1_truth$ID)

# Find overlaps
overlap_hits <- findOverlaps(gr, gr, ignore.strand = TRUE)

# Initialize a column for the overlap category
A1_truth$overlap_category <- "no_overlap"

# Process the overlaps
for (i in seq_along(gr)) {
  overlaps <- subjectHits(overlap_hits[queryHits(overlap_hits) == i])
  if (length(overlaps) > 1) { # Check if there are overlapping rows
    overlap_TEs <- unique(gr[overlaps]$TE)
    if (length(overlap_TEs) == 1) {
      A1_truth$overlap_category[i] <- "nested_same_TE"
    } else {
      A1_truth$overlap_category[i] <- "nested_different_TE"
    }
  }
}


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
    f1_score = unlist(map2(nonreference_freq_filter, A1_truth_nonref, ~ f1_score(.x, .y))),
    nonref_false_negatives = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE, invert=TRUE)),
    nonref_false_negatives = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE, invert=TRUE)),
    nonref_true_positives = map2(nonreference_freq_filter, A1_truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)),
    nonreference_freq_filter_extend = map(nonreference_freq_filter, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    nonref_true_positives_truth_dataset_coords = map2(A1_truth_nonref_noextend, nonreference_freq_filter_extend, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)), #find tp coords in truth dataset ranges
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

path <- "/nas/longleaf/home/adaigle/TEforest/workflow/scripts/"
save(benchmark_mapping_results, file = paste0(path, "benchmark_mapping_results_150_50X.RData"))

benchmark_mapping_results_ref <- mcclintock_results %>% mutate(
    A1_truth_forTE = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x), keep.extra.columns = T )),
    A1_truth_ref = map(A1_truth_forTE,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type="equal", ignore.strand=TRUE)),
    #nonref_gr_filter = map2(nonref_gr_filter, TE, ~ GRanges(seqnames = seqnames(.x), 
    #    ranges = ranges(.x),
    #    TE = .y, heterozygosity=nonref_gr_filter$heterozygosity)),
    #f1_score = unlist(map2(nonreference_freq_filter, A1_truth_nonref, ~ f1_score(.x, .y))),
    ref_false_negatives = map2(ref_gr, A1_truth_ref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE, invert=TRUE, type="equal")),
    ref_true_positives = map2(ref_gr, A1_truth_ref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, type="equal")),
    ref_false_positives = map2(ref_gr, A1_truth_ref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, invert = TRUE, type="equal")),
    ref_true_positives_length = unlist(map(A1_truth_ref,~ length(.x))),
    ref_calls_true_positives_length = unlist(map(ref_true_positives,~ length(.x))),
    ref_calls_false_positives_length = unlist(map(ref_false_positives,~ length(.x))),
    ref_calls_false_negatives_length = unlist(map(ref_false_negatives,~ length(.x))),
    ref_calls_length= unlist(map(ref_gr,~ length(.x))),
    #stats = map2(ref_gr_filter, A1_truth_ref, ~ tp_fp_fn(.x, .y))
    precision = ref_calls_true_positives_length / (ref_calls_true_positives_length + ref_calls_false_positives_length),
    recall = ref_calls_true_positives_length / (ref_calls_true_positives_length + ref_calls_false_negatives_length),
    f1_score = 2 * (precision * recall) / (precision + recall)
)
#benchmark_mapping_results_newcaller <- benchmark_mapping_results %>% filter(caller=='TEforest_classifier_filter')
#benchmark_mapping_results_oldcaller <- benchmark_mapping_results %>% filter(caller=='TEforest_classifier_old')
#
#newnof <- benchmark_mapping_results_newcaller |> filter(TE=="NOF")

#benchmark_mapping_results_ref_subset <- benchmark_mapping_results_ref %>% filter(caller=='teflon', TE=="roo")
#benchmark_mapping_results_ref_subset <- benchmark_mapping_results_ref %>% filter(caller=='TEforest_classifier_bps')
#benchmark_mapping_results_ref_subset <- benchmark_mapping_results_ref %>% filter(caller=='temp2', TE=="roo")


benchmark_mapping_results2 <- benchmark_mapping_results %>%
  mutate(
    distances_vector = map2(nonref_true_positives, A1_truth_nonref_noextend, ~ distanceToNearest(center(.x), .y, ignore.strand = TRUE)),
    distances_vector2 = map(distances_vector, ~ .x@elementMetadata[[1]]),
    distances_vector_mean = unlist(map(distances_vector, ~ mean(.x@elementMetadata[[1]]))),
    distances_vector_sd = unlist(map(distances_vector, ~ sd(.x@elementMetadata[[1]])))
  )

#of my true positives, find which ones are homozygotes called as heterozygotes, etc
benchmark_mapping_results_het <- benchmark_mapping_results2 %>%
  mutate(
    A1_truth_forTE_het = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x, heterozygosity >= 0.25, heterozygosity <= 0.75),  keep.extra.columns = T )),
    A1_truth_ref_het = map(A1_truth_forTE_het,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref_het = map(A1_truth_forTE_het, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    A1_truth_forTE_homo = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x, , heterozygosity >= 0.75),  keep.extra.columns = T )),
    A1_truth_ref_homo = map(A1_truth_forTE_homo,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref_homo = map(A1_truth_forTE_homo, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    nonref_true_positives_homo = map2(nonref_true_positives, A1_truth_nonref_homo, 
        ~ as.data.frame(subsetByOverlaps(.x, .y, ignore.strand = TRUE))),
    homo_tp = map(nonref_true_positives_homo, ~ .x %>% filter(heterozygosity >= 0.75)),
    homo_fp_het = map(nonref_true_positives_homo, ~ .x %>% filter(heterozygosity >= 0.25, heterozygosity <= 0.75)),
    homo_tp_length = unlist(map(homo_tp, ~ nrow(.x))),
    homo_fp_het_length = unlist(map(homo_fp_het, ~ nrow(.x))),
    homo_fn = map(nonref_false_negatives, ~ .x %>% as.data.frame %>% filter(heterozygosity==1.0)),
    homo_fn_length = unlist(map(homo_fn, ~ sum(.x$heterozygosity==1.0) )),
    nonref_true_positives_het = map2(nonref_true_positives, A1_truth_nonref_het, 
        ~ as.data.frame(subsetByOverlaps(.x, .y, ignore.strand = TRUE))),
    het_tp = map(nonref_true_positives_het, ~ .x %>% filter(heterozygosity >= 0.25, heterozygosity <= 0.75)),
    #eventually need to differentiate FPs that are below 0.25, as these are a different class of error
    het_fp_homo = map(nonref_true_positives_het, ~ .x %>% filter(heterozygosity > 0.75)),
    het_tp_length = unlist(map(het_tp, ~ nrow(.x))),
    het_fp_homo_length = unlist(map(het_fp_homo, ~ nrow(.x))),
    het_fn = map(nonref_false_negatives, ~ .x %>% as.data.frame %>% filter(heterozygosity==0.5)),
    het_fn_length = unlist(map(het_fn, ~ sum(.x$heterozygosity==0.5))),
    homo_fp_zero = map(nonref_true_positives_homo, ~ .x %>% filter(heterozygosity < 0.25)),
    het_fp_zero = map(nonref_true_positives_het, ~ .x %>% filter(heterozygosity < 0.25)),
    homo_fp_zero_length = unlist(map(homo_fp_zero, ~ nrow(.x))),
    het_fp_zero_length = unlist(map(het_fp_zero, ~ nrow(.x))),
    fp = map(nonref_false_positives, ~ .x %>% as.data.frame),
    homo_fp_length = unlist(map(fp, ~ sum(.x$heterozygosity>=0.75))),
    het_fp_length = unlist(map(fp, ~ sum(.x$heterozygosity<0.75))),
  )


benchmark_mapping_results_het <- benchmark_mapping_results_het %>% mutate(
  true_negatives_truth = map(
      TE, 
      ~ {
        df <- training_csv_negative_downsampled %>% filter(TE == .x) 
        if (nrow(df) == 0) {
          return(GRanges())  # Return empty grange if no true negatives for this TE fam
        } else {
          return(makeGRangesFromDataFrame(df, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE))
        }
      }
    ),
  true_negatives = map2(true_negatives_truth, nonreference_freq_filter,
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, invert=TRUE)),
  true_negatives_length = unlist(map(true_negatives, ~ length(.x))),
)
#test <- benchmark_mapping_results_het %>% filter(caller=="TEforest_classifier")
#test2 <- benchmark_mapping_results_het %>% filter(caller=="2L_test")

benchmark_mapping_results_frequency_plot <- benchmark_mapping_results_het %>%
  mutate(
    homo_tp_labeled = map2(caller, homo_tp, ~ .y %>% mutate(caller = .x, truth = 1.0)),
    het_tp_labeled = map2(caller, het_tp, ~ .y %>% mutate(caller = .x, truth = 0.5)),
    homo_fp_het_labeled = map2(caller, homo_fp_het, ~ .y %>% mutate(caller = .x, truth = 1.0)),
    het_fp_homo_labeled = map2(caller, het_fp_homo, ~ .y %>% mutate(caller = .x, truth = 0.5)),
    homo_fn_labeled = map2(caller, homo_fn, ~ .y %>% mutate(caller = .x, truth = 1.0, ID = NULL, heterozygosity=0, length = NULL)),
    het_fn_labeled = map2(caller, het_fn, ~ .y %>% mutate(caller = .x, truth = 0.5, heterozygosity=0, ID = NULL, length = NULL)),
    homo_fp_zero_labeled = map2(caller, homo_fp_zero, ~ .y %>% mutate(caller = .x, truth = 1.0)),
    het_fp_zero_labeled = map2(caller, het_fp_zero, ~ .y %>% mutate(caller = .x, truth = 0.5)),
    fp = map2(caller, fp, ~ .y %>% mutate(caller = .x, truth = 0)),
  )

benchmark_mapping_results_frequency_plot <- benchmark_mapping_results_frequency_plot %>% mutate(
  tn_df = map(true_negatives, # reduce to get rid of nested TEs of same type
        ~ as.data.frame(.x)),
  tn = map2(caller, tn_df, ~ .y %>% mutate(caller = .x, truth = 0, heterozygosity=0, true = NULL, pred=NULL, cntrl_score=NULL, genome=NULL, class=NULL)))

# Combine the labeled columns into a single data frame
frequency_plot <- do.call(rbind, c(
  benchmark_mapping_results_frequency_plot$homo_tp_labeled,
  benchmark_mapping_results_frequency_plot$het_tp_labeled,
  benchmark_mapping_results_frequency_plot$homo_fp_het_labeled,
  benchmark_mapping_results_frequency_plot$het_fp_homo_labeled,
  benchmark_mapping_results_frequency_plot$homo_fn_labeled,
  benchmark_mapping_results_frequency_plot$het_fn_labeled,
  benchmark_mapping_results_frequency_plot$homo_fp_zero_labeled,
  benchmark_mapping_results_frequency_plot$het_fp_zero_labeled,
  benchmark_mapping_results_frequency_plot$fp,
  benchmark_mapping_results_frequency_plot$tn
))


# Calculate R-squared value
calculate_r_squared <- function(data) {
  lm_model <- lm(truth ~ as.numeric(heterozygosity), data = data)
  return(data.frame(R_squared = summary(lm_model)$r.squared))
}
r_squared_results <- frequency_plot %>%
  group_by(caller) %>%
  do(calculate_r_squared(.))

frequency_plot <- left_join(frequency_plot, r_squared_results, by = "caller") %>% mutate(heterozygosity = as.numeric(heterozygosity))

# Create scatter plot with separate linear regression lines for each group
#ggplot(frequency_plot, aes(x = truth, y = heterozygosity, color = caller)) +
#  geom_point() +
#  geom_smooth(method = "lm", se = FALSE, aes(group = caller), color = "blue") +
#  labs(title = "Truth vs Prediction Scatter Plot with Separate Regression Lines",
#       x = "Truth", y = "Prediction") +
#  facet_wrap(~caller) +
#  theme_minimal()

freqplt <- ggplot(frequency_plot, aes(x = truth, y = heterozygosity, color = caller)) +
  geom_jitter(alpha = 0.2, width = 0.1) +
  #geom_boxplot(aes(group=truth), width = 0.05, fill = "lightgray", outlier.shape = NA) +
  geom_violin(aes(group=truth)) + 
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~caller) +
  labs(title = "Scatterplot with Boxplot of True vs Predicted Values",
       x = "True Values",
       y = "Predicted Values") 

frequency_plot_nozero <- frequency_plot %>% filter(heterozygosity!=0)
freqpltnozero <-ggplot(frequency_plot_nozero, aes(x = truth, y = heterozygosity, color = caller)) +
  geom_jitter(alpha = 0.2, width = 0.1) +
  geom_boxplot(aes(group=truth), width = 0.05, fill = "lightgray", outlier.shape = NA) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~caller) +
  labs(title = "Scatterplot with Boxplot of True vs Predicted Values",
       x = "True Values",
       y = "Predicted Values") 
#check using roo
#tmp <- benchmark_mapping_results2 %>% select(caller, TE, nonref_calls_true_positives_length, distances_vector_mean, distances_vector_sd)
#as.data.frame(tmp %>% filter(nonref_calls_true_positives_length!=0) %>% filter(TE=="roo"))


benchmark_mapping_results_nodata <- benchmark_mapping_results_het[c(1,7,8,20,24:30,33, 45:48, 52:53, 55, 58:59, 61:62, 65)]

# Check for NaN values and replace with 0
benchmark_mapping_results_nodata[is.na(benchmark_mapping_results_nodata)] <- 0


#length of reduced nonref DF
#will mess up reduced nested TEs of different classes
#but necessary for proper analysis 
nonref_te_number <- length(unique(truth_nonreference_extend_gr))
truth_het_nonreference_gr <- subsetByOverlaps(truth_het_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert=TRUE)
nonref_te_number_het <- length(unique(truth_het_nonreference_gr))

truth_homo_nonreference_gr <- subsetByOverlaps(truth_homo_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert=TRUE)
nonref_te_number_homo <- length(unique(truth_homo_nonreference_gr))

benchmark_mapping_results_summary <- benchmark_mapping_results_nodata %>%
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
    distance_mean = mean(unlist(distances_vector2)),
    distance_sd = sd(unlist(distances_vector2)),
    distance_vectors = list(c(distances_vector2)),
    homo_tp_length = sum(homo_tp_length),
    het_tp_length = sum(het_tp_length),
    homo_fp_het_length = sum(homo_fp_het_length),
    het_fp_homo_length = sum(het_fp_homo_length),
    #homo_fn_length = sum(homo_fn_length),
    homo_fn_length = nonref_te_number_homo - homo_tp_length - homo_fp_het_length,
    #het_fn_length = sum(het_fn_length),
    het_fn_length = nonref_te_number_het - het_tp_length - het_fp_homo_length,
    homo_fp_zero_length = sum(homo_fp_zero_length),
    het_fp_zero_length = sum(het_fp_zero_length),
    homo_fp_length = sum(homo_fp_length),
    het_fp_length = sum(het_fp_length),
    true_negatives_length = sum(true_negatives_length),
  )
benchmark_mapping_results_summary$recalc_false_negatives <- ifelse(benchmark_mapping_results_summary$recalc_false_negatives < 0, 0, benchmark_mapping_results_summary$recalc_false_negatives)
benchmark_mapping_results_summary$homo_fn_length <- ifelse(benchmark_mapping_results_summary$homo_fn_length < 0, 0, benchmark_mapping_results_summary$homo_fn_length)
benchmark_mapping_results_summary$het_fn_length <- ifelse(benchmark_mapping_results_summary$het_fn_length < 0, 0, benchmark_mapping_results_summary$het_fn_length)

benchmark_mapping_results_summary <- benchmark_mapping_results_summary %>% mutate(
    precision = sum_true_positives / (sum_true_positives + sum_false_positives),
    recall = sum_true_positives / (sum_true_positives + recalc_false_negatives),
    f1_score = 2 * (precision * recall) / (precision + recall)
)

benchmark_mapping_results_summary[is.na(benchmark_mapping_results_summary)] <- 0

benchmark_mapping_results_summary <- benchmark_mapping_results_summary %>%
  mutate(distance_vectors = map(distance_vectors, ~ unlist(unname(.x))),
  distance_vectors_df = map(distance_vectors, ~ data.frame(distance = .x)))

het_conf_matrix <- function(caller_name) {
num <- as.numeric(row.names(benchmark_mapping_results_summary)[benchmark_mapping_results_summary$caller == caller_name])

conf_matrix_mycaller <- benchmark_mapping_results_summary[num,]
#heterozygosity classification accuracy
conf_matrix_values <- matrix(c(0, conf_matrix_mycaller$het_fn_length, conf_matrix_mycaller$homo_fn_length + conf_matrix_mycaller$homo_fp_zero_length, # need to grab fp hets and homos
                                conf_matrix_mycaller$het_fp_length, conf_matrix_mycaller$het_tp_length, conf_matrix_mycaller$homo_fp_het_length, 
                                conf_matrix_mycaller$homo_fp_length, conf_matrix_mycaller$het_fp_homo_length, conf_matrix_mycaller$homo_tp_length),
                             nrow = 3, byrow = TRUE)
conf_matrix <- confusionMatrix(conf_matrix_values)
cm <- conf_mat(conf_matrix_values)

plt <- autoplot(cm, type = "heatmap") +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1") + 
  theme(legend.position = "right")+
  scale_x_discrete(labels = c("absent", "heterozygote", "homozygote")) +
  scale_y_discrete(labels = c("homozygote", "heterozygote", "absent"))

#summary(cm)[1,]$.estimate
return(plt)
}
#conf_matrix_mycaller <- het_conf_matrix("2L_test")
#het_conf_matrix("temp")
#et_conf_matrix("temp2")
#et_conf_matrix("teflon")
#et_conf_matrix("retroseq")
#et_conf_matrix("popoolationte")
#et_conf_matrix("popoolationte2")



het_conf_matrix <- benchmark_mapping_results_summary %>%
  rowwise() %>%
  mutate(
    conf_matrix = list(matrix(c(nrow(training_csv_negative_downsampled) - het_fp_length - homo_fp_length, het_fn_length, homo_fn_length,
                                 het_fp_length, het_tp_length, homo_fp_het_length,
                                 homo_fp_length, het_fp_homo_length, homo_tp_length),
                              nrow = 3, byrow = TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    conf_matrix_summary = map(conf_matrix, ~confusionMatrix(.x)),
    het_stats = map(conf_matrix_summary, ~ c(precision=.x$byClass["Class: B","Precision"], 
      recall=.x$byClass["Class: B","Recall"], f1=.x$byClass["Class: B","F1"])),
    homo_stats = map(conf_matrix_summary, ~ c(precision=.x$byClass["Class: C","Precision"], 
      recall=.x$byClass["Class: C","Recall"], f1=.x$byClass["Class: C","F1"])),
  )

generate_confusion_matrix <- function(het_conf_matrix, caller_name) {
  # Extract the confusion matrix for the specified caller
  cm <- het_conf_matrix %>% filter(caller == caller_name) %>% pull(conf_matrix)
  if (length(cm) == 0) {
    # If caller is not present, return NA and a message or create a blank plot with a warning text
    warning(sprintf("Caller '%s' not found in the dataframe. Returning NA.", caller_name))
    return(NA)
  }
  cm <- cm[[1]]

  # Calculate column-wise proportions for coloring
  cm_proportion <- sweep(cm, 2, colSums(cm), FUN = "/")

  # Convert the original matrix to a data frame
  cm_df <- as.data.frame(cm)
  cm_df$Actual <- factor(c("0", "0.5", "1"))

  # Convert the proportion matrix to a data frame
  cm_proportion_df <- as.data.frame(cm_proportion)
  cm_proportion_df$Actual <- factor(c("0", "0.5", "1"))

  # Rename the columns to reflect Predicted classes
  colnames(cm_df) <- c("0", "0.5", "1", "Actual")
  colnames(cm_proportion_df) <- c("0", "0.5", "1", "Actual")

  # Reshape the data frames to long format
  cm_long <- cm_df %>% pivot_longer(cols = -Actual, names_to = "Predicted", values_to = "Value")
  cm_proportion_long <- cm_proportion_df %>% pivot_longer(cols = -Actual, names_to = "Predicted", values_to = "Proportion")

  # Combine the data frames
  cm_combined <- cm_long %>%
    mutate(Proportion = cm_proportion_long$Proportion, 
           Actual = factor(Actual, levels = c(0, 0.5, 1), labels = c("absent", "heterozygote", "homozygote")),
           Predicted = factor(Predicted, levels = c(0, 0.5, 1), labels = c("absent", "heterozygote", "homozygote"))
    )

  # Generate the confusion matrix plot
  confusion_matrix <- ggplot(cm_combined, aes(x = Predicted, y = Actual, fill = Proportion)) +
    geom_tile(color = "grey", size = 0.5) +
    scale_fill_gradient(low = "white", high = "#0066ff") +
    geom_text(aes(label = sprintf("%.2f", Proportion)), size = 8) +
    labs(x = "True Genotype", y = "Predicted Genotype", fill = "Col proportion") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(confusion_matrix)
}

confusion_matrix_TEforest <- generate_confusion_matrix(het_conf_matrix, "TEforest_classifier_filter_bps")
confusion_matrix_retroseq <- generate_confusion_matrix(het_conf_matrix, "retroseq")
confusion_matrix_temp2 <- generate_confusion_matrix(het_conf_matrix, "temp2")

# Expand the list columns into separate columns
het_stats_plot <- het_conf_matrix %>% select(caller, het_stats) %>%
  tidyr::unnest_wider(het_stats) %>%
  gather(metric, value, precision:f1) %>%
  group_by(caller) 
# Calculate the average F1 score for each caller and order by this score
f1_order <- het_stats_plot %>%
  filter(metric == "f1") %>%
  group_by(caller) %>%
  summarise(mean_f1 = mean(value)) %>%
  arrange(desc(mean_f1)) %>%
  pull(caller)

# Reorder the caller factor levels based on F1 score
het_stats_plot <- het_stats_plot %>%
  mutate(caller = factor(caller, levels = f1_order))

# Create the bar plot
hetplt <- ggplot(het_stats_plot, aes(x = caller, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Caller", y = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),  # Adjust legend title size
        legend.text = element_text(size = 12),   # Adjust legend text size
        axis.text = element_text(size = 12),     # Adjust axis text size
        axis.title = element_text(size = 14))  +  # Adjust axis title size
  expand_limits(y=c(0,1)) +
# Add F1 score numbers above the bars
 geom_text(aes(label = sprintf("%.2f", value), y = value, group = metric),
              position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)


# Expand the list columns into separate columns
homo_stats_plot <- het_conf_matrix %>% select(caller, homo_stats) %>%
  tidyr::unnest_wider(homo_stats) %>%
  gather(metric, value, precision:f1) %>%
  group_by(caller)

f1_order_homo <- homo_stats_plot %>%
  filter(metric == "f1") %>%
  group_by(caller) %>%
  summarise(mean_f1 = mean(value)) %>%
  arrange(desc(mean_f1)) %>%
  pull(caller)

# Reorder the caller factor levels based on F1 score
homo_stats_plot <- homo_stats_plot %>%
  mutate(caller = factor(caller, levels = f1_order_homo))


# Create the bar plot
homoplt <- ggplot(homo_stats_plot, aes(x = caller, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Caller", y = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),  # Adjust legend title size
        legend.text = element_text(size = 12),   # Adjust legend text size
        axis.text = element_text(size = 12),     # Adjust axis text size
        axis.title = element_text(size = 14))  +  # Adjust axis title size
  expand_limits(y=c(0,1)) +
# Add F1 score numbers above the bars
 geom_text(aes(label = sprintf("%.2f", value), y = value, group = metric),
              position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

# Create the ggplot
benchmark_mapping_results_summary_plot <- benchmark_mapping_results_summary #%>% #filter(caller !="ngs_te_mapper")%>% filter(caller !="ngs_te_mapper2") %>%
    #filter(caller !=c("allseeingeye", "allseeingeye2", "allseeingeye3"))
benchmark_mapping_results_summary_plot <- benchmark_mapping_results_summary_plot #%>% 
  #mutate(caller = case_when(
  #  caller == "allseeingeye4" ~ "mystuff",
  #  caller == "allseeingeye5" ~ "mystuff_with_bps",
  #  TRUE ~ caller  # Keep other values unchanged
  #))
benchmark_mapping_results_summary_plot <- benchmark_mapping_results_summary_plot %>%
  arrange(desc(f1_score))



# Reshape the data to long format
df_long <- pivot_longer(benchmark_mapping_results_summary_plot, cols = c(precision, recall, f1_score), names_to = "Metric", values_to = "Score")

df_long$caller <- factor(df_long$caller, levels = benchmark_mapping_results_summary_plot$caller)
df_long_save <- df_long
df_long_save$covg <- paste0(coverage, "X")
# Create a grouped bar plot
p <- ggplot(df_long, aes(x = caller, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  #scale_fill_manual(values = c("precision" = "blue", "recall" = "green", "f1_score" = "red")) +
  labs(y = "Score",
       x = "Caller") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),  # Adjust legend title size
        legend.text = element_text(size = 12),   # Adjust legend text size
        axis.text = element_text(size = 12),     # Adjust axis text size
        axis.title = element_text(size = 14))  +  # Adjust axis title size
  expand_limits(y=c(0,1)) +
# Add F1 score numbers above the bars
 geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = Metric),
              position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)#return(benchmark_mapping_results_summary)


#breakpoints benchmarking

make_bp_plot <- function(caller_name, benchmark_mapping_results_summary) {
#if (caller_name == "TEforest_classifier_filter_bps") {
#  caller_name <- "TEforest"
#}

colors <- c(
    "TEforest_classifier_filter_bps" = "#F8766D",
    "temp2" = "#CD9600",
    "retroseq" = "#7CAE00",
    "temp" = "#00BE67",
    "popoolationte" = "#00BFC4",
    "teflon" = "#00A9FF",
    "PopoolationTE2" = "#C77CFF",
    "tepid" = "#FF61CC",
    "TEforest_nolengthfilter" = "black"
)

num <- as.numeric(row.names(benchmark_mapping_results_summary)[benchmark_mapping_results_summary$caller == caller_name])
# Extract the vector directly
original_distance_data  <- benchmark_mapping_results_summary$distance_vectors_df[[num]][[1]]



# Calculate mean and median from the original data
original_mean <- mean(original_distance_data)
original_median <- median(original_distance_data)

# Modify the data: replace values greater than 200 with 200 for the final bin
modified_distance_data <- ifelse(original_distance_data > 25, 25, original_distance_data)

# Convert the modified vector to a data frame for plotting
data <- data.frame(distance = modified_distance_data)

# Plotting
plt <- ggplot(data, aes(x=distance, fill = caller_name)) + 
  stat_bin(aes(y=..count../sum(..count..)), geom="col", bins = 500, breaks = c(seq(0, 24, by = 1), 25)) + 
  scale_fill_manual(values = colors) +  # Apply your custom colors
  geom_vline(aes(xintercept=original_mean, color="Mean"), linetype="dashed") +
  geom_vline(aes(xintercept=original_median, color="Median"), linetype="dashed") +
  scale_color_manual(values=c("red", "blue"), labels=c("Mean", "Median")) + 
  
  # Customize x-axis to include "25+" label
  scale_x_continuous(
    breaks = c(seq(0, 20, by = 10), 24),
    labels = c("0", "10", "20", "25+"),
    limits = c(-2, 26)  # Adjust limits to accommodate "25+"
  ) +

  labs(x = "Distance from true TE insertion (bp)", 
       y = "Proportion of calls", 
       title = caller_name,
       color = NULL  # Remove legend title for Mean and Median
  ) + 
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), 
        plot.title=element_text(size=12), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +
  expand_limits(x=c(-2, 25), y=c(0, 0.6)) +
  annotate(
    "text",
    x = 12.5,
    y = 0.375,  # Adjust the y-coordinate as needed
    label = paste("Mean: ", round(original_mean, 2)),
    size = 5
  ) +
  annotate(
    "text",
    x = 12.5,
    y = 0.45,  # Adjust the y-coordinate as needed
    label = paste("Median: ", round(original_median, 2)),
    size = 5
  ) + theme_minimal()
  return(plt)
}

#temp_bp_plot <- make_bp_plot("temp", benchmark_mapping_results_summary)
#temp2_bp_plot <- make_bp_plot("temp2", benchmark_mapping_results_summary)
#retroseq_bp_plot <- make_bp_plot("retroseq", benchmark_mapping_results_summary)
#teflon_bp_plot <- make_bp_plot("teflon", benchmark_mapping_results_summary)
#popoolationte_bp_plot <- make_bp_plot("popoolationte", benchmark_mapping_results_summary)
#popoolationte2_bp_plot <- make_bp_plot("popoolationte2")
#teforest_bp_plot <- make_bp_plot(paste0(caller_name2, "_bps"), benchmark_mapping_results_summary)

#p2 <- ggarrange(teforest_bp_plot, temp2_bp_plot, teflon_bp_plot, retroseq_bp_plot,
#                    labels = c("A", "B", "C", "D"),
#                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
#                    ncol = 2, nrow = 2,
#                    common.legend = TRUE, legend = "right", vjust=1)
#p2

create_plot_dir <- function(dir) {
if (!dir.exists(paste0(plt_dir, dir))) {
  # Create the directory
  dir.create(paste0(plt_dir, "/", dir, "/"))
}}
plot_dirs <- c("performance_plots", "breakpoints_plots", "hetplt", "homoplt", "freqpltnozero", "freqplt")
mapply(create_plot_dir, plot_dirs)

ggsave(paste0(plt_dir, "/performance_plots/", genome1, "_", genome2, ".pdf"), p, width=14)
#ggsave(paste0(plt_dir, "/breakpoints_plots/",genome1, "_", genome2, ".png"), p2, width=8.5, height=8.5)
ggsave(paste0(plt_dir, "/hetplt/",genome1, "_", genome2, ".pdf"), hetplt, width=14)
ggsave(paste0(plt_dir, "/homoplt/",genome1, "_", genome2, ".pdf"), homoplt, width=14)
ggsave(paste0(plt_dir, "/freqpltnozero/",genome1, "_", genome2, ".pdf"), freqpltnozero, width=14, height=10)
ggsave(paste0(plt_dir, "/freqplt/",genome1, "_", genome2, ".pdf"), freqplt, width=14, height=10)

#ggsave(paste0(plt_dir, "/confusion_matrices/", genome, "mycaller", ".pdf"), conf_matrix_mycaller, width=10, height=8.5)

benchmark_mapping_results_summary_plot$genome <- genome
print(paste0(genome, " done!"))
#return(benchmark_mapping_results_summary_plot)
#}
#genome1s <- c("JUT-008", "A2", "A4", "A6", "B1", "B3")
#genome2s <- c("MUN-009","A3", "A5", "A7", "B2", "B4")
#genome1s <- c("A4", "A6", "B1", "B3")
#genome2s <- c("A5", "A7", "B2", "B4")
#genome1s <- c("JUT-008")
#genome2s <- c("MUN-009")
#summary <- mapply(benchmark, genome1s, genome2s, SIMPLIFY = FALSE)
#summary_combined <- bind_rows(summary)
#takes awhile to run, so saving a copy of the data for easy access to the results
#results_dir <- "/nas/longleaf/home/adaigle/work/mcclintock_stuff/find_candidate_regions/heterozygosity_experiments/outputs/plots/"
#saveRDS(summary_combined, file=paste0(results_dir, "benchmark_summary.rds"))

#loaded_summary <- readRDS(paste0(results_dir, "benchmark_summary.rds")) %>% 
#  filter(!caller %in% c("allseeingeye","allseeingeye2","allseeingeye3", "mystuff", "mystuff_with_bps")) %>%
#  group_by(caller) 
#
#
#df_long <- pivot_longer(loaded_summary, cols = c(precision, recall, f1_score), names_to = "Metric", values_to = "Score")
#
#df_long$caller <- factor(df_long$caller, levels = benchmark_mapping_results_summary_plot$caller)
#
##loaded_summary_summarized <- loaded_summary %>%
##  summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
#
#average_scores <- loaded_summary %>%
#  group_by(caller) %>%
#  summarize(avg_f1_score = mean(f1_score, na.rm = TRUE))
#
#ggplot(loaded_summary, aes(x = reorder(caller, -f1_score), y = f1_score, fill = caller)) +
#  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
#  geom_errorbar(stat = "summary", fun.data = "mean_sd", width = 0.2, position = position_dodge(width = 0.9)) +
#  labs(y = "Average F1 Score", x = "Caller") +
#  geom_text(data = average_scores, aes(x = caller, y = avg_f1_score, label = round(avg_f1_score, 4)), vjust = 5)
#
#ggplot(loaded_summary, aes(x = reorder(caller, -f1_score), y = f1_score, fill = caller)) +
#  #geom_bar(stat = "summary", fun = "mean", position = "dodge") +
#  #geom_errorbar(stat = "summary", fun.data = "mean_sd", width = 0.2, position = position_dodge(width = 0.9)) +
#  geom_boxplot() +
#  labs(y = "Average F1 Score", x = "Caller") +
#  geom_text(data = average_scores, aes(x = caller, y = avg_f1_score, label = round(avg_f1_score, 4)), vjust = -3)
#
#
##boxplot with f1, precision, recall
#ggplot(df_long, aes(x = caller, y = Score, fill = Metric)) +
#  #geom_bar(stat = "summary", fun = "mean", position = "dodge") +
#  #geom_errorbar(stat = "summary", fun.data = "mean_sd", width = 0.2, position = position_dodge(width = 0.9)) +
#    #geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
#  geom_boxplot(position = position_dodge(width = 0.8)) +
#  labs(y = "Average F1 Score", x = "Caller") 
#  #geom_text(data = average_scores, aes(x = caller, y = avg_f1_score, label = round(avg_f1_score, 4)), vjust = -3)
#
#ggplot(df_long, aes(x = caller, y = Score, fill = Metric)) +
#  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
#  geom_errorbar(stat = "summary", fun.data = "mean_sd", width = 0.2, position = position_dodge(width = 0.9)) +
#    #geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
#  #geom_boxplot(position = position_dodge(width = 0.8)) +
#  labs(y = "Average F1 Score", x = "Caller")

#START REFERENCE 


#of my true positives, find which ones are homozygotes called as heterozygotes, etc
benchmark_mapping_results_het_ref <- benchmark_mapping_results_ref %>%
  mutate(
    A1_truth_forTE_het = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x, heterozygosity >= 0.25, heterozygosity <= 0.75),  keep.extra.columns = T )),
    A1_truth_ref_het = map(A1_truth_forTE_het,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref_het = map(A1_truth_forTE_het, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    A1_truth_forTE_homo = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x, , heterozygosity >= 0.75),  keep.extra.columns = T )),
    A1_truth_ref_homo = map(A1_truth_forTE_homo,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref_homo = map(A1_truth_forTE_homo, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    ref_true_positives_homo = map2(ref_true_positives, A1_truth_ref_homo, 
        ~ as.data.frame(subsetByOverlaps(.x, .y, ignore.strand = TRUE, type=c("equal")))),
    homo_tp = map(ref_true_positives_homo, ~ .x %>% filter(heterozygosity >= 0.75)),
    homo_fp_het = map(ref_true_positives_homo, ~ .x %>% filter(heterozygosity >= 0.25, heterozygosity <= 0.75)),
    homo_tp_length = unlist(map(homo_tp, ~ nrow(.x))),
    homo_fp_het_length = unlist(map(homo_fp_het, ~ nrow(.x))),
    homo_fn = map(ref_false_negatives, ~ .x %>% as.data.frame %>% filter(heterozygosity==1.0)),
    homo_fn_length = unlist(map(homo_fn, ~ sum(.x$heterozygosity==1.0))),
    ref_true_positives_het = map2(ref_true_positives, A1_truth_ref_het, 
        ~ as.data.frame(subsetByOverlaps(.x, .y, ignore.strand = TRUE, type=c("equal")))),
    het_tp = map(ref_true_positives_het, ~ .x %>% filter(heterozygosity >= 0.25, heterozygosity <= 0.75)),
    #eventually need to differentiate FPs that are below 0.25, as these are a different class of error
    het_fp_homo = map(ref_true_positives_het, ~ .x %>% filter(heterozygosity > 0.75)),
    het_tp_length = unlist(map(het_tp, ~ nrow(.x))),
    het_fp_homo_length = unlist(map(het_fp_homo, ~ nrow(.x))),
    het_fn = map(ref_false_negatives, ~ .x %>% as.data.frame %>% filter(heterozygosity==0.5)),
    het_fn_length = unlist(map(het_fn, ~ sum(.x$heterozygosity==0.5))),
    homo_fp_zero = map(ref_true_positives_homo, ~ .x %>% filter(heterozygosity < 0.25)),
    het_fp_zero = map(ref_true_positives_het, ~ .x %>% filter(heterozygosity < 0.25)),
    homo_fp_zero_length = unlist(map(homo_fp_zero, ~ nrow(.x))),
    het_fp_zero_length = unlist(map(het_fp_zero, ~ nrow(.x))),
    fp = map(ref_false_positives, ~ .x %>% as.data.frame),
    homo_fp_length = unlist(map(fp, ~ sum(.x$heterozygosity<=0.75))),
    het_fp_length = unlist(map(fp, ~ sum(.x$heterozygosity>0.75))),
  )

#test <- benchmark_mapping_results_het %>% filter(caller=="TEforest_classifier")
#test2 <- benchmark_mapping_results_het %>% filter(caller=="2L_test")

benchmark_mapping_results_frequency_plot_ref <- benchmark_mapping_results_het_ref %>%
  mutate(
    homo_tp_labeled = map2(caller, homo_tp, ~ .y %>% mutate(caller = .x, truth = 1.0)),
    het_tp_labeled = map2(caller, het_tp, ~ .y %>% mutate(caller = .x, truth = 0.5)),
    homo_fp_het_labeled = map2(caller, homo_fp_het, ~ .y %>% mutate(caller = .x, truth = 1.0)),
    het_fp_homo_labeled = map2(caller, het_fp_homo, ~ .y %>% mutate(caller = .x, truth = 0.5)),
    homo_fn_labeled = map2(caller, homo_fn, ~ .y %>% mutate(caller = .x, truth = 1.0, ID = NULL, heterozygosity=0, length = NULL)),
    het_fn_labeled = map2(caller, het_fn, ~ .y %>% mutate(caller = .x, truth = 0.5, heterozygosity=0, ID = NULL, length = NULL)),
    homo_fp_zero_labeled = map2(caller, homo_fp_zero, ~ .y %>% mutate(caller = .x, truth = 1.0)),
    het_fp_zero_labeled = map2(caller, het_fp_zero, ~ .y %>% mutate(caller = .x, truth = 0.5)),
    fp = map2(caller, fp, ~ .y %>% mutate(caller = .x, truth = 0)),
  )

# Combine the labeled columns into a single data frame
frequency_plot_ref <- do.call(rbind, c(
  benchmark_mapping_results_frequency_plot_ref$homo_tp_labeled,
  benchmark_mapping_results_frequency_plot_ref$het_tp_labeled,
  benchmark_mapping_results_frequency_plot_ref$homo_fp_het_labeled,
  benchmark_mapping_results_frequency_plot_ref$het_fp_homo_labeled,
  benchmark_mapping_results_frequency_plot_ref$homo_fn_labeled,
  benchmark_mapping_results_frequency_plot_ref$het_fn_labeled,
  benchmark_mapping_results_frequency_plot_ref$homo_fp_zero_labeled,
  benchmark_mapping_results_frequency_plot_ref$het_fp_zero_labeled,
  benchmark_mapping_results_frequency_plot_ref$fp
))


# Calculate R-squared value
calculate_r_squared <- function(data) {
  lm_model <- lm(truth ~ as.numeric(heterozygosity), data = data)
  return(data.frame(R_squared = summary(lm_model)$r.squared))
}
r_squared_results_ref <- frequency_plot_ref %>%
  group_by(caller) %>%
  do(calculate_r_squared(.))

frequency_plot_ref <- left_join(frequency_plot_ref, r_squared_results_ref, by = "caller") %>% mutate(heterozygosity = as.numeric(heterozygosity))

freqplt_ref <- ggplot(frequency_plot_ref, aes(x = truth, y = heterozygosity, color = caller)) +
  geom_jitter(alpha = 0.2, width = 0.1) +
  geom_boxplot(aes(group=truth), width = 0.05, fill = "lightgray", outlier.shape = NA) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~caller) +
  labs(title = "Scatterplot with Boxplot of True vs Predicted Values",
       x = "True Values",
       y = "Predicted Values") 

frequency_plot_nozero_ref <- frequency_plot_ref %>% filter(heterozygosity!=0)
freqpltnozero_ref <-ggplot(frequency_plot_nozero_ref, aes(x = truth, y = heterozygosity, color = caller)) +
  geom_jitter(alpha = 0.2, width = 0.1) +
  geom_boxplot(aes(group=truth), width = 0.05, fill = "lightgray", outlier.shape = NA) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~caller) +
  labs(title = "Scatterplot with Boxplot of True vs Predicted Values",
       x = "True Values",
       y = "Predicted Values") 

benchmark_mapping_results_ref_nodata <- benchmark_mapping_results_het_ref[c(1,7,8, 22:29, 39:40, 42:43, 46:47, 49, 52:53, 55:56)]

# Check for NaN values and replace with 0
benchmark_mapping_results_ref_nodata[is.na(benchmark_mapping_results_ref_nodata)] <- 0

#length of reduced nonref DF
#will mess up reduced nested TEs of different classes
#but necessary for proper analysis 
ref_te_number_het <- length(unique(truth_het_reference_gr))

ref_te_number_homo <- length(unique(truth_homo_reference_gr))

benchmark_mapping_results_summary_ref <- benchmark_mapping_results_ref_nodata %>%
  group_by(caller) %>%
  summarize(
    #per_fam_mean_f1_score = mean(f1_score, na.rm = TRUE),
    #per_fam_mean_precision = mean(precision, na.rm = TRUE),
    #per_fam_mean_recall = mean(recall, na.rm = TRUE),
    sum_ref_true_positives_length = sum(ref_true_positives_length),
    sum_true_positives = sum(ref_calls_true_positives_length),
    sum_false_positives = sum(ref_calls_false_positives_length),
    sum_false_negatives = sum(ref_calls_false_negatives_length),
    recalc_false_negatives = num_true_ref - sum_true_positives,
    precision = sum_true_positives / (sum_true_positives + sum_false_positives),
    recall = sum_true_positives / (sum_true_positives + recalc_false_negatives),
    f1_score = 2 * (precision * recall) / (precision + recall), 
    homo_tp_length = sum(homo_tp_length),
    het_tp_length = sum(het_tp_length),
    homo_fp_het_length = sum(homo_fp_het_length),
    het_fp_homo_length = sum(het_fp_homo_length),
    #homo_fn_length = sum(homo_fn_length),
    homo_fn_length = ref_te_number_homo - homo_tp_length,
    #het_fn_length = sum(het_fn_length),
    het_fn_length = ref_te_number_het - het_tp_length,
    homo_fp_zero_length = sum(homo_fp_zero_length),
    het_fp_zero_length = sum(het_fp_zero_length),
    homo_fp_length = sum(homo_fp_length),
    het_fp_length = sum(het_fp_length),
  )

benchmark_mapping_results_summary_ref[is.na(benchmark_mapping_results_summary_ref)] <- 0

het_conf_matrix_ref <- function(caller_name) {
num <- as.numeric(row.names(benchmark_mapping_results_summary_ref)[benchmark_mapping_results_summary_ref$caller == caller_name])

conf_matrix_mycaller <- benchmark_mapping_results_summary_ref[num,]
#heterozygosity classification accuracy
conf_matrix_values <- matrix(c(0, conf_matrix_mycaller$het_fn_length + conf_matrix_mycaller$het_fp_zero_length, conf_matrix_mycaller$homo_fn_length + conf_matrix_mycaller$homo_fp_zero_length, # need to grab fp hets and homos
                                conf_matrix_mycaller$het_fp_length, conf_matrix_mycaller$het_tp_length, conf_matrix_mycaller$homo_fp_het_length, 
                                conf_matrix_mycaller$homo_fp_length, conf_matrix_mycaller$het_fp_homo_length, conf_matrix_mycaller$homo_tp_length),
                             nrow = 3, byrow = TRUE)
conf_matrix <- confusionMatrix(conf_matrix_values)
cm <- conf_mat(conf_matrix_values)

plt <- autoplot(cm, type = "heatmap") +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1") + 
  theme(legend.position = "right")+
  scale_x_discrete(labels = c("absent", "heterozygote", "homozygote")) +
  scale_y_discrete(labels = c("homozygote", "heterozygote", "absent"))

#summary(cm)[1,]$.estimate
return(plt)
}
#conf_matrix_mycaller <- het_conf_matrix("2L_test")
#et_conf_matrix("temp")
#et_conf_matrix("temp2")
#et_conf_matrix("teflon")
#et_conf_matrix("retroseq")
#et_conf_matrix("popoolationte")
#et_conf_matrix("popoolationte2")



het_conf_matrix_ref <- benchmark_mapping_results_summary_ref %>%
  rowwise() %>%
  mutate(
    conf_matrix = list(matrix(c(0, het_fn_length + het_fp_zero_length, homo_fn_length + homo_fp_zero_length,
                                 het_fp_length, het_tp_length, homo_fp_het_length,
                                 homo_fp_length, het_fp_homo_length, homo_tp_length),
                              nrow = 3, byrow = TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    conf_matrix_summary = map(conf_matrix, ~confusionMatrix(.x)),
    het_stats = map(conf_matrix_summary, ~ c(precision=.x$byClass["Class: B","Precision"], 
      recall=.x$byClass["Class: B","Recall"], f1=.x$byClass["Class: B","F1"])),
    homo_stats = map(conf_matrix_summary, ~ c(precision=.x$byClass["Class: C","Precision"], 
      recall=.x$byClass["Class: C","Recall"], f1=.x$byClass["Class: C","F1"])),
  )

# Expand the list columns into separate columns
het_stats_plot_ref <- het_conf_matrix_ref %>% select(caller, het_stats) %>%
  tidyr::unnest_wider(het_stats) %>%
  gather(metric, value, precision:f1) %>%
  group_by(caller) 
# Calculate the average F1 score for each caller and order by this score
f1_order <- het_stats_plot_ref %>%
  filter(metric == "f1") %>%
  group_by(caller) %>%
  summarise(mean_f1 = mean(value)) %>%
  arrange(desc(mean_f1)) %>%
  pull(caller)

# Reorder the caller factor levels based on F1 score
het_stats_plot_ref <- het_stats_plot_ref %>%
  mutate(caller = factor(caller, levels = f1_order))

# Create the bar plot
hetplt_ref <- ggplot(het_stats_plot_ref, aes(x = caller, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Caller", y = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),  # Adjust legend title size
        legend.text = element_text(size = 12),   # Adjust legend text size
        axis.text = element_text(size = 12),     # Adjust axis text size
        axis.title = element_text(size = 14))  +  # Adjust axis title size
  expand_limits(y=c(0,1)) +
# Add F1 score numbers above the bars
 geom_text(aes(label = sprintf("%.2f", value), y = value, group = metric),
              position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)


# Expand the list columns into separate columns
homo_stats_plot_ref <- het_conf_matrix_ref %>% select(caller, homo_stats) %>%
  tidyr::unnest_wider(homo_stats) %>%
  gather(metric, value, precision:f1) %>%
  group_by(caller)

f1_order_homo <- homo_stats_plot_ref %>%
  filter(metric == "f1") %>%
  group_by(caller) %>%
  summarise(mean_f1 = mean(value)) %>%
  arrange(desc(mean_f1)) %>%
  pull(caller)

# Reorder the caller factor levels based on F1 score
homo_stats_plot_ref <- homo_stats_plot_ref %>%
  mutate(caller = factor(caller, levels = f1_order_homo))


# Create the bar plot
homoplt_ref <- ggplot(homo_stats_plot_ref, aes(x = caller, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Caller", y = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),  # Adjust legend title size
        legend.text = element_text(size = 12),   # Adjust legend text size
        axis.text = element_text(size = 12),     # Adjust axis text size
        axis.title = element_text(size = 14))  +  # Adjust axis title size
  expand_limits(y=c(0,1)) +
# Add F1 score numbers above the bars
 geom_text(aes(label = sprintf("%.2f", value), y = value, group = metric),
              position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

#benchmark_mapping_results_summary_plot_ref <- benchmark_mapping_results_summary_ref #%>% #filter(caller !="ngs_te_mapper")%>% filter(caller !="ngs_te_mapper2") %>%
    #filter(caller !=c("allseeingeye", "allseeingeye2", "allseeingeye3"))
#benchmark_mapping_results_summary_plot <- benchmark_mapping_results_summary_plot #%>% 
  #mutate(caller = case_when(
  #  caller == "allseeingeye4" ~ "mystuff",
  #  caller == "allseeingeye5" ~ "mystuff_with_bps",
  #  TRUE ~ caller  # Keep other values unchanged
  #))
benchmark_mapping_results_summary_plot_ref <- benchmark_mapping_results_summary_ref %>%
  arrange(desc(f1_score))



# Reshape the data to long format
df_long_ref <- pivot_longer(benchmark_mapping_results_summary_plot_ref, cols = c(precision, recall, f1_score), names_to = "Metric", values_to = "Score")

df_long_ref$caller <- factor(df_long_ref$caller, levels = benchmark_mapping_results_summary_plot_ref$caller)
df_long_ref_save <- df_long_ref
df_long_ref_save$covg <- paste0(coverage, "X")
# Create a grouped bar plot
p_ref <- ggplot(df_long_ref, aes(x = caller, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  #scale_fill_manual(values = c("precision" = "blue", "recall" = "green", "f1_score" = "red")) +
  labs(y = "Score",
       x = "Caller") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),  # Adjust legend title size
        legend.text = element_text(size = 12),   # Adjust legend text size
        axis.text = element_text(size = 12),     # Adjust axis text size
        axis.title = element_text(size = 14))  +  # Adjust axis title size
  expand_limits(y=c(0,1)) +
# Add F1 score numbers above the bars
 geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = Metric),
              position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)#return(benchmark_mapping_results_summary)


ggsave(paste0(plt_dir, "/performance_plots/", genome1, "_", genome2, "_ref.pdf"), p_ref, width=14)
ggsave(paste0(plt_dir, "/hetplt/",genome1, "_", genome2, "_ref.pdf"), hetplt_ref, width=14)
ggsave(paste0(plt_dir, "/homoplt/",genome1, "_", genome2, "_ref.pdf"), homoplt_ref, width=14)
ggsave(paste0(plt_dir, "/freqpltnozero/",genome1, "_", genome2, "_ref.pdf"), freqpltnozero_ref, width=14, height=10)
ggsave(paste0(plt_dir, "/freqplt/",genome1, "_", genome2, "_ref.pdf"), freqplt_ref, width=14, height=10)

saveRDS(df_long_save, file = paste0(plt_dir, "/performance_plots/", genome1, "_", genome2, ".rds"))
saveRDS(df_long_ref_save, file = paste0(plt_dir, "/performance_plots/", genome1, "_", genome2, "_ref", ".rds"))


## Finding a common set of nonref tp calls

# Function to find common GRanges across multiple callers for each TE
find_common_granges <- function(gr_list) {
  #if (length(gr_list) < 2) return(gr_list[[1]])  # If there's only one or no GRanges, return as is
  
  # Start with the first GRanges object
  common_gr <- gr_list[[1]]
  
  # Intersect with the remaining GRanges objects
  for (i in 2:length(gr_list)) {
    common_gr <- subsetByOverlaps(common_gr, gr_list[[i]], ignore.strand = TRUE, type='equal')
  }
  
  return(common_gr)
}

common_calls_df <- benchmark_mapping_results_frequency_plot %>% mutate(  
    nonref_true_positives_truthrange = map2(A1_truth_nonref, nonref_true_positives, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE))) 


common_calls_benchmarks <- function(callers) {


# choose callers
test <- common_calls_df %>% filter(caller %in% callers)
#test <- common_calls_df %>% filter(caller %in% c('truth', 'TEforest_classifier_filter_bps', 'retroseq', 'temp2', 'teflon'))
#test <- common_calls_df %>% filter(caller %in% c('temp2', 'truth', 'teflon', 'TEforest_classifier_filter_bps'))

total_callers <- test %>% 
  distinct(caller) %>% 
  nrow()
# Group by TE and find common GRanges across callers
common_granges_per_TE <- test %>%
  group_by(TE) %>%
  filter(n() == total_callers) %>%  # Filter to ensure the TE is present for all callers
  summarise(common_granges = list(find_common_granges(nonref_true_positives_truthrange)), .groups = 'drop') %>%
  mutate(common_granges_count = map_int(common_granges, length))

# Print the results
#print(common_granges_per_TE)
#print(sum(common_granges_per_TE$common_granges_count))

benchmark_with_common <- test %>%
  left_join(common_granges_per_TE, by = "TE") %>%
  filter(!is.na(common_granges_count) & common_granges_count != 0) %>% 
  mutate(  
    common_granges_callercoordinates = map2(nonref_true_positives, common_granges, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE))) 

distance_calcs_common_calls <- benchmark_with_common %>%
  mutate(
    distances_vector = map2(common_granges_callercoordinates, A1_truth_nonref_noextend, ~ distanceToNearest(center(.x), .y, ignore.strand = TRUE)),
    distances_vector2 = map(distances_vector, ~ .x@elementMetadata[[1]]),
    distances_vector_mean = unlist(map(distances_vector, ~ mean(.x@elementMetadata[[1]]))),
    distances_vector_sd = unlist(map(distances_vector, ~ sd(.x@elementMetadata[[1]])))
  )

#of my true positives, find which ones are homozygotes called as heterozygotes, etc
benchmark_mapping_results_het_common_calls <- distance_calcs_common_calls %>%
  mutate(
    A1_truth_forTE_het = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x, heterozygosity >= 0.25, heterozygosity <= 0.75),  keep.extra.columns = T )),
    A1_truth_ref_het = map(A1_truth_forTE_het,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref_het = map(A1_truth_forTE_het, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    A1_truth_forTE_homo = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(A1_truth %>% filter(TE==.x, , heterozygosity >= 0.75),  keep.extra.columns = T )),
    A1_truth_ref_homo = map(A1_truth_forTE_homo,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    A1_truth_nonref_homo = map(A1_truth_forTE_homo, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    nonref_true_positives_homo = map2(common_granges_callercoordinates, A1_truth_nonref_homo, 
        ~ as.data.frame(subsetByOverlaps(.x, .y, ignore.strand = TRUE))),
    homo_tp = map(nonref_true_positives_homo, ~ .x %>% filter(heterozygosity >= 0.75)),
    homo_fp_het = map(nonref_true_positives_homo, ~ .x %>% filter(heterozygosity >= 0.25, heterozygosity <= 0.75)),
    homo_tp_length = unlist(map(homo_tp, ~ nrow(.x))),
    homo_fp_het_length = unlist(map(homo_fp_het, ~ nrow(.x))),
    homo_fn = map(nonref_false_negatives, ~ .x %>% as.data.frame %>% filter(heterozygosity==1.0)),
    homo_fn_length = unlist(map(homo_fn, ~ sum(.x$heterozygosity==1.0) )),
    nonref_true_positives_het = map2(common_granges_callercoordinates, A1_truth_nonref_het, 
        ~ as.data.frame(subsetByOverlaps(.x, .y, ignore.strand = TRUE))),
    het_tp = map(nonref_true_positives_het, ~ .x %>% filter(heterozygosity >= 0.25, heterozygosity <= 0.75)),
    #eventually need to differentiate FPs that are below 0.25, as these are a different class of error
    het_fp_homo = map(nonref_true_positives_het, ~ .x %>% filter(heterozygosity > 0.75)),
    het_tp_length = unlist(map(het_tp, ~ nrow(.x))),
    het_fp_homo_length = unlist(map(het_fp_homo, ~ nrow(.x))),
    het_fn = map(nonref_false_negatives, ~ .x %>% as.data.frame %>% filter(heterozygosity==0.5)),
    het_fn_length = unlist(map(het_fn, ~ sum(.x$heterozygosity==0.5))),
    homo_fp_zero = map(nonref_true_positives_homo, ~ .x %>% filter(heterozygosity < 0.25)),
    het_fp_zero = map(nonref_true_positives_het, ~ .x %>% filter(heterozygosity < 0.25)),
    homo_fp_zero_length = unlist(map(homo_fp_zero, ~ nrow(.x))),
    het_fp_zero_length = unlist(map(het_fp_zero, ~ nrow(.x))),
    fp = map(nonref_false_positives, ~ .x %>% as.data.frame),
    homo_fp_length = unlist(map(fp, ~ sum(.x$heterozygosity<=0.75))),
    het_fp_length = unlist(map(fp, ~ sum(.x$heterozygosity>0.75))),
  )

benchmark_mapping_results_het_common_calls <- benchmark_mapping_results_het_common_calls %>% mutate(
  true_negatives_truth = map(
      TE, 
      ~ {
        df <- training_csv_negative_downsampled %>% filter(TE == .x) 
        if (nrow(df) == 0) {
          return(GRanges())  # Return empty grange if no true negatives for this TE fam
        } else {
          return(makeGRangesFromDataFrame(df, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE))
        }
      }
    ),
  true_negatives = map2(true_negatives_truth, nonreference_freq_filter,
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, invert=TRUE)),
  true_negatives_length = unlist(map(true_negatives, ~ length(.x))),
)

benchmark_mapping_results_common_calls_nodata <- benchmark_mapping_results_het_common_calls[c(1,7,8,20,24:30,33, 45:48, 52:53, 55, 58:59, 61:62, 65)]

# Check for NaN values and replace with 0
benchmark_mapping_results_common_calls_nodata[is.na(benchmark_mapping_results_common_calls_nodata)] <- 0

debug_breakpoints <- benchmark_mapping_results_het_common_calls %>% filter(caller=="TEforest_classifier_filter_bps") %>% filter(TE=="roo")
debug_breakpointstmp2 <- benchmark_mapping_results_het_common_calls %>% filter(caller=="temp2") %>% filter(TE=="roo")

#length of reduced nonref DF
#will mess up reduced nested TEs of different classes
#but necessary for proper analysis 
nonref_te_number <- length(unique(truth_nonreference_extend_gr)) #SWAP
truth_het_nonreference_gr <- subsetByOverlaps(truth_het_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert=TRUE)
nonref_te_number_het <- length(unique(truth_het_nonreference_gr))

truth_homo_nonreference_gr <- subsetByOverlaps(truth_homo_gr, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert=TRUE)
nonref_te_number_homo <- length(unique(truth_homo_nonreference_gr))

benchmark_mapping_results_common_calls_summary <- benchmark_mapping_results_common_calls_nodata %>%
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
    distance_mean = mean(unlist(distances_vector2)),
    distance_sd = sd(unlist(distances_vector2)),
    distance_vectors = list(c(distances_vector2)),
    homo_tp_length = sum(homo_tp_length),
    het_tp_length = sum(het_tp_length),
    homo_fp_het_length = sum(homo_fp_het_length),
    het_fp_homo_length = sum(het_fp_homo_length),
    #homo_fn_length = sum(homo_fn_length),
    homo_fn_length = nonref_te_number_homo - homo_tp_length,
    #het_fn_length = sum(het_fn_length),
    het_fn_length = nonref_te_number_het - het_tp_length,
    homo_fp_zero_length = sum(homo_fp_zero_length),
    het_fp_zero_length = sum(het_fp_zero_length),
    homo_fp_length = sum(homo_fp_length),
    het_fp_length = sum(het_fp_length),
    het_fp_length = sum(het_fp_length),
    true_negatives_length = sum(true_negatives_length),
  )

benchmark_mapping_results_common_calls_summary$recalc_false_negatives <- ifelse(benchmark_mapping_results_common_calls_summary$recalc_false_negatives < 0, 0, benchmark_mapping_results_common_calls_summary$recalc_false_negatives)
benchmark_mapping_results_common_calls_summary$homo_fn_length <- ifelse(benchmark_mapping_results_common_calls_summary$homo_fn_length < 0, 0, benchmark_mapping_results_common_calls_summary$homo_fn_length)
benchmark_mapping_results_common_calls_summary$het_fn_length <- ifelse(benchmark_mapping_results_common_calls_summary$het_fn_length < 0, 0, benchmark_mapping_results_common_calls_summary$het_fn_length)

benchmark_mapping_results_common_calls_summary <- benchmark_mapping_results_common_calls_summary %>% mutate(
    precision = sum_true_positives / (sum_true_positives + sum_false_positives),
    recall = sum_true_positives / (sum_true_positives + recalc_false_negatives),
    f1_score = 2 * (precision * recall) / (precision + recall)
)

benchmark_mapping_results_common_calls_summary[is.na(benchmark_mapping_results_common_calls_summary)] <- 0

benchmark_mapping_results_common_calls_summary <- benchmark_mapping_results_common_calls_summary %>%
  mutate(distance_vectors = map(distance_vectors, ~ unlist(unname(.x))),
  distance_vectors_df = map(distance_vectors, ~ data.frame(distance = .x)))

benchmark_mapping_results_common_calls_summary <- benchmark_mapping_results_common_calls_summary %>%
  mutate(
    # F1 Score for Genotyping (Class 1)
    genotyping_recall = homo_tp_length / (homo_tp_length + homo_fp_het_length + homo_fp_zero_length),
    genotyping_precision = homo_tp_length / (homo_tp_length + het_fp_homo_length + het_fp_zero_length), 
    genotyping_f1_score = if_else(
      (genotyping_recall + genotyping_precision) == 0,
      0,
      (2 * genotyping_recall * genotyping_precision) / (genotyping_recall + genotyping_precision)
    ),
    
    # F1 Score for Het (Class 2)
    het_recall = het_tp_length / (het_tp_length + het_fp_homo_length + het_fp_zero_length),
    het_precision = het_tp_length / (het_tp_length + homo_fp_het_length + het_fp_zero_length),
    het_f1_score = if_else(
      (het_recall + het_precision) == 0,
      0,
      (2 * het_recall * het_precision) / (het_recall + het_precision)
    ),
    
    # Macro-Averaged F1 Score
    macro_f1_score = (genotyping_f1_score + het_f1_score) / 2
  )

#retroseq_commoncalls <- make_bp_plot("retroseq", benchmark_mapping_results_common_calls_summary)
#TEforest_commoncalls <- make_bp_plot("TEforest_regressor_bps", benchmark_mapping_results_common_calls_summary)
#temp2_commoncalls <- make_bp_plot("temp2", benchmark_mapping_results_common_calls_summary)
#teflon_commoncalls <- make_bp_plot("teflon", benchmark_mapping_results_common_calls_summary)
#popoo_commoncalls <- make_bp_plot("popoolationte", benchmark_mapping_results_common_calls_summary)
#temp_commoncalls <- make_bp_plot("temp", benchmark_mapping_results_common_calls_summary)

#ggarrange(TEforest_commoncalls, temp_commoncalls, temp2_commoncalls, teflon_commoncalls, retroseq_commoncalls, popoo_commoncalls,
#                    labels = c("A", "B", "C", "D", "E", "F"),
#                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
#                    ncol = 2, nrow = 3,
#                    common.legend = TRUE, legend = "right", vjust=1)

#ggarrange(TEforest_commoncalls, temp2_commoncalls, retroseq_commoncalls,
#                    labels = c("A", "B"),
#                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
#                    ncol = 3, nrow = 1,
#                    common.legend = TRUE, legend = "right", vjust=1)
benchmark_mapping_results_common_calls_summary_clean <- benchmark_mapping_results_common_calls_summary
benchmark_mapping_results_common_calls_summary_clean$distance_vectors <- NULL
benchmark_mapping_results_common_calls_summary_clean$distance_vectors_df <- NULL

return(benchmark_mapping_results_common_calls_summary_clean)
}


six_callers <- c('temp2', 'truth', 'TEforest_classifier_filter_bps', 'retroseq', 'teflon', 'temp', 'popoolationte')
three_callers <- c('temp2', 'truth', 'retroseq', 'TEforest_classifier_filter_bps')
two_callers <- c('temp2', 'truth', 'TEforest_classifier_filter_bps')

three_callers_df <- common_calls_benchmarks(three_callers)
#two_callers_df <- common_calls_benchmarks(two_callers)

#write.csv(six_callers_df, file = paste0(plt_dir, "/", genome1, "_", genome2, "_six_common_callers", ".txt"), quote=F, row.names=F)
write.csv(three_callers_df, file = paste0(plt_dir, "/", genome1, "_", genome2, "_three_common_callers", ".txt"), quote=F, row.names=F)

#if (coverage > 10) {
#
#six_callers_df <- common_calls_benchmarks(six_callers)
#write.csv(six_callers_df, file = paste0(plt_dir, "/", genome1, "_", genome2, "_six_common_callers", ".txt"), quote=F, row.names=F)
#
#}
###start of conference specific plots
#load("/nas/longleaf/home/adaigle/work/test_TEforest/test_basenorm_feats/2L_2R_plotsJUT-008_MUN-009.RData")
#load("/nas/longleaf/home/adaigle/work/test_TEforest/test_basenorm_feats/2L_2R_plotsAKA-017_GIM-024.RData")
#
#benchmark_mapping_results_mycaller <- benchmark_mapping_results %>% filter(caller=="TEforest_classifier_filter_bps")
#false_negatives <- do.call(c, unname(benchmark_mapping_results_mycaller$nonref_false_negatives))
#
#colors <- c(
#  "TEforest" = "#F8766D", 
#  "temp2" = "#CD9600", 
#  "retroseq" = "#7CAE00", 
#  "temp" = "#00BE67", 
#  "popoolationte" = "#00BFC4", 
#  "teflon" = "#00A9FF", 
#  "popoolationte2" = "#C77CFF", 
#  "tepid" = "#FF61CC"
#)
#
## Filter and recode the data
#df_long_filtered <- df_long %>% 
#  filter(caller %in% c("TEforest_classifier_filter_bps", "temp", "temp2", "teflon", "retroseq", "popoolationte", "popoolationte2", "tepid")) %>%
#  mutate(caller = recode(caller, 
#                         "TEforest_classifier_filter_bps" = "TEforest",
#                         "temp2" = "temp2",
#                         "retroseq" = "retroseq",
#                         "temp" = "temp",
#                         "popoolationte" = "popoolationte",
#                         "teflon" = "teflon",
#                         "popoolationte2" = "popoolationte2",
#                         "tepid" = "tepid"),
#         Metric = recode(Metric, "precision" = "Precision", "recall" = "Recall", "f1_score" = "F1 Score"))
#
## Create the plot
#nonref_performance <- ggplot(df_long_filtered, aes(x = caller, y = Score, fill = caller)) +
#  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
#  scale_fill_manual(values = colors) +
#  labs(y = "Score",
#       x = "Caller") +
#  theme_minimal() +
#  theme(axis.text.x = element_blank(),
#        legend.title = element_blank(), 
#        legend.text = element_text(size = 12),
#        legend.position = "bottom",
#        axis.text = element_text(size = 12), 
#        axis.title = element_text(size = 14), 
#        strip.text = element_text(size = 16)) + 
#  expand_limits(y = c(0, 1)) +
#  guides(fill = guide_legend(nrow = 1)) +
#  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = Metric),
#            position = position_dodge(width = 0.8), size = 3.75, vjust = -0.5) +
#  facet_wrap(~ Metric)
#
## Save the plot as a PNG file
#output_path <- "/nas/longleaf/home/adaigle/TEforest/plots/nonref_performance_30X.png"
#ggsave(output_path, plot = nonref_performance, width = 10, height = 5)
#
#ggsave("/nas/longleaf/home/adaigle/TEforest/plots/nonref_performance_125_30.png", plot = nonref_performance, width = 12, height = 6, dpi = 300)
#
#df_long_ref_filtered <- df_long_ref %>% 
#  filter(caller %in% c("TEforest_classifier_bps", "temp", "temp2", "teflon", "popoolationte2", "tepid")) %>%
#  mutate(caller = recode(caller, "TEforest_classifier_bps" = "TEforest"), 
#         Metric = recode(Metric, "precision" = "Precision", "recall" = "Recall", "f1_score" = "F1 Score"))
#
#ref_performance <- ggplot(df_long_ref_filtered, aes(x = caller, y = Score, fill = caller)) +
#  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
#  #scale_fill_manual(values = c("precision" = "blue", "recall" = "green", "f1_score" = "red")) +
#  labs(y = "Score",
#       x = "Caller") +
#  theme_minimal() +
#  #scale_fill_viridis_d(option="inferno") +
#  theme(axis.text.x = element_blank(),
#        legend.title = element_blank(), 
#        legend.text = element_text(size = 12),
#        legend.position = "bottom",
#        axis.text = element_text(size = 12), 
#        axis.title = element_text(size = 14), 
#        strip.text = element_text(size=16)) + 
#  expand_limits(y=c(0,1)) +
## Add F1 score numbers above the bars
#  guides(fill = guide_legend(nrow = 1)) +
#  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = Metric),
#              position = position_dodge(width = 0.8), size = 3.75, vjust = -0.5) +
#  facet_wrap(~ Metric)
#ggsave("/nas/longleaf/home/adaigle/TEforest/plots/ref_performance_150_30.png", plot = ref_performance, width = 12, height = 6, dpi = 300)
#


benchmark_mapping_results_genotypes <- benchmark_mapping_results_summary %>%
  mutate(
    # F1 Score for Genotyping (Class 1)
    genotyping_recall = homo_tp_length / (homo_tp_length + homo_fp_het_length + homo_fp_zero_length),
    genotyping_precision = homo_tp_length / (homo_tp_length + het_fp_homo_length + het_fp_zero_length), 
    genotyping_f1_score = if_else(
      (genotyping_recall + genotyping_precision) == 0,
      0,
      (2 * genotyping_recall * genotyping_precision) / (genotyping_recall + genotyping_precision)
    ),
    
    # F1 Score for Het (Class 2)
    het_recall = het_tp_length / (het_tp_length + het_fp_homo_length + het_fp_zero_length),
    het_precision = het_tp_length / (het_tp_length + homo_fp_het_length + het_fp_zero_length),
    het_f1_score = if_else(
      (het_recall + het_precision) == 0,
      0,
      (2 * het_recall * het_precision) / (het_recall + het_precision)
    ),
    
    # Macro-Averaged F1 Score
    macro_f1_score = (genotyping_f1_score + het_f1_score) / 2
  )


unique_callers <- unique(benchmark_mapping_results_het$caller)

# Initialize an empty list to store numeric vectors
fp_heterozygosity_list <- list()

# Loop through each caller
for (caller in unique_callers) {
  
  # Subset the dataframe for the current caller
  subset_df <- benchmark_mapping_results_het[benchmark_mapping_results_het$caller == caller, ]
  
  # Extract the 'value' column (list of GRanges) for the subset
  granges_list <- subset_df$nonref_false_negatives
  
  # Extract 'heterozygosity' from each GRanges object
  # Assuming the metadata column is named 'heterozygosity'
  het_values <- lapply(granges_list, function(gr) mcols(gr)$heterozygosity)
  
  # Flatten the list of heterozygosity values into a single numeric vector
  het_numeric <- unlist(het_values)
  
  # Ensure the vector is numeric
  het_numeric <- as.numeric(het_numeric)
  
  # Store in the list with caller as the name
  fp_heterozygosity_list[[caller]] <- het_numeric
}

# Initialize an empty list to store numeric vectors
tphomo_heterozygosity_list <- list()

# Loop through each caller
for (caller in unique_callers) {
  
  # Subset the dataframe for the current caller
  subset_df <- benchmark_mapping_results_het[benchmark_mapping_results_het$caller == caller, ]
  
  # Extract the 'value' column (list of GRanges) for the subset
  granges_list <- subset_df$nonref_true_positives_homo
  
  # Extract 'heterozygosity' from each GRanges object
  # Assuming the metadata column is named 'heterozygosity'
  het_values <- lapply(granges_list, function(gr) gr$heterozygosity)
  
  # Flatten the list of heterozygosity values into a single numeric vector
  het_numeric <- unlist(het_values)
  
  # Ensure the vector is numeric
  het_numeric <- as.numeric(het_numeric)
  
  # Store in the list with caller as the name
  tphomo_heterozygosity_list[[caller]] <- het_numeric
}

# Initialize an empty list to store numeric vectors
tphet_heterozygosity_list <- list()

# Loop through each caller
for (caller in unique_callers) {
  
  # Subset the dataframe for the current caller
  subset_df <- benchmark_mapping_results_het[benchmark_mapping_results_het$caller == caller, ]
  
  # Extract the 'value' column (list of GRanges) for the subset
  granges_list <- subset_df$nonref_true_positives_het
  
  # Extract 'heterozygosity' from each GRanges object
  # Assuming the metadata column is named 'heterozygosity'
  het_values <- lapply(granges_list, function(gr) gr$heterozygosity)
  
  # Flatten the list of heterozygosity values into a single numeric vector
  het_numeric <- unlist(het_values)
  
  # Ensure the vector is numeric
  het_numeric <- as.numeric(het_numeric)
  
  # Store in the list with caller as the name
  tphet_heterozygosity_list[[caller]] <- het_numeric
}

genotype_freq_df <- tibble(
  caller = unique_callers,
  tphomo = map(unique_callers, ~ if (.x %in% names(tphomo_heterozygosity_list)) tphomo_heterozygosity_list[[.x]] else NA),
  tphet = map(unique_callers, ~ if (.x %in% names(tphet_heterozygosity_list)) tphet_heterozygosity_list[[.x]] else NA),
  fp = map(unique_callers, ~ if (.x %in% names(fp_heterozygosity_list)) fp_heterozygosity_list[[.x]] else NA)
)

benchmark_mapping_results_genotypes_final <- left_join(benchmark_mapping_results_genotypes, genotype_freq_df)

saveRDS(benchmark_mapping_results_genotypes_final, file = paste0(plt_dir, "/", genome1, "_", genome2, "_genotypes.rds"))
saveRDS(confusion_matrix_TEforest, file = paste0(plt_dir, "/", genome1, "_", genome2, "_TEforest_confusion_matrix.rds"))
saveRDS(confusion_matrix_retroseq, file = paste0(plt_dir, "/", genome1, "_", genome2, "_retroseq_confusion_matrix.rds"))
saveRDS(confusion_matrix_temp2, file = paste0(plt_dir, "/", genome1, "_", genome2, "_temp2_confusion_matrix.rds"))
saveRDS(frequency_plot, file = paste0(plt_dir, "/", genome1, "_", genome2, "_freqplot_data.rds"))


save.image(file = paste0(plt_dir, "/", genome1, "_", genome2, ".RData"))

