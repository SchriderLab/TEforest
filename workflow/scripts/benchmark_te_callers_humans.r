#!/usr/bin/env Rscript
# the purpose of this script is to benchmark all TE callers and my caller 
# in an identical fashion. 
# I will benchmark each TE family individually to achive maximum accuracy.
# And calculate the avg distance of the calls from breakpoints. 
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
#library(viridis)
#library(yardstick)

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
genome <- "NA12878"
euchromatin_coordinates_path <- "/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_inputs/chroms.txt"
caller_name <- "TEforest"
plt_dir <- "/nas/longleaf/home/adaigle/work/human_TEs/benchmarks"
basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_5X"
coverage <- "30"
result_path <- paste0("/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_runs/highcov_repeatmasker/NA12878_highcov_1/results/")

#training_csv_human <- read.csv("/nas/longleaf/home/adaigle/work/human_TEs/human_truth_data_NA12878.csv") %>% filter(DET!="DEL")
#training_csv_human <- read.csv("/nas/longleaf/home/adaigle/work/human_TEs/phase3_data/1kg_sv_grch38/NA12878_MEI_nonref.GRCh38.csv")

#training_csv_human_ref <- training_csv_human %>% filter(DET=="DEL")

human_truth <- read.table("/nas/longleaf/home/adaigle/work/human_TEs/12878.valudated.hg38.bed")
truth <- human_truth %>%
  dplyr::rename(
    seqnames  = V1,
    start     = V2,
    end       = V3,
    TE        = V4,
    zygosity  = V5,
    strand    = V6
  ) %>%
  mutate(
    length         = end - start,
    heterozygosity = case_when(
      zygosity == "0/1" ~ 1,
      zygosity == "0/0" ~ 0.5,
      TRUE               ~ NA_real_
    )
  ) %>%
  group_by(TE) %>%
  mutate(
    ID = paste0(TE, row_number())
  ) %>%
  ungroup() %>%
  select(seqnames, start, end, ID, length, strand, TE, heterozygosity)

truth <- truth %>%
  mutate(end = ifelse(end == ".", start +1, end))
truth <- truth %>%
  mutate(strand = ifelse(strand %in% c("+", "-"), strand, "*"))

ISO1_og <- read.table("/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_inputs/human_ref_tes.bed") #%>%
ISO1_og_gr <- GRanges(
    seqnames = ISO1_og$V1,
    ranges = IRanges(start = ISO1_og$V2, end = ISO1_og$V3),
    mcols=ISO1_og$V7
    )

read_mcclintock_het_format <- function(path) {
  #this function reads in the alternative format of mcclintock files. 
  #it needs to be a separate funciton so it outputs a df to a tibble (readLines is tricky)
  #has been modified in this script to only keep calls on chrs. 2L and 2R
  lines <- readLines(path)
  data <- strsplit(lines, "\t|\\|")
  df <- data.frame(do.call(rbind, data))
  df <- df[-1,]
  #df <- df %>% filter(df$X1%in%c("3L","3R", "X"))
  return(df)
}

mcclintock_results <- tibble(
  caller = list.files(result_path)[!grepl("summary|trimgalore", list.files(result_path))],
  result_file = unlist(map(paste0(result_path, caller), ~list.files(.x, pattern = "_nonredundant.bed"))),
  fullpath = paste0(result_path, caller, "/", result_file),
  data = lapply(fullpath, read_mcclintock_het_format),
  df_list = map(data, ~split(.x, .x[["X4"]]))  # split dataframes by TE
) %>%
  unnest(cols = df_list) %>%
  mutate(
    genome = "NA12878",
    TE = map_chr(df_list, ~.x$X4[1]),
    length = map_int(df_list, ~length(.x$X1))
  )

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
  # I gave mcclintock a 1-based bed file for refs, so don't need to correct this
  # This means the programs were checking with a start +1 of where it needed to be
  reference_freq_filter = map(reference, ~ .x %>% filter(X6 >= 0.10)),
  ref_gr = map(reference_freq_filter, ~ GRanges(
    seqnames = .x$X1,
    #teflon ref coords are off by 1 for unknown reasons (mcclintock error?). Adding 1 to correct this. 
    ranges = IRanges(start = ifelse(.x$X8 == "teflon", as.numeric(.x$X2)+1, as.numeric(.x$X2)), end = as.numeric(.x$X3)),
    TE = .x$X4, heterozygosity = ifelse(.x$X6 == "NA", '1', .x$X6))))

mcclintock_results <- mcclintock_results %>% mutate(
    nonreference_freq_filter = map(nonref_gr, ~ .x %>% as.data.frame() %>% filter(heterozygosity >= 0.10) %>% makeGRangesFromDataFrame(keep.extra.columns=T)),
)

benchmark_mapping_results <- mcclintock_results %>% mutate(
    truth_forTE = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(truth %>% filter(TE==.x), keep.extra.columns = T )),
    truth_ref = map(truth_forTE,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    truth_nonref = map(truth_forTE, #extended bc a tp is within 500 bp of a tp
        ~ extend(
            subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE),
            500, 500)),
    truth_nonref_noextend = map(truth_forTE, #extended bc a tp is within 500 bp of a tp
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE)),
    #nonref_gr_filter = map2(nonref_gr_filter, TE, ~ GRanges(seqnames = seqnames(.x), 
    #    ranges = ranges(.x),
    #    TE = .y, heterozygosity=nonref_gr_filter$heterozygosity)),
   # f1_score = unlist(map2(nonreference_freq_filter, truth_nonref, ~ f1_score(.x, .y))),
    nonref_false_negatives = map2(nonreference_freq_filter, truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE, invert=TRUE)),
    nonref_true_positives = map2(nonreference_freq_filter, truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)),
    nonref_false_positives = map2(nonreference_freq_filter, truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, invert = TRUE)),
    nonref_true_positives_length = unlist(map(truth_nonref,~ length(.x))),
    nonref_calls_true_positives_length = unlist(map(nonref_true_positives,~ length(.x))),
    nonref_calls_false_positives_length = unlist(map(nonref_false_positives,~ length(.x))),
    nonref_calls_false_negatives_length = unlist(map(nonref_false_negatives,~ length(.x))),
    nonref_calls_length= unlist(map(nonreference_freq_filter,~ length(.x))),
    #stats = map2(nonref_gr_filter, truth_nonref, ~ tp_fp_fn(.x, .y))
    precision = nonref_calls_true_positives_length / (nonref_calls_true_positives_length + nonref_calls_false_positives_length),
    recall = nonref_calls_true_positives_length / (nonref_calls_true_positives_length + nonref_calls_false_negatives_length),
    f1_score = 2 * (precision * recall) / (precision + recall)
)

# pivot to long form
df_long <- benchmark_mapping_results %>%
  pivot_longer(
    cols      = c(precision, recall, f1_score),
    names_to  = "metric",
    values_to = "value"
  )

caller_order <- c("TEforest", "retroseq", "temp2", "megane_builtin", "megane_custom", "deepmei")

# 2. Apply this order to the dataframe
#    (Assuming your dataframe is named df_long)
#    We force the 'caller' column to be a factor with the specific levels.
#    Any callers in the data not in this list will become NA.
#    If you have other callers you want to keep at the end, add them to the vector above.
df_long$caller <- factor(df_long$caller, levels = caller_order)

# Define the colors
colors <- c(
  "TEforest" = "#F8766D",
  "temp2" = "#CD9600",
  "retroseq" = "#7CAE00",
  "temp" = "#00BE67",
  "popoolationte" = "#00BFC4",
  "teflon" = "#00A9FF",
  "PopoolationTE2" = "#C77CFF",
  "tepid" = "#FF61CC",
  "deepmei" = "#00A9FF",
  "megane_builtin" = "#C77CFF",
  "megane_custom" = "#FF61CC"
)

# 3. Make the plot
p <- ggplot(df_long, aes(x = metric, y = value, fill = caller)) +
  # Added color = "black" for dark borders
  # Added size = 0.3 (optional) to make the border lines slightly thinner/crisper
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.3) +
  
  scale_fill_manual(values = colors) +
  facet_wrap(~ TE, ncol = 4) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  ) +
  labs(
    title = "",
    x     = "",
    y     = "Score",
    fill  = ""
  )

# save to file
ggsave(
  filename = "/nas/longleaf/home/adaigle/TEforest/revised_plots/human_performance.png",
  plot     = p,
  width    = 12,
  height   = 8,
  dpi      = 300
)

ggsave(
  filename = "/nas/longleaf/home/adaigle/TEforest/revised_plots/human_performance.jpg",
  plot     = p,
  width    = 12,
  height   = 8,
  dpi      = 300
)

# 1) Collect per-TE counts and sum across TE per caller
combined_counts <- benchmark_mapping_results %>%
  transmute(
    caller,
    TE,
    TP = nonref_calls_true_positives_length,
    FP = nonref_calls_false_positives_length,
    FN = nonref_calls_false_negatives_length
  ) %>%
  group_by(caller) %>%
  summarise(
    TP = sum(TP, na.rm = TRUE),
    FP = sum(FP, na.rm = TRUE),
    FN = sum(FN, na.rm = TRUE),
    support = TP + FN,                   # total true nonref sites across TEs
    .groups = "drop"
  ) %>%
  mutate(
    precision = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
    recall    = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
    # micro-F1 (same as 2*P*R/(P+R), but numerically stable):
    f1_score  = if_else((2*TP + FP + FN) > 0, 2*TP / (2*TP + FP + FN), NA_real_)
  )

# 2) Long format for plotting
df_combined_long <- combined_counts %>%
  select(caller, precision, recall, f1_score) %>%
  pivot_longer(
    cols = c(precision, recall, f1_score),
    names_to = "metric",
    values_to = "value"
  )

# 3) Plot: Combined (micro-averaged) scores across all TE families
p_combined <- ggplot(df_combined_long, aes(x = metric, y = value, fill = caller)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  ) +
  labs(
    title = "Combined (Micro-averaged) Precision, Recall & F1 across TE families",
    x     = "Metric",
    y     = "Score",
    fill  = "Caller"
  )

ggsave(
  filename = "/nas/longleaf/home/adaigle/TEforest/plots/human_performance_repeatmasker_consensusonly_all.png",
  plot     = p_combined,
  width    = 12,
  height   = 8,
  dpi      = 300
)


## MEGAnE -> McClintock converter (filtered to LINE1/ALU/SVA)
## -----------------------------------------
library(utils)

## Map MEGAnE family labels to McClintock TE names
map_family_to_te <- function(fam) {
  f <- tolower(trimws(fam))
  if (grepl("line/l1", f, fixed = TRUE)) return("LINE1")
  if (grepl("sine/alu", f, fixed = TRUE)) return("ALU")
  if (grepl("retroposon/sva", f, fixed = TRUE)) return("SVA")
  return(NA_character_)
}

## Extract subfield like "subfamily_pred:status=...,MEI=...,<len>,+/-"
get_subfamily_pred_block <- function(x) {
  m <- regexpr("subfamily_pred:[^\t]+", x, perl = TRUE)
  ifelse(m > 0, regmatches(x, m), NA_character_)
}

## Strand from "+/-" token inside subfamily_pred; default "."
extract_strand <- function(s) {
  if (is.na(s)) return(".")
  m <- regexec("([+-])/([+-])", s, perl = TRUE)
  reg <- regmatches(s, m)
  if (length(reg) && length(reg[[1]]) >= 2) substr(reg[[1]][2], 1, 1) else "."
}

## Genotype number from a token like "1;PASS;..." or "2;PASS;..."
extract_geno_from_tokens <- function(row_tokens) {
  toks <- row_tokens[!is.na(row_tokens) & nzchar(row_tokens)]
  hit  <- grep("^[0-9]+;PASS;", toks, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(sub(";.*", "", hit[1])) / 2  # 1 -> 0.5 (het), 2 -> 1 (hom)
}

## Main converter
convert_megane_to_mcclintock <- function(in_bed, out_bed,
                                         genome_name = "NA12878_highcov_1",
                                         caller = "megane_builtin",
                                         rp_tag = "rp",
                                         track_label = NULL) {
  if (!file.exists(in_bed)) stop("Input file not found: ", in_bed)
  if (file.exists(out_bed)) stop("Output already exists, not overwriting: ", out_bed)

  df <- read.delim(in_bed, header = FALSE, sep = "\t",
                   quote = "", comment.char = "", fill = TRUE,
                   colClasses = "character")
  if (ncol(df) < 4) stop("Unexpected BED format: fewer than 4 columns in: ", in_bed)

  # Keep only the three requested families
  keep <- grepl("LINE/L1|SINE/Alu|Retroposon/SVA", df[[4]], ignore.case = TRUE)
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) stop("No records remain after filtering to LINE/L1, SINE/Alu, Retroposon/SVA in: ", in_bed)

  chr   <- df[[1]]
  start <- suppressWarnings(as.integer(df[[2]]))
  end   <- suppressWarnings(as.integer(df[[3]]))

  # Determine TE name from family column (4)
  te_name <- vapply(df[[4]], map_family_to_te, character(1))

  # Collapse trailing columns per row to parse flexible fields
  extra_mat <- if (ncol(df) > 4) df[, 5:ncol(df), drop = FALSE] else df[, 0, drop = FALSE]
  collapse_row <- function(v) {
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v)) paste(v, collapse = "\t") else ""
  }
  extra <- if (ncol(extra_mat)) apply(extra_mat, 1, collapse_row) else rep("", nrow(df))

  # Strand from subfamily_pred (+/-)
  subfield <- vapply(extra, get_subfamily_pred_block, character(1))
  strand   <- vapply(subfield, extract_strand, character(1))
  strand[!strand %in% c("+","-")] <- "."  # keep '.' if unknown

  # Genotype (float)
  geno <- if (ncol(extra_mat)) apply(extra_mat, 1, extract_geno_from_tokens) else rep(NA_real_, nrow(df))

  # Row numbering after filtering
  rownum <- seq_len(nrow(df))

  # Build McClintock ID
  id <- paste(
    te_name,
    "non-reference",
    ifelse(is.na(geno), "NA", as.character(geno)),
    genome_name,
    caller,
    rp_tag,
    rownum,
    sep = "|"
  )

  mcclintock <- data.frame(
    chr, start, end, id,
    score = 0,
    strand,
    stringsAsFactors = FALSE
  )

  # Prepare output
  out_dir <- dirname(out_bed)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(track_label)) {
    # Derive a readable track from caller
    base <- paste0(genome_name, "_", caller, "_nonredundant")
    track_label <- gsub("[^A-Za-z0-9_\\-]+", "_", base)
  }

  con <- file(out_bed, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(sprintf('track name="%s" description="%s"', track_label, track_label), con)
  write.table(mcclintock, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  message("Wrote: ", out_bed, "  (", nrow(mcclintock), " records)")
}

## ----------------- Run for both inputs -----------------

genome_name <- "NA12878_highcov_1"

# 1) Built-in MEGAnE
in_builtin  <- "/work/users/a/d/adaigle/human_TEs/megane_output/results_builtin/MEI_final_gaussian_genotyped.bed"
out_builtin <- "/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_runs/highcov_repeatmasker/NA12878_highcov_1/results/megane_builtin/NA12878_highcov_megane_builtin_nonredundant.bed"
convert_megane_to_mcclintock(
  in_bed      = in_builtin,
  out_bed     = out_builtin,
  genome_name = genome_name,
  caller      = "megane_builtin",
  rp_tag      = "rp",
  track_label = "NA12878_highcov_megane_builtin_nonredundant"
)

# 2) Custom MEGAnE
in_custom  <- "/work/users/a/d/adaigle/human_TEs/megane_output/results_custom/MEI_final_gaussian_genotyped.bed"
out_custom <- "/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_runs/highcov_repeatmasker/NA12878_highcov_1/results/megane_custom/NA12878_highcov_megane_custom_nonredundant.bed"
convert_megane_to_mcclintock(
  in_bed      = in_custom,
  out_bed     = out_custom,
  genome_name = genome_name,
  caller      = "megane_custom",
  rp_tag      = "rp",
  track_label = "NA12878_highcov_megane_custom_nonredundant"
)


### Create Deepmei mcclintock style file
## Map DeepMEI family (from INFO/ME second token) to McClintock name
map_family_to_te <- function(fam) {
  f <- tolower(trimws(fam))
  if (grepl("^line/l1$", f)) return("LINE1")
  if (grepl("^sine/alu$", f)) return("ALU")
  if (grepl("^retroposon/sva$", f)) return("SVA")
  return(NA_character_)
}

## Parse ME=Subfamily:Family:... from INFO
extract_me_family <- function(info_str) {
  m <- regexec("ME=([^;]+)", info_str, perl = TRUE)
  got <- regmatches(info_str, m)
  if (!length(got) || length(got[[1]]) < 2) return(NA_character_)
  parts <- strsplit(got[[1]][2], ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  parts[2]  # Family token like "SINE/Alu"
}

## Extract integer INFO values safely
extract_info_int <- function(info_str, key) {
  m <- regexec(paste0(";", key, "=(-?\\d+)"), paste0(";", info_str), perl = TRUE)
  got <- regmatches(paste0(";", info_str), m)
  if (!length(got) || length(got[[1]]) < 2) return(NA_integer_)
  as.integer(got[[1]][2])
}

## Strand from polyN_direction (1 -> '+', -1 -> '-', 0/NA -> '.')
extract_strand <- function(info_str) {
  v <- extract_info_int(info_str, "polyN_direction")
  if (is.na(v)) return(".")
  if (v > 0) return("+")
  if (v < 0) return("-")
  "."
}

## Genotype from a single-sample GT
gt_to_val <- function(gt) {
  if (is.na(gt) || gt == "." || gt == "./.") return(NA_real_)
  alleles <- unlist(strsplit(gt, "[/|]"))
  if (!length(alleles)) return(NA_real_)
  # count nonzero alleles as 'alt'
  alt_ct <- sum(suppressWarnings(as.integer(alleles)) > 0, na.rm = TRUE)
  if (alt_ct == 1) return(0.5)
  if (alt_ct >= 2) return(1.0)  # treat multi-alt as homo for our purposes
  NA_real_  # 0/0 -> NA (we only want non-reference)
}

## Choose BED interval from POS + INFO
calc_bed_coords <- function(pos, info_str) {
  l <- extract_info_int(info_str, "clipPosLeft")
  r <- extract_info_int(info_str, "clipPosRight")
  tsd <- extract_info_int(info_str, "TSD_len")
  if (!is.na(l) && !is.na(r) && !is.na(tsd) && tsd > 0 && r > l) {
    # TSD interval; DeepMEI encodes r as 1-based exclusive end -> BED end = r - 1
    start <- l - 1L
    end   <- r - 1L
  } else {
    # Fallback to 1-bp site at POS
    start <- pos - 1L
    end   <- pos
  }
  c(start, end)
}

convert_deepmei_vcf_to_mcclintock <- function(in_vcf, out_bed,
                                              genome_name = "NA12878_highcov_1",
                                              caller = "deepmei",
                                              rp_tag = "rp",
                                              track_label = NULL) {
  if (!file.exists(in_vcf)) stop("Input VCF not found: ", in_vcf)
  if (file.exists(out_bed)) stop("Output already exists, not overwriting: ", out_bed)

  # Read VCF body (skip headers)
  v <- read.delim(in_vcf, header = FALSE, sep = "\t", quote = "",
                  comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(v) < 10) stop("VCF appears not to have a FORMAT and sample column: ", in_vcf)

  colnames(v)[1:10] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

  # Keep insertions that passed filters
  v <- v[v$FILTER == "PASS", , drop = FALSE]

  # Extract family and map to TE; keep only LINE1/ALU/SVA
  fam  <- vapply(v$INFO, extract_me_family, character(1))
  te   <- vapply(fam, map_family_to_te, character(1))
  keep <- !is.na(te)
  v    <- v[keep, , drop = FALSE]; te <- te[keep]

  if (!nrow(v)) stop("No LINE1/ALU/SVA records remain after filtering in: ", in_vcf)

  # Coordinates (BED) and strand
  coords <- t(mapply(calc_bed_coords, v$POS, v$INFO))
  bed_start <- as.integer(coords[,1])
  bed_end   <- as.integer(coords[,2])
  strand    <- vapply(v$INFO, extract_strand, character(1))

  # Genotype
  # Extract GT from SAMPLE using FORMAT
  fmt_keys <- strsplit(v$FORMAT, ":", fixed = TRUE)
  smp_vals <- strsplit(v$SAMPLE, ":", fixed = TRUE)
  get_gt <- function(keys, vals) {
    idx <- match("GT", keys)
    if (is.na(idx) || idx < 1 || idx > length(vals)) return(NA_character_)
    vals[idx]
  }
  gt_str <- mapply(get_gt, fmt_keys, smp_vals, USE.NAMES = FALSE)
  geno   <- vapply(gt_str, gt_to_val, numeric(1))

  # Keep only non-reference (contains '1')
  nonref <- grepl("1", gt_str)
  v      <- v[nonref, , drop = FALSE]
  te     <- te[nonref]
  bed_start <- bed_start[nonref]; bed_end <- bed_end[nonref]
  strand <- strand[nonref]; geno <- geno[nonref]

  # Row numbering after filtering
  rownum <- seq_len(nrow(v))

  # Build McClintock ID
  id <- paste(
    te,
    "non-reference",
    ifelse(is.na(geno), "NA", as.character(geno)),
    genome_name,
    caller,
    rp_tag,
    rownum,
    sep = "|"
  )

  mcclintock <- data.frame(
    v$CHROM,
    bed_start,
    bed_end,
    id,
    score = 0,
    strand,
    stringsAsFactors = FALSE
  )

  # Output
  out_dir <- dirname(out_bed)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(track_label)) {
    base <- paste0(genome_name, "_", caller, "_nonredundant")
    track_label <- gsub("[^A-Za-z0-9_\\-]+", "_", base)
  }

  con <- file(out_bed, open = "wt"); on.exit(close(con), add = TRUE)
  writeLines(sprintf('track name="%s" description="%s"', track_label, track_label), con)
  write.table(mcclintock, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  message("Wrote: ", out_bed, "  (", nrow(mcclintock), " records)")
}

## ---- Run it ----
in_vcf  <- "/work/users/a/d/adaigle/deepmei/runs/11355357/out/DeepMEI_output/NA12878_highcov.bam/NA12878_highcov.bam.vcf"  # adjust path if needed
out_bed <- "/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_runs/highcov_repeatmasker/NA12878_highcov_1/results/deepmei/NA12878_highcov_deepmei_nonredundant.bed"

convert_deepmei_vcf_to_mcclintock(
  in_vcf    = in_vcf,
  out_bed   = out_bed,
  genome_name = "NA12878_highcov_1",
  caller      = "deepmei",
  rp_tag      = "rp",
  track_label = "NA12878_highcov_deepmei_nonredundant"
)