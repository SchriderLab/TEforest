#!/usr/bin/env Rscript
# Load required packages
library(GenomicAlignments)
library(tidyverse)
library(GenomicRanges)

# collect command line args in a file
# in order genome1 name,base dir path (ouputs)
args <- commandArgs(TRUE)
print(args)
genome1 <- args[1]
genome2 <- args[2]
euchromatin_coordinates_path <- args[3]

basedir_outputs_path <- getwd()

#genome1 <- "A2"
#genome2 <- "A3"
#basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/A2_improved_consensus_revertedfiltering"

combine_id <- paste0(genome1,"_", genome2)
output_path <- paste0(basedir_outputs_path, "/candidate_regions_data_het/", combine_id)
featvec_csv_path <- paste0(basedir_outputs_path, "/featvec_csvs_hettrain")

dir.create(file.path(featvec_csv_path))

#need to replace this path once I fix this... 
ISO1_og <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/ISO-1_Ref_Coord.bed")
ISO1_og$V7 <- gsub("-", "_", ISO1_og$V7)

euchromatin_coordinates <- makeGRangesFromDataFrame(read.table("/nas/longleaf/home/adaigle/work/mcclintock_stuff/euchromatin.txt"), seqnames.field="V1", start.field="V2", end.field="V3")
#euchromatin_coordinates <- makeGRangesFromDataFrame(read.table(euchromatin_coordinates_path), seqnames.field="V1", start.field="V2", end.field="V3")

nonreference_genome1_truth <- read.table(paste0("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/", genome1, "_Ref_Coord.bed"))
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
#test <- distanceToNearest(makeGRangesFromDataFrame(truth))
#find indices of truth that are within 0 bp of another call.
#zeros <- as.data.frame(test) %>% filter(distance==0) 
#truth[as.integer(unlist(zeros[1])),] %>% arrange(seqnames, start) %>% mutate(width=end-start)

# version of truth without heterozygosity info, so I can run the pipeline normally
nonreference_genome_truth <- truth %>% select(seqnames, start, end, ID, length, strand, TE)
heterozygosity_2R3R_onlyhets <- truth %>% filter(heterozygosity==0.5) %>% select(seqnames, start, end, ID, length, strand, TE)
heterozygosity_2R3R_nohets <- truth %>% filter(heterozygosity==1.0) %>% select(seqnames, start, end, ID, length, strand, TE)

#functions
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

mapping_results <- tibble(fullpath = list.files(path = output_path, pattern="to_ISO1.bam$", full.names = TRUE),
    TE = str_extract(fullpath, "\\w+(?=_to_ISO1)") , #?= is a positive lookahead assertion
    bam = lapply(fullpath, readGAlignments),
    granges = map(bam, ~ GenomicRanges::reduce(GRanges(seqnames = seqnames(.x),
                                 ranges = ranges(.x),
                                 strand = strand(.x)), ignore.strand = TRUE))
  )

mapping_results <- mapping_results %>% mutate(
    split_reads = map(bam, ~ .x[grepl("S", .x@cigar)]),
    split_reads_ranges = map(split_reads, ~ GenomicRanges::reduce(GRanges(seqnames = seqnames(.x),
                                 ranges = ranges(.x),
                                 strand = strand(.x)), ignore.strand = TRUE)),
    regions_with_split_reads = map2(granges, split_reads_ranges, ~ subsetByOverlaps(.x, .y))
)
print("mapping results loaded")


#ISO1_filtered <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/DeNovoCoordinates/ISO1.bed")
ISO1_og_gr <- GRanges(
    seqnames = ISO1_og$V1,
    ranges = IRanges(start = ISO1_og$V2, end = ISO1_og$V3),
    mcols=ISO1_og$V7
    )

#family specific reference TEs

mapping_results_filter <- mapping_results %>% mutate(
    reference_seqs_df = map(TE, ~ ISO1_og %>% filter(V7==.x)),
    reference_seqs_granges = map(reference_seqs_df, ~ GRanges(seqnames = .x$V1,
                                 ranges = IRanges(start = .x$V2, end = .x$V3),
                                 ignore.strand = TRUE)),
    euchromatin_calls = map(granges, ~ subsetByOverlaps(.x, euchromatin_coordinates, ignore.strand = TRUE)),
    euchromatin_calls_not_inrefte_expand = map(euchromatin_calls, # toggle numbers in ranges to change expansion
        ~ GenomicRanges::reduce(GRanges(seqnames = seqnames(.x), 
            ranges = IRanges(start = start(.x) - 200, end = end(.x) + 200)))), #expand to capture deleted ref tes
    #euchromatin_calls_touching_called_refte = map2(euchromatin_calls_not_inrefte, reference_seqs_granges, 
    #    ~ subsetByOverlaps(.x,.y,type="within",ignore.strand = TRUE)),
    euchromatin_calls_wsplitreads = map(regions_with_split_reads, ~ subsetByOverlaps(.x, euchromatin_coordinates, ignore.strand = TRUE)),
    ref_tes_sametype_within_candidate_regions = map2(reference_seqs_granges, euchromatin_calls_not_inrefte_expand,  
        ~ subsetByOverlaps(.x,.y,type="within",ignore.strand = TRUE)),
    euchromatin_calls_without_nested_reftes = map2(euchromatin_calls_wsplitreads, ref_tes_sametype_within_candidate_regions,  
        ~ subsetByOverlaps(.x,.y,ignore.strand = TRUE, invert=T)),
    euchromatin_calls_inrefte = map(euchromatin_calls_without_nested_reftes, 
        ~ subsetByOverlaps(.x,ISO1_og_gr,type="within",ignore.strand = TRUE)),
    euchromatin_calls_not_inrefte = map2(euchromatin_calls_without_nested_reftes, euchromatin_calls_inrefte, 
        ~ GenomicRanges::setdiff(.x,.y)),
    #euchromatin_calls_not_touching_called_refte = map2(euchromatin_calls_not_inrefte, euchromatin_calls_touching_called_refte, 
    #    ~ GenomicRanges::setdiff(.x,.y)),
    #euchromatin_calls_without_nested_reftes_shrunk = map(euchromatin_calls_without_nested_reftes, # toggle numbers in ranges to change expansion
    #    ~ GRanges(seqnames = seqnames(.x), 
    #        ranges = IRanges(start = start(.x) + 200, end = end(.x) - 200))),
    early_width_filter = map(euchromatin_calls_without_nested_reftes, 
        ~ subset(.x, width(.x) >= 154)), 
    expand1 = map(early_width_filter, # toggle numbers in ranges to change expansion
        ~ GenomicRanges::reduce(GRanges(seqnames = seqnames(.x), 
            ranges = IRanges(start = start(.x) - 400, end = end(.x) + 400)))),
    expand = map(expand1, # toggle numbers in ranges to change expansion
        ~ GRanges(seqnames = seqnames(.x), 
            ranges = IRanges(start = start(.x) + 200, end = end(.x) - 200)))
)
print("filtering complete")

benchmark_mapping_results <- mapping_results_filter %>% mutate(
    nonreference_genome_truth_forTE = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(nonreference_genome_truth %>% filter(TE==.x), keep.extra.columns = T )),
    nonreference_genome_truth_ref = map(nonreference_genome_truth_forTE,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    nonreference_genome_truth_nonref = map(nonreference_genome_truth_forTE,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE)),
    nonreference_genome_truth_forTE_het = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(heterozygosity_2R3R_onlyhets %>% filter(TE==.x), keep.extra.columns = T )),
    nonreference_genome_truth_ref_het = map(nonreference_genome_truth_forTE_het,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    nonreference_genome_truth_nonref_het = map(nonreference_genome_truth_forTE_het,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE)),
    nonreference_genome_truth_forTE_nohet = map(TE, # reduce to get rid of nested TEs of same type
        ~ makeGRangesFromDataFrame(heterozygosity_2R3R_nohets %>% filter(TE==.x), keep.extra.columns = T )),
    nonreference_genome_truth_ref_nohet = map(nonreference_genome_truth_forTE_nohet,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE)),
    nonreference_genome_truth_nonref_nohet = map(nonreference_genome_truth_forTE_nohet,
        ~ subsetByOverlaps(.x, ISO1_og_gr, type=c("equal"), ignore.strand=TRUE, invert = TRUE)),
    #nonreference_genome_truth_nonref = map(nonreference_genome_truth_nonref,
    #    ~ extend(.x, 500,500)),
    expand = map2(expand, TE, ~ GRanges(seqnames = seqnames(.x), 
        ranges = ranges(.x),
        TE = .y)),
    f1_score = unlist(map2(expand, nonreference_genome_truth_forTE, ~ f1_score(.x, .y))),
    nonref_false_negatives = map2(expand, nonreference_genome_truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE, invert=TRUE)),
    nonref_true_positives_true_coords = map2(expand, nonreference_genome_truth_nonref, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE)),
    nonref_true_positives = map2(expand, nonreference_genome_truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)),
    nonref_false_positives = map2(expand, nonreference_genome_truth_nonref, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE, invert = TRUE)),
    nonref_true_positives_true_coords_het = map2(expand, nonreference_genome_truth_nonref_het, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE)),
    nonref_true_positives_het = map2(expand, nonreference_genome_truth_nonref_het, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)),
    nonref_true_positives_true_coords_nohet = map2(expand, nonreference_genome_truth_nonref_nohet, 
        ~ subsetByOverlaps(.y, .x, ignore.strand = TRUE)),
    nonref_true_positives_nohet = map2(expand, nonreference_genome_truth_nonref_nohet, 
        ~ subsetByOverlaps(.x, .y, ignore.strand = TRUE)),
    nonref_true_positives_length = unlist(map(nonreference_genome_truth_nonref,~ length(.x))),
    stats = map2(expand, nonreference_genome_truth_nonref, ~ tp_fp_fn(.x, .y))
) %>%
  unnest_wider(stats)
print("benchmarking complete")

#take all my true coordinates from the benchmarking DF
#not necessary to reduce because any candidate range overlapping this dataset
#will count as a TP
false_negatives <- do.call(c, c(benchmark_mapping_results$nonref_false_negatives))
true_positives_true_coords <- do.call(c, c(benchmark_mapping_results$nonref_true_positives_true_coords))

true_positives <- do.call(c, c(benchmark_mapping_results$nonref_true_positives))
false_positives <- do.call(c, c(benchmark_mapping_results$nonref_false_positives))

true_positives_het <- do.call(c, c(benchmark_mapping_results$nonref_true_positives_het))
true_positives_nohet <- do.call(c, c(benchmark_mapping_results$nonref_true_positives_nohet))

truth_gr <- makeGRangesFromDataFrame(truth,keep.extra.columns = T )
subsetByOverlaps(truth_gr, true_positives)
subsetByOverlaps(truth_gr, true_positives_true_coords)
fn_hets <- subsetByOverlaps(truth_gr, false_negatives) %>% as.data.frame
table(fn_hets$heterozygosity)
fn_hets_2R3R <- fn_hets %>% filter(seqnames %in% c("2R", "3R"))
table(fn_hets_2R3R$heterozygosity)
table(fn_hets_2R3R$grange1)
table(fn_hets_2R3R$grange2)

tp_hets <- subsetByOverlaps(truth_gr, true_positives_true_coords) %>% as.data.frame
table(tp_hets$heterozygosity)
tp_hets_2R3R <- tp_hets %>% filter(seqnames %in% c("2R", "3R"))
table(tp_hets_2R3R$heterozygosity)
table(tp_hets_2R3R$grange1)
table(tp_hets_2R3R$grange2)

# filter preds > 2k bp from data

true_positives_het <- true_positives_het[width(true_positives_het)<2000]
true_positives_nohet <- true_positives_nohet[width(true_positives_nohet)<2000]
false_positives <- false_positives[width(false_positives)<2000]

#Here I make csv's for neural net training
# I also subtract 1 from beginning to make coordinates 0 based
tp_df_het <- true_positives_het %>%
  as.data.frame() %>%
  mutate(Class = 1, Sample = combine_id, start = start - 1) %>%
  select(Sample, seqnames, start, end, Class, TE) %>%
  rename(Chrom = seqnames, Ref_begin = start, Ref_end = end)

tp_df_nohet <- true_positives_nohet %>%
  as.data.frame() %>%
  mutate(Class = 2, Sample = combine_id, start = start - 1) %>%
  select(Sample, seqnames, start, end, Class, TE) %>%
  rename(Chrom = seqnames, Ref_begin = start, Ref_end = end)

fp_df <- false_positives %>%
  as.data.frame() %>%
  mutate(Class = 0, Sample = combine_id, start = start - 1) %>%
  select(Sample, seqnames, start, end, Class, TE) %>%
  rename(Chrom = seqnames, Ref_begin = start, Ref_end = end)

#use this code to do quick tests with tremolo outputs
#read.table("/nas/longleaf/home/adaigle/work/tremolo_output3/POS_TE_OUTSIDER_ON_REF.bed")
#tremolo_outsider <- makeGRangesFromDataFrame(read.table("/nas/longleaf/home/adaigle/work/tremolo_output3/POS_TE_OUTSIDER_ON_REF.bed"), seqnames.field="V1", start.field="V2", end.field="V3", keep.extra.columns = T)
#subsetByOverlaps(tremolo_outsider, true_positives)
#subsetByOverlaps(true_positives, tremolo_outsider)
#subsetByOverlaps(tremolo_outsider, false_positives)
#subsetByOverlaps(false_positives, tremolo_outsider)
##read.table("/nas/longleaf/home/adaigle/work/tremolo_output3/POS_TE_INSIDER_ON_REF.bed")
#tremolo_insider <- makeGRangesFromDataFrame(read.table("/nas/longleaf/home/adaigle/work/tremolo_output3/POS_TE_INSIDER_ON_REF.bed"), seqnames.field="V1", start.field="V2", end.field="V3", keep.extra.columns = T)
#subsetByOverlaps(tremolo_insider, true_positives)
#subsetByOverlaps(true_positives, tremolo_insider)
#subsetByOverlaps(tremolo_insider, false_positives)
#subsetByOverlaps(false_positives, tremolo_insider)
#
#tremolo <- c(tremolo_outsider, tremolo_insider)
#subsetByOverlaps(true_positives, tremolo)
#subsetByOverlaps(extend(false_negatives,50,50), tremolo)

featvec_csv <- rbind(tp_df_het, tp_df_nohet, fp_df)
write.csv(featvec_csv, file = paste0(featvec_csv_path, "/", combine_id, ".csv"), quote = FALSE, row.names = FALSE)
print(paste(combine_id, "csv written!"))


#telabel_df <- rbind(as.data.frame(true_positives), as.data.frame(false_positives)) %>%
#    arrange(seqnames,start) %>% 
#    mutate(start = start - 1) %>%
#    mutate(TE= as.character(TE))
#    #might be better to do this after neural net filtering...
#    #mutate(TE_string=paste(TE, "non-reference", "NA", "A4_1", "allseeingeye", "rp", row_number(), sep = "|"),
#
##using TE labels path from the input
#
#write.table(telabel_df, file = paste0(telabels_path, "/", genome, "_labels.txt"), quote = FALSE, row.names = FALSE)
#print(paste(genome, "ID table written!"))

