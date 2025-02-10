#!/usr/bin/env Rscript
# Load required packages
library(GenomicAlignments)
library(tidyverse)
library(GenomicRanges)

args <- commandArgs(TRUE)
print(args)
genome <- args[1]
euchromatin_coordinates_path <- args[2]
ISO1_og_path <- args[3]
#basedir_outputs_path <- args[2]

basedir_outputs_path <- getwd()
#genome <- "RW7000"
#basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/test_celegans2"
#euchromatin_coordinates_path <- "/nas/longleaf/home/adaigle/work/celegans/chroms.txt"
#ISO1_og_path <- "/nas/longleaf/home/adaigle/work/celegans/tc1_improved.bed"

output_path <- paste0(basedir_outputs_path, "/candidate_regions_data/", genome)
featvec_csv_path <- paste0(basedir_outputs_path, "/featvec_csvs")
dir.create(file.path(featvec_csv_path))


ISO1_og <- read.table(ISO1_og_path)
ISO1_og$V7 <- gsub("-", "_", ISO1_og$V7)

euchromatin_coordinates <- makeGRangesFromDataFrame(read.table(euchromatin_coordinates_path), seqnames.field="V1", start.field="V2", end.field="V3")

#nonreference_genome_truth <- read.table(paste0("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/", genome, "_Ref_Coord.bed"))
#colnames(nonreference_genome_truth) <- c("seqnames", "start", "end", "ID", "length", "strand", "TE")
#nonreference_genome_truth$TE <- gsub("-", "_", nonreference_genome_truth$TE)


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
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

ISO1_og_gr <- GRanges(
    seqnames = ISO1_og$V1,
    ranges = IRanges(start = ISO1_og$V2, end = ISO1_og$V3),
    mcols=ISO1_og$V7
    )

ref_IDs <- ISO1_og$V4

ISO1_og$heterozygosity <- 0

reference_df_labled <- ISO1_og %>% arrange(V1) %>%
  rename(Chrom = V1, Ref_begin = V2, Ref_end = V3, TE= V7, Class=heterozygosity) %>%
  #expand candidate region 50 bp, convert from 1-based to 0-based
  mutate(Sample = genome, Ref_begin = Ref_begin - 51, Ref_end = Ref_end + 50) %>% 
  select(Sample, Chrom, Ref_begin, Ref_end, Class, TE)

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
            ranges = IRanges(start = start(.x) + 200, end = end(.x) - 200))),
    expand = map2(expand, TE, ~ GRanges(seqnames = seqnames(.x), 
            ranges = ranges(.x),
            TE = .y))
)
print("filtering complete")

calls <- do.call(c, c(mapping_results_filter$expand))

calls_df <- calls %>%
  as.data.frame() %>%
  mutate(Class = 0, Sample = genome, start = start - 1) %>%
  select(Sample, seqnames, start, end, Class, TE) %>%
  rename(Chrom = seqnames, Ref_begin = start, Ref_end = end)

write.csv(calls_df, file = paste0(featvec_csv_path, "/", genome, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(reference_df_labled, file = paste0(featvec_csv_path, "/", genome, "_reference.csv"), quote = FALSE, row.names = FALSE)

print(paste(genome, "csv written!"))
