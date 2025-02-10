#!/usr/bin/env Rscript
# the purpose of this script is to take the ouputs of the random forest
# and convert to a mcclintock style BED file
# currently outputs one with and without breakpoint calls
library(GenomicRanges)
library(tidyverse)
library(GenomicAlignments)

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

extend <- function(x, upstream=0, downstream=0) {
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}


args <- commandArgs(TRUE)
print(args)
genomes <- as.list(args[1])
csv_path <- as.list(args[2])
results_path <- as.list(args[3])
aligned_dir <- as.list(args[4])
csv_path_reference <- as.list(args[5])
candidate_regions_data_dir <- as.list(args[6])

basedir_outputs_path <- getwd()
csv_path <- paste0(getwd(), "/", csv_path)
print(csv_path)
csv_path_reference <- paste0(getwd(), "/", csv_path_reference)
print(csv_path_reference)

#genomes <- "AKA-017_GIM-024"
#csv_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X_redo2/2L2R_classifer_reference/AKA-017_GIM-024.csv"
#results_path <- "2L2R_model_validation_output_classifier/"
#aligned_dir <- "aligned_het/"
#basedir_outputs_path <- "/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X_redo2/"
#csv_path_reference <- "/work/users/a/d/adaigle/test_TEforest/basenorm_feats_50X_redo2/2L2R_classifer_reference/AKA-017_GIM-024.csv"
#candidate_regions_data_dir <- "/candidate_regions_data/"

apply_labels_reference <- function (genome) {
# this function also returns ref coords to correct length
genome_name_length <- length(strsplit(genome, "-")[[1]])
# Set the variable 'names' based on the length
if (genome_name_length == 1) {
    names <- c("genome", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path_reference) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(start=as.numeric(start)+51, end = as.numeric(end)-50, class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%
    GRanges()
} else if (genome_name_length == 2) {
    names <- c("genome1", "genome2", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path_reference) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    mutate(start=as.numeric(start)+51, end = as.numeric(end)-50, class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges()
} else if (genome_name_length == 3) {
    names <- c("genome1", "genome2", "genome3", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path_reference) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  genome3, sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    mutate(start=as.numeric(start)+51, end = as.numeric(end)-50, class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges() 
} else {
  names <- "Error"
}
tp_calls <- tp_calls[,1:2]


#labels <- GRanges(read.table(labels_path, header=TRUE))

#pred here assumes truth labels are 2x the frequeuncy
tp_labeled <- as.data.frame(tp_calls) %>%
    mutate(TE_string=paste(TE, "reference", pred*0.5, genome, "TEforest", "rp", row_number(), sep = "|"))

mcclintock_format_df <- tp_labeled[c(1:3, 8)]
mcclintock_format_df$score <- 0
mcclintock_format_df$strand <- "."

#results_path <- paste0(basedir_outputs_path, "/output/", genome, "/")

#dir.create(file.path(results_path))
#write.table(mcclintock_format_df, 
#    file = paste0(results_path, genome, "_TEforest_nonredundant_reference.bed"),
#    quote = F, sep = "\t", row.names = F
#)

return(mcclintock_format_df)
}


formatted_reference <- lapply(genomes, apply_labels_reference)[[1]]


apply_labels <- function (genome) {
#csv_path <- paste0(basedir_outputs_path, "/", outdir, "/", genome, "/predictions.csv")


#if a genome has hyphens in its name it needs to be split up differently
#if it exceeds three hyphens this will break
#TODO make genome names have no hyphens, or make this code more flexible
genome_name_length <- length(strsplit(genome, "-")[[1]])
# Set the variable 'names' based on the length
if (genome_name_length == 1) {
    names <- c("genome", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%
    GRanges() 
} else if (genome_name_length == 2) {
    names <- c("genome1", "genome2", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    mutate(class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges()
} else if (genome_name_length == 3) {
    names <- c("genome1", "genome2", "genome3", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  genome3, sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    mutate(class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges() 
} else {
  names <- "Error"
}
tp_calls <- tp_calls[,1:2]


#labels <- GRanges(read.table(labels_path, header=TRUE))

#pred here assumes truth labels are 2x the frequeuncy
tp_labeled <- as.data.frame(tp_calls) %>%
    mutate(TE_string=paste(TE, "non-reference", pred*0.5, genome, "TEforest", "rp", row_number(), sep = "|"))

mcclintock_format_df <- tp_labeled[c(1:3, 8)]
mcclintock_format_df$score <- 0
mcclintock_format_df$strand <- "."

#results_path <- paste0(basedir_outputs_path, "/output/", genome, "/")
mcclintock_format_df <- rbind(mcclintock_format_df, formatted_reference) %>%
    arrange(seqnames,start)

dir.create(file.path(results_path))
write.table(mcclintock_format_df, 
    file = paste0(results_path, genome, "_TEforest_nonredundant.bed"),
    quote = F, sep = "\t", row.names = F
)

#return(mcclintock_format_df)
}


lapply(genomes, apply_labels)


find_breakpoint <- function(candidate_region, genome) {
    #print(candidate_region)
# shrink candidate region down to avoid spurious split reads
#based on original size of region
#shrink_size <- (width(candidate_region) / 1.5)/2
#shrink_size <- 0
shrink_size <- 200
candidate_region <- shrink(candidate_region, shrink_size, shrink_size)
tesbp <- ScanBamParam(which = IRangesList(c(candidate_region)), what = c("qname", "seq"))
#bamfile <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/synthetic_hets_mcclintock/", genome, "_1/intermediate/mapped_reads/", genome, "_1.sorted.bam")
bamfile <- paste0(basedir_outputs_path, "/", aligned_dir, genome, ".bam")
tebam <- BamFile(bamfile)
tereads <- readGAlignments(tebam, param = tesbp)

# Regular expression to match cigars ending in "S" or starting with a digit followed by "S"
# Goal is to grab truly soft clipped reads rather than short inserts
pattern <- "^(\\d+S)|.*(S)$"

softclips <- subset(tereads, grepl(pattern, cigar))
#check for hardclipped reads, none in this dataset
#hardclips <- subset(tereads, grepl("H", cigar))

split_lengths <- gsub("[^0-9]+", "", softclips@cigar)
non_split_lengths <- gsub("([0-9]+)S.*", "\\1", softclips@cigar)
split_seq <- data.frame(qname = softclips@elementMetadata[[1]], 
                         rname = softclips@seqnames,
                         start = start(softclips),
                         end = end(softclips),
                         strand = strand(softclips),
                         split = substr(softclips@elementMetadata[[2]], 1, as.integer(non_split_lengths)), # non-split part of sequence
                         nonsplit = substr(softclips@elementMetadata[[2]], as.integer(non_split_lengths)+1, nchar(softclips@elementMetadata[[2]]))) # split part of sequence

split_seq <- split_seq[complete.cases(split_seq), ] #remove rows with NA values

#deals with cases where M comes before S in split read
#notice different in row naming-- important 
split_lengths2 <- gsub("[^0-9]+", "", softclips@cigar)
non_split_lengths2 <- gsub("([0-9]+)M.*", "\\1", softclips@cigar)
#edited from original, see comments to right
split_seq2 <- data.frame(qname = softclips@elementMetadata[[1]], 
                         rname = softclips@seqnames,
                         start = start(softclips),
                         end = end(softclips), 
                         strand = strand(softclips),
                         nonsplit = substr(softclips@elementMetadata[[2]], 1, as.integer(non_split_lengths2)), #changed this bc te is on other side!
                         split = substr(softclips@elementMetadata[[2]], as.integer(non_split_lengths2)+1, nchar(softclips@elementMetadata[[2]]))) #changed this bc te is on other side!

split_seq2 <- split_seq2[complete.cases(split_seq2), ] #remove rows with NA values

#merge two cases together
split_seq_final <- rbind(split_seq, split_seq2)

#filter rows where the shortest subsequence is less than 10 bp
# Create a new column for the length of the shortest sequence
split_seq_final$shortest_seq <- apply(split_seq_final[, c("split", "nonsplit")], 1, function(x) min(nchar(x)))
# Filter the r table to remove rows where the shortest sequence is less than 5
split_seq_final <- subset(split_seq_final, shortest_seq >= 5)

if (nrow(split_seq_final) == 0) {
    #print(candidate_region)
    #print(shrink_size)
    return()  # Exit the loop
  }

#use supporting reads to filter Galignments object
softclips_filter <- softclips[softclips@elementMetadata[[1]] %in% split_seq_final$qname]

cigar_ops <- strsplit(as.character(softclips_filter@cigar), "(?<=\\D)(?=\\d)", perl=TRUE) 
cigar_first <- lapply(cigar_ops, function(x)  gsub("[0-9]+", "", x[1]))


# Extract either the start or end position based on the first cigar operation
# If the first string is M, the end coordinate is the
supported_breakpoints <- c()
supported_breakpoints <- ifelse(cigar_first == "S", 
    start(softclips_filter), end(softclips_filter))
seqname <- as.character(seqnames(softclips[1]))

  # Check if supported_breakpoints is NULL
  if (is.null(supported_breakpoints[1])) {
    #print("PURGED")
    #print(candidate_region)
    return()  # Exit the loop
  }

bps <- table(supported_breakpoints)
#taking top 2 rn, can be modified to only pick one if it is clearly larger
# also could think of a rule to select close bps if two are far from each other?
# just need to find good examples...
top_breakpoints <- sort(as.numeric(names(sort(bps, decreasing = TRUE)[1:2])))

if (length(top_breakpoints) == 1) {
    gr <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = top_breakpoints[1], end = top_breakpoints[1]),
    TE = candidate_region$TE[[1]],
    pred = candidate_region$pred[[1]]
    )
    return(gr)  # Exit the loop
  }


gr <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = top_breakpoints[1], end = top_breakpoints[2]),
    TE = candidate_region$TE[[1]],
    pred = candidate_region$pred[[1]]
    )

return(gr)
}

#candidate_region <- split(tp_calls)[[46]]

find_breakpoint_tefirstversion <- function(candidate_region, genome) {
    #print(candidate_region)
# shrink candidate region down to avoid spurious split reads
#based on original size of region
#shrink_size <- (width(candidate_region) / 1.5)/2
#shrink_size <- 0
TE_specific_bam <- paste0(basedir_outputs_path, candidate_regions_data_dir, genome, "/", candidate_region$TE, "_to_ISO1.bam")

shrink_size <- 200
candidate_region <- shrink(candidate_region, shrink_size, shrink_size)
tesbp <- ScanBamParam(which = IRangesList(c(candidate_region)), what = c("qname", "seq"))
#bamfile <- paste0("/nas/longleaf/home/adaigle/work/mcclintock_stuff/synthetic_hets_mcclintock/", genome, "_1/intermediate/mapped_reads/", genome, "_1.sorted.bam")
bamfile <- paste0(basedir_outputs_path, "/", aligned_dir, genome, ".bam")

tebam <- BamFile(TE_specific_bam)
tereads <- readGAlignments(tebam, param = tesbp)

# Regular expression to match cigars ending in "S" or starting with a digit followed by "S"
# Goal is to grab truly soft clipped reads rather than short inserts
pattern <- "^(\\d+S)|.*(S)$"

softclips <- subset(tereads, grepl(pattern, cigar))
#check for hardclipped reads, none in this dataset
#hardclips <- subset(tereads, grepl("H", cigar))

split_lengths <- gsub("[^0-9]+", "", softclips@cigar)
non_split_lengths <- gsub("([0-9]+)S.*", "\\1", softclips@cigar)
split_seq <- data.frame(qname = softclips@elementMetadata[[1]], 
                         rname = softclips@seqnames,
                         start = start(softclips),
                         end = end(softclips),
                         strand = strand(softclips),
                         split = substr(softclips@elementMetadata[[2]], 1, as.integer(non_split_lengths)), # non-split part of sequence
                         nonsplit = substr(softclips@elementMetadata[[2]], as.integer(non_split_lengths)+1, nchar(softclips@elementMetadata[[2]]))) # split part of sequence

split_seq <- split_seq[complete.cases(split_seq), ] #remove rows with NA values

#deals with cases where M comes before S in split read
#notice different in row naming-- important 
split_lengths2 <- gsub("[^0-9]+", "", softclips@cigar)
non_split_lengths2 <- gsub("([0-9]+)M.*", "\\1", softclips@cigar)
#edited from original, see comments to right
split_seq2 <- data.frame(qname = softclips@elementMetadata[[1]], 
                         rname = softclips@seqnames,
                         start = start(softclips),
                         end = end(softclips), 
                         strand = strand(softclips),
                         nonsplit = substr(softclips@elementMetadata[[2]], 1, as.integer(non_split_lengths2)), #changed this bc te is on other side!
                         split = substr(softclips@elementMetadata[[2]], as.integer(non_split_lengths2)+1, nchar(softclips@elementMetadata[[2]]))) #changed this bc te is on other side!

split_seq2 <- split_seq2[complete.cases(split_seq2), ] #remove rows with NA values

#merge two cases together
split_seq_final <- rbind(split_seq, split_seq2)

#filter rows where the shortest subsequence is less than 10 bp
# Create a new column for the length of the shortest sequence
split_seq_final$shortest_seq <- apply(split_seq_final[, c("split", "nonsplit")], 1, function(x) min(nchar(x)))
# Filter the r table to remove rows where the shortest sequence is less than 5
split_seq_final <- subset(split_seq_final, shortest_seq >= 5)

if (nrow(split_seq_final) == 0) {
    #print(candidate_region)
    #print(shrink_size)
    return()  # Exit the loop
  }

#use supporting reads to filter Galignments object
softclips_filter <- softclips[softclips@elementMetadata[[1]] %in% split_seq_final$qname]

cigar_ops <- strsplit(as.character(softclips_filter@cigar), "(?<=\\D)(?=\\d)", perl=TRUE) 
cigar_first <- lapply(cigar_ops, function(x)  gsub("[0-9]+", "", x[1]))


# Extract either the start or end position based on the first cigar operation
# If the first string is M, the end coordinate is the
supported_breakpoints <- c()
supported_breakpoints <- ifelse(cigar_first == "S", 
    start(softclips_filter), end(softclips_filter))
seqname <- as.character(seqnames(softclips[1]))

  # Check if supported_breakpoints is NULL
  if (is.null(supported_breakpoints[1])) {
    #print("PURGED")
    #print(candidate_region)
    return()  # Exit the loop
  }

bps <- table(supported_breakpoints)
#taking top 2 rn, can be modified to only pick one if it is clearly larger
# also could think of a rule to select close bps if two are far from each other?
# just need to find good examples...
top_breakpoints <- sort(as.numeric(names(sort(bps, decreasing = TRUE)[1:2])))

if (length(top_breakpoints) == 0) {
print(candidate_region)
print("no te specific breakpoints, looking for split reads")
tebam <- BamFile(bamfile)
tereads <- readGAlignments(tebam, param = tesbp)

# Regular expression to match cigars ending in "S" or starting with a digit followed by "S"
# Goal is to grab truly soft clipped reads rather than short inserts
pattern <- "^(\\d+S)|.*(S)$"

softclips <- subset(tereads, grepl(pattern, cigar))
#check for hardclipped reads, none in this dataset
#hardclips <- subset(tereads, grepl("H", cigar))

split_lengths <- gsub("[^0-9]+", "", softclips@cigar)
non_split_lengths <- gsub("([0-9]+)S.*", "\\1", softclips@cigar)
split_seq <- data.frame(qname = softclips@elementMetadata[[1]], 
                         rname = softclips@seqnames,
                         start = start(softclips),
                         end = end(softclips),
                         strand = strand(softclips),
                         split = substr(softclips@elementMetadata[[2]], 1, as.integer(non_split_lengths)), # non-split part of sequence
                         nonsplit = substr(softclips@elementMetadata[[2]], as.integer(non_split_lengths)+1, nchar(softclips@elementMetadata[[2]]))) # split part of sequence

split_seq <- split_seq[complete.cases(split_seq), ] #remove rows with NA values

#deals with cases where M comes before S in split read
#notice different in row naming-- important 
split_lengths2 <- gsub("[^0-9]+", "", softclips@cigar)
non_split_lengths2 <- gsub("([0-9]+)M.*", "\\1", softclips@cigar)
#edited from original, see comments to right
split_seq2 <- data.frame(qname = softclips@elementMetadata[[1]], 
                         rname = softclips@seqnames,
                         start = start(softclips),
                         end = end(softclips), 
                         strand = strand(softclips),
                         nonsplit = substr(softclips@elementMetadata[[2]], 1, as.integer(non_split_lengths2)), #changed this bc te is on other side!
                         split = substr(softclips@elementMetadata[[2]], as.integer(non_split_lengths2)+1, nchar(softclips@elementMetadata[[2]]))) #changed this bc te is on other side!

split_seq2 <- split_seq2[complete.cases(split_seq2), ] #remove rows with NA values

#merge two cases together
split_seq_final <- rbind(split_seq, split_seq2)

#filter rows where the shortest subsequence is less than 10 bp
# Create a new column for the length of the shortest sequence
split_seq_final$shortest_seq <- apply(split_seq_final[, c("split", "nonsplit")], 1, function(x) min(nchar(x)))
# Filter the r table to remove rows where the shortest sequence is less than 5
split_seq_final <- subset(split_seq_final, shortest_seq >= 5)

if (nrow(split_seq_final) == 0) {
    #print(candidate_region)
    #print(shrink_size)
    return()  # Exit the loop
  }

#use supporting reads to filter Galignments object
softclips_filter <- softclips[softclips@elementMetadata[[1]] %in% split_seq_final$qname]

cigar_ops <- strsplit(as.character(softclips_filter@cigar), "(?<=\\D)(?=\\d)", perl=TRUE) 
cigar_first <- lapply(cigar_ops, function(x)  gsub("[0-9]+", "", x[1]))


# Extract either the start or end position based on the first cigar operation
# If the first string is M, the end coordinate is the
supported_breakpoints <- c()
supported_breakpoints <- ifelse(cigar_first == "S", 
    start(softclips_filter), end(softclips_filter))
seqname <- as.character(seqnames(softclips[1]))

  # Check if supported_breakpoints is NULL
  if (is.null(supported_breakpoints[1])) {
    #print("PURGED")
    #print(candidate_region)
    return()  # Exit the loop
  }

bps <- table(supported_breakpoints)
#taking top 2 rn, can be modified to only pick one if it is clearly larger
# also could think of a rule to select close bps if two are far from each other?
# just need to find good examples...
top_breakpoints <- sort(as.numeric(names(sort(bps, decreasing = TRUE)[1:2])))

}
if (length(top_breakpoints) == 1) {
    gr <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = top_breakpoints[1], end = top_breakpoints[1]),
    TE = candidate_region$TE[[1]],
    pred = candidate_region$pred[[1]]
    )
    return(gr)  # Exit the loop
  }


gr <- GRanges(
    seqnames = seqname,
    ranges = IRanges(start = top_breakpoints[1], end = top_breakpoints[2]),
    TE = candidate_region$TE[[1]],
    pred = candidate_region$pred[[1]]
    )

return(gr)
}

apply_labels_breakpoints <- function (genome) {
#if a genome has hyphens in its name it needs to be split up differently
#if it exceeds three hyphens this will break
#TODO make genome names have no hyphens, or make this code more flexible
genome_name_length <- length(strsplit(genome, "-")[[1]])
# Set the variable 'names' based on the length
if (genome_name_length == 1) {
    names <- c("genome", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges() 
} else if (genome_name_length == 2) {
    names <- c("genome1", "genome2", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    mutate(class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges()
} else if (genome_name_length == 3) {
    names <- c("genome1", "genome2", "genome3", "seqnames", "start", "end", "class", "TE")
    tp_calls <- read.csv(csv_path) %>%
    filter(pred != 0) %>% # remove 0 entries
    separate_wider_delim(file,names = names, delim = "-") %>%
    mutate(genome = paste(genome1,  genome2,  genome3, sep = "-"),
             genome1 = NULL, genome2 = NULL, genome3 = NULL) %>% 
    mutate(class = NULL,true = NULL, cntrl_score = NULL, genome = NULL) %>%    
    GRanges() 
} else {
  names <- "Error"
}
tp_calls <- tp_calls[,1:2]

mapping_results_filter_breakpoints <- lapply(split(tp_calls), function(candidate_region) find_breakpoint_tefirstversion(candidate_region, genome))

print("breakpoint finding complete!")
remove_null_elements <- function(granges_list) {
    granges_list[!sapply(granges_list, is.null)]
}
mapping_results_filter_breakpoints <- remove_null_elements(mapping_results_filter_breakpoints)

mapping_results_filter_breakpoints <- do.call(c, GRangesList(mapping_results_filter_breakpoints))

add_mc_format <- as.data.frame(mapping_results_filter_breakpoints) %>%
    mutate(TE_string=paste(TE, "non-reference", pred*0.5, genome, "TEforest", "rp", row_number(), sep = "|"))

mcclintock_format_df <- add_mc_format[c(1:3, 8)]
mcclintock_format_df$score <- 0
mcclintock_format_df$strand <- "."
#results_path <- paste0(basedir_outputs_path, "/output/", genome, "/")

mcclintock_format_df <- rbind(mcclintock_format_df, formatted_reference) %>%
    arrange(seqnames,start)

dir.create(file.path(results_path))
write.table(mcclintock_format_df, 
    file = paste0(results_path, genome, "_TEforest_bps_nonredundant.bed"),
    quote = F, sep = "\t", row.names = F
)

}

lapply(genomes, apply_labels_breakpoints)
