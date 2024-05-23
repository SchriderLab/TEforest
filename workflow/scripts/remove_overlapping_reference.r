#!/usr/bin/env Rscript

# the purpose of this script is to remove any reference TEs that overlap with other 
library(tidyverse)
library(GenomicRanges)

args <- commandArgs(TRUE)
print(args)
refinput <- args[1]

basedir_outputs_path <- getwd()

output_path <- paste0(basedir_outputs_path, "/ref_genome/")
dir.create(file.path(output_path))

#ref <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/DeNovoCoordinates/ISO1.bed")
#find any regions that do not overlap with other regions

ref <- read.table(refinput)

ref_gr <- GRanges(
    seqnames = ref$V1,
    ranges = IRanges(start = ref$V2, end = ref$V3),
    mcols=ref$V7
    )

ref_gr <- makeGRangesFromDataFrame(ref,
    seqnames.field = "V1",
    start.field = "V2", 
    end.field = "V3",
    keep.extra.columns = T)

ref_gr_edit <- GRanges(
    seqnames = ref$V1,
    ranges = IRanges(start = ref$V2 +1, end = ref$V3),
    mcols=ref$V7
    )

nested_tes <- subsetByOverlaps(ref_gr,ref_gr_edit,type="within")
nested_filter <- subsetByOverlaps(ref_gr,nested_tes,invert=T)

reduced_ref <- GenomicRanges::reduce(nested_filter)

overlap_filter <- subsetByOverlaps(ref_gr, reduced_ref, type="equal")

#looking for nested TEs, changing coordinates one bp to make sure coords aren't identical




filtered_ref <- as.data.frame(overlap_filter)
filtered_ref$width <- NULL
filtered_ref$strand <- NULL
write.table(filtered_ref, file= paste0(output_path, "filtered_ISO1.bed"), row.name=F,col.name=F, quote=F)
