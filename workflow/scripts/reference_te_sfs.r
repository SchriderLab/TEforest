library(tidyverse)
library(GenomicRanges)

ISO1 <- read.table("/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates/ISO-1_Ref_Coord.bed") #%>%

dir_path <- "/nas/longleaf/home/adaigle/Rech_updated_supplemental/ReferenceCoordinates"

# List all files in the directory
files <- list.files(path = dir_path, pattern = "_Ref_Coord.bed$", full.names = TRUE)

# Function to read a file and return a list with file name and data frame
read_bed_file <- function(file_path) {
  # Extract the base name without the suffix
  base_name <- sub("_Ref_Coord.bed$", "", basename(file_path))
  # Read the BED file into a data frame
  bed_df <- read.table(file_path, header = FALSE)
  # Return a list with the name and the data frame
  list(name = base_name, file_path = file_path, df = bed_df)
}

# Apply the function to all files and combine the results into a data frame
dfs_list <- lapply(files, read_bed_file)

# Create a data frame of data frames
dfs_df <- tibble(
  name = sapply(dfs_list, function(x) x$name),
  file_path = sapply(dfs_list, function(x) x$file_path),
  df = I(sapply(dfs_list, function(x) x$df, simplify = FALSE)),
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
  labs(title = "Histogram of Reference ID Overlaps",
       x = "Number of Overlaps",
       y = "Frequency")

unique_ISO1 <- ISO1_ids[!(ISO1_ids %in% reference_ids_vector)]

(length(unique_ISO1) + 550) / (length(ref_id_counts$frequency) + length(unique_ISO1)) 


unique_ISO1_coords <- ISO1 %>% filter(V4 %in% unique_ISO1)
