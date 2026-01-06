library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
library(viridis)
library(scales)
library(cowplot)

# basepath_nolengthfilter <- "/nas/longleaf/home/adaigle/work/test_TEforest/nolengthfilter_"
basepath_basenorm_feats <- "/nas/longleaf/home/adaigle/users/test_teforest/hetrefte_"

# Modified function to accept a vector of genome names and sum them
read_df <- function(basepath, covg, genomes) {
  path <- paste0(basepath, covg, "X_otherscript/2L_2R_plots/performance_plots/")
  
  # 1. Read all genomes into a list
  df_list <- lapply(genomes, function(g) {
    d <- readRDS(paste0(path, g, "_ref.rds"))
    return(d)
  })
  
  # 2. Combine, Group, and Sum
  # This assumes all Non-Numeric columns are identifiers (e.g., 'caller', 'family')
  # and all Numeric columns are data to be summed.
  df <- bind_rows(df_list) %>%
    group_by(across(where(negate(is.numeric)))) %>% # Group by all text/factor columns
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
  
  # 3. Proceed with standard processing
  df$covg <- covg
  df <- df %>%
    mutate(
      caller = case_when(
        caller == "TEforest_classifier_filter_bps" ~ "TEforest",
        caller == "temp" ~ "TEMP",
        caller == "temp2" ~ "TEMP2",
        caller == "teflon" ~ "TEFLoN",
        caller == "retroseq" ~ "RetroSeq",
        caller == "popoolationte" ~ "PopoolationTE",
        caller == "popoolationte2" ~ "PopoolationTE2",
        TRUE ~ caller 
      )
    ) %>%
    filter(caller %in% c("TEforest", "TEMP2", "TEFLoN", "PopoolationTE2"))
  
  return(df)
}

# Updated read_df_ref to match the new logic (in case it is needed)
read_df_ref <- function(basepath, covg, genomes) {
  path <- paste0(basepath, covg, "X_otherscript/2L_2R_plots/performance_plots/")
  
  df <- readRDS(paste0(path, genomes[1], "_ref.rds"))
  
  if (length(genomes) > 1) {
    for (g in genomes[-1]) {
      df_next <- readRDS(paste0(path, g, "_ref.rds"))
      numeric_cols <- sapply(df, is.numeric)
      df[numeric_cols] <- df[numeric_cols] + df_next[numeric_cols]
    }
  }
  
  df$covg <- covg
  df <- df %>%
    mutate(
      caller = case_when(
        caller == "TEforest_classifier_filter_bps" ~ "TEforest",
        caller == "temp" ~ "TEMP",
        caller == "temp2" ~ "TEMP2",
        caller == "teflon" ~ "TEFLoN",
        caller == "retroseq" ~ "RetroSeq",
        caller == "popoolationte" ~ "PopoolationTE",
        caller == "popoolationte2" ~ "PopoolationTE2",
        TRUE ~ caller 
      )
    ) %>%
    filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))
  
  return(df)
}

colors <- c(
  "TEforest" = "#F8766D",
  "TEMP2" = "#CD9600",
  "RetroSeq" = "#7CAE00",
  "TEMP" = "#00BE67",
  "PopoolationTE" = "#00BFC4",
  "TEFLoN" = "#00A9FF",
  "PopoolationTE2" = "#C77CFF",
  "tepid" = "#FF61CC",
  "TEforest_nolengthfilter" = "black"
)

plot_coverage <- function(df) {
  # Reorder the levels of the caller factor to match the order in colors
  df$caller <- factor(df$caller, levels = names(colors))
  cov <- ggplot(df, aes(x = covg, y = Score, color = caller)) +
    geom_line(size = 1.5, alpha = 0.5) +
    geom_point(size = 3.5) +
    facet_wrap(~Metric) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50), labels = c(5, 10, 20, 30, 40, 50)) +
    labs(
      title = "",
      x = "Coverage",
      y = "Value"
    ) +
    theme_minimal() +
    scale_color_manual(values = colors) +
    theme(
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 18)
    )
  return(cov)
}

# Function to plot each metric separately
plot_coverage_threeplot <- function(df) {
  df$caller <- factor(df$caller, levels = names(colors))
  metrics <- unique(df$Metric)
  plot_list <- list()
  
  for (metric in metrics) {
    d <- df[df$Metric == metric, ]
    d_top    <- d[d$caller == "TEforest", ]
    d_bottom <- d[d$caller != "TEforest", ]
    
    # callers present in this panel, with TEforest first
    present <- intersect(names(colors), unique(as.character(d$caller)))
    present <- c("TEforest", setdiff(present, "TEforest"))
    p <- ggplot() +
      # draw others first (underneath)
      geom_line( data = d_bottom, aes(covg, Score, color = caller), size = 1.5, alpha = 0.5) +
      geom_point(data = d_bottom, aes(covg, Score, color = caller), size = 3.5) +
      # draw TEforest last (on top)
      geom_line( data = d_top, aes(covg, Score, color = caller), size = 1.5, alpha = 0.5) +
      geom_point(data = d_top, aes(covg, Score, color = caller), size = 3.5) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(breaks = c(5,10,20,30,40,50)) +
      labs(x = "Coverage", y = metric) +
      theme_minimal() +
      scale_color_manual(
        values = colors[present],
        breaks = present,
        limits = present       # ensures both order and inclusion == present
      ) +
      theme(
        legend.title = element_blank(),
        legend.text  = element_text(size = 16),
        axis.text    = element_text(size = 12),
        axis.title   = element_text(size = 12),
        strip.text   = element_text(size = 18)
      )
    
    plot_list[[metric]] <- p
  }
  plot_list
}


### EXECUTION BLOCK ###

coverage_values <- c(5, 10, 20, 30, 40, 50)

# 1. 54bp Data: Combining original A2_A3 with ISO-154_A2
basenorm_feats_54bp <- do.call(rbind, lapply(coverage_values, function(covg) {
  read_df(basepath_basenorm_feats, covg, c("A2_A3", "ISO-154_A2"))
}))
basenorm_feats_54bp$Score <- basenorm_feats_54bp$Score / 2
basenorm_feats_54bp_plot <- basenorm_feats_54bp %>% select(caller, covg, Metric, Score)


# 2. 125bp Data: Combining original AKA-017_GIM-024 with ISO-1125_AKA-017
basenorm_feats_125bp <- do.call(rbind, lapply(coverage_values, function(covg) {
  read_df(basepath_basenorm_feats, covg, c("AKA-017_GIM-024", "ISO-1125_AKA-017"))
}))
basenorm_feats_125bp$Score <- basenorm_feats_125bp$Score / 2
basenorm_feats_125bp_plot <- basenorm_feats_125bp %>% select(caller, covg, Metric, Score)


# 3. 150bp (154) Data: Combining original JUT-008_MUN-009 with ISO-1_JUT-008
basenorm_feats_150bp <- do.call(rbind, lapply(coverage_values, function(covg) {
  read_df(basepath_basenorm_feats, covg, c("JUT-008_MUN-009", "ISO-1_JUT-008"))
}))
basenorm_feats_150bp$Score <- basenorm_feats_150bp$Score / 2
basenorm_feats_150bp_plot <- basenorm_feats_150bp %>% select(caller, covg, Metric, Score)

performance_150 <- plot_coverage_threeplot(basenorm_feats_150bp_plot)

performance_125 <- plot_coverage_threeplot(basenorm_feats_125bp_plot)

performance_54 <- plot_coverage_threeplot(basenorm_feats_54bp_plot)

merged_performance <- ggarrange(
    performance_150$f1_score, performance_125$f1_score, performance_54$f1_score, 
    performance_150$precision, performance_125$precision, performance_54$precision, 
    performance_150$recall, performance_125$recall, performance_54$recall, 
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    ncol = 3, nrow = 3, common.legend = T, legend="bottom"
  )

titles_row <- ggarrange(
  text_grob("151 bp", face = "bold", size = 14, just = "center"),
  text_grob("125 bp", face = "bold", size = 14, just = "center"),
  text_grob("54 bp",  face = "bold", size = 14, just = "center"),
  ncol = 3
)
final_with_col_labels <- ggarrange(
  titles_row, merged_performance,
  ncol = 1, heights = c(0.06, 1)  # tweak the header height as you like
)

ggsave("/nas/longleaf/home/adaigle/TEforest/revised_plots/merged_performance_reference.svg", final_with_col_labels, dpi=300, width = 8.5, height = 9)
ggsave("/nas/longleaf/home/adaigle/TEforest/revised_plots/merged_performance_reference.jpg", final_with_col_labels, dpi=300, width = 8.5, height = 9)

read_genotype_data <- function(basepath, covg, genome) {
    # Initialize an empty list to store all loaded objects
    loaded_data <- list()
    
    # Iterate over the covg values (c(5, 30))
    for (coverage in covg) {
        # Construct the path based on basepath, coverage, and genome names
        path <- paste0(basepath, coverage, "X_otherscript/2L_2R_plots/")
        
        # Load each of the RDS files for the given coverage
        genotype_data <- readRDS(paste0(path, genome, "_genotypes.rds"))
        teforest_conf_matrix <- readRDS(paste0(path, genome, "_TEforest_confusion_matrix.rds"))
        retroseq_conf_matrix <- readRDS(paste0(path, genome, "_retroseq_confusion_matrix.rds"))
        temp2_conf_matrix <- readRDS(paste0(path, genome, "_temp2_confusion_matrix.rds"))
        freqplot_data <- readRDS(paste0(path, genome, "_freqplot_data.rds"))
        
        # Add coverage information to the genotype data frame
        genotype_data$covg <- coverage
        
        # Apply mutations and filtering on genotype data frame
        genotype_data <- genotype_data %>%
            mutate(
                caller = case_when(
                    caller == "TEforest_classifier_filter_bps" ~ "TEforest",
                    caller == "temp" ~ "TEMP",
                    caller == "temp2" ~ "TEMP2",
                    caller == "teflon" ~ "TEFLoN",
                    caller == "retroseq" ~ "RetroSeq",
                    caller == "popoolationte" ~ "PopoolationTE",
                    caller == "popoolationte2" ~ "PopoolationTE2",
                    TRUE ~ caller  # Keep the original value for other cases
                )
            ) %>%
            filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))
        
        # Store all loaded data in the list
        loaded_data[[paste0("genotype_data_", coverage, "X")]] <- genotype_data
        loaded_data[[paste0("TEforest_conf_matrix_", coverage, "X")]] <- teforest_conf_matrix
        loaded_data[[paste0("retroseq_conf_matrix_", coverage, "X")]] <- retroseq_conf_matrix
        loaded_data[[paste0("temp2_conf_matrix_", coverage, "X")]] <- temp2_conf_matrix
        loaded_data[[paste0("freqplot_data_", coverage, "X")]] <- freqplot_data
    }
    
    # Return the list of loaded data
    return(loaded_data)
}

read_genotype_data <- function(basepath, covg_list, genomes) {
    loaded_data <- list()
    
    for (coverage in covg_list) {
        path <- paste0(basepath, coverage, "X_otherscript/2L_2R_plots/")
        
        # Inner function to read and clean a single genome's genotype file
        read_single_genome <- function(g_name) {
            df <- readRDS(paste0(path, g_name, "_genotypes.rds"))
            df <- df %>%
                mutate(
                    caller = case_when(
                        caller == "TEforest_classifier_filter_bps" ~ "TEforest",
                        caller == "temp" ~ "TEMP",
                        caller == "temp2" ~ "TEMP2",
                        caller == "teflon" ~ "TEFLoN",
                        caller == "retroseq" ~ "RetroSeq",
                        caller == "popoolationte" ~ "PopoolationTE",
                        caller == "popoolationte2" ~ "PopoolationTE2",
                        TRUE ~ caller
                    )
                ) %>%
                filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))
            return(df)
        }
        
        # Read all genomes provided (e.g., the pair)
        list_of_dfs <- lapply(genomes, read_single_genome)
        
        # Bind them together and summarize
        # Numeric columns (F1 score) get averaged ((A+B)/2)
        # List columns (tphomo/tphet) get concatenated so unlist() works on the full set later
        combined_df <- bind_rows(list_of_dfs) %>%
            group_by(caller) %>%
            summarise(
                # Average the simple numeric metrics (F1, Precision, etc)
                across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                # Concatenate the list columns (tphomo, tphet) containing the raw data
                across(where(is.list), ~ list(unlist(.x))),
                .groups = "drop"
            ) %>%
            mutate(covg = coverage)
        
        # Store in list
        loaded_data[[paste0("genotype_data_", coverage, "X")]] <- combined_df
    }
    
    return(loaded_data)
}

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
library(viridis)
library(scales)
library(cowplot)

basepath_basenorm_feats <- "/nas/longleaf/home/adaigle/users/test_teforest/hetrefte_"

colors <- c(
    "TEforest" = "#F8766D",
    "TEMP2" = "#CD9600",
    "RetroSeq" = "#7CAE00",
    "TEMP" = "#00BE67",
    "PopoolationTE" = "#00BFC4",
    "TEFLoN" = "#00A9FF",
    "PopoolationTE2" = "#C77CFF",
    "tepid" = "#FF61CC"
)

# Modified to accept a vector of genomes, combine them, and remove unused matrix/freq files
read_genotype_data <- function(basepath, covg_list, genomes) {
    loaded_data <- list()
    
    for (coverage in covg_list) {
        path <- paste0(basepath, coverage, "X_otherscript/2L_2R_plots/")
        
        # Inner function to read and clean a single genome's genotype file
        read_single_genome <- function(g_name) {
            df <- readRDS(paste0(path, g_name, "_genotypes.rds"))
            df <- df %>%
                mutate(
                    caller = case_when(
                        caller == "TEforest_classifier_filter_bps" ~ "TEforest",
                        caller == "temp" ~ "TEMP",
                        caller == "temp2" ~ "TEMP2",
                        caller == "teflon" ~ "TEFLoN",
                        caller == "retroseq" ~ "RetroSeq",
                        caller == "popoolationte" ~ "PopoolationTE",
                        caller == "popoolationte2" ~ "PopoolationTE2",
                        TRUE ~ caller
                    )
                ) %>%
                filter(caller %in% c("TEforest", "TEMP2", "TEFLoN", "PopoolationTE2"))
            return(df)
        }
        
        # Read all genomes provided (e.g., the pair)
        list_of_dfs <- lapply(genomes, read_single_genome)
        
        # Bind them together and summarize
        # Numeric columns (F1 score) get averaged ((A+B)/2)
        # List columns (tphomo/tphet) get concatenated so unlist() works on the full set later
        combined_df <- bind_rows(list_of_dfs) %>%
            group_by(caller) %>%
            summarise(
                # Average the simple numeric metrics (F1, Precision, etc)
                across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                # Concatenate the list columns (tphomo, tphet) containing the raw data
                across(where(is.list), ~ list(unlist(.x))),
                .groups = "drop"
            ) %>%
            mutate(covg = coverage)
        
        # Store in list
        loaded_data[[paste0("genotype_data_", coverage, "X")]] <- combined_df
    }
    
    return(loaded_data)
}

plot_genotyping <- function(df) {
    # Reorder the levels of the caller factor to match the order in colors
    df$caller <- factor(df$caller, levels = names(colors))
    df <- df[order(df$caller, decreasing = TRUE), ]
    cov <- ggplot(df, aes(x = covg, y = macro_f1_score, color = caller)) +
        geom_line(size = 1.5, alpha = 0.5) +
        geom_point(size = 3.5, alpha=0.8) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50), labels = c(5, 10, 20, 30, 40, 50)) +
        labs(
            title = "",
            x = "Coverage",
            y = "f1_score"
        ) +
        theme_minimal() +
        scale_color_manual(values = colors) +
        theme(
            legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 18)
        )
    return(cov)
}

plot_combined <- function(df) {
    set.seed(123)  # For reproducibility
    df <- df %>%
        mutate(jittered_covg = jitter(covg, amount = 1.2)) 
    df$caller <- factor(df$caller, levels = names(colors))
    df <- df[order(df$caller, decreasing = TRUE), ]
    
    ggplot(df) +
        # Plot for Tphomo
        geom_line(aes(x = jittered_covg, y = tphomo_mean, color = caller), size = 1, alpha = 0.8) +
        geom_point(aes(x = jittered_covg, y = tphomo_mean, color = caller), size = 3) +
        # Plot for Tphet with lighter color
        geom_line(aes(x = jittered_covg, y = tphet_mean, color = caller), size = 1, alpha=0.5) +
        geom_point(aes(x = jittered_covg, y = tphet_mean, color = caller), size = 3, alpha=0.6, shape = 17) +
        labs(
            x = "Coverage",
            y = "Predicted prevalence",
            title = ""
        ) +
        theme_minimal() +
        scale_color_manual(values = colors) +
        theme(
            legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 18)
        )
}

combine_genotype_data <- function(all_data, covg) {
    combined_genotype_data <- bind_rows(
        lapply(covg, function(c) all_data[[paste0("genotype_data_", c, "X")]])
    )
    return(combined_genotype_data)
}


### EXECUTION ###

covg <- c(5, 10, 20, 30, 40, 50)

# 150bp (154): Combine JUT-008_MUN-009 and ISO-1_JUT-008
all_data <- read_genotype_data(basepath_basenorm_feats, covg, c("JUT-008_MUN-009", "ISO-1_JUT-008"))

# 125bp: Combine AKA-017_GIM-024 and ISO-1125_AKA-017
all_data_125 <- read_genotype_data(basepath_basenorm_feats, covg, c("AKA-017_GIM-024", "ISO-1125_AKA-017"))

# 54bp: Combine A2_A3 and ISO-154_A2
all_data_54 <- read_genotype_data(basepath_basenorm_feats, covg, c("A2_A3", "ISO-154_A2"))


# Combine across coverages
combined_genotype_data <- combine_genotype_data(all_data, covg)
combined_genotype_data_125 <- combine_genotype_data(all_data_125, covg)
combined_genotype_data_54 <- combine_genotype_data(all_data_54, covg)

# Process Statistics for Heterozygosity Plots
# (This uses unlist() on the concatenated lists created in read_genotype_data)
combined_genotype_data <- combined_genotype_data %>%
  rowwise() %>%
  mutate(
    tphomo_mean = mean(unlist(tphomo), na.rm = TRUE),
    tphomo_sd = sd(unlist(tphomo), na.rm = TRUE),
    tphet_mean = mean(unlist(tphet), na.rm = TRUE),
    tphet_sd = sd(unlist(tphet), na.rm = TRUE)
  )

combined_genotype_data_125 <- combined_genotype_data_125 %>%
  rowwise() %>%
  mutate(
    tphomo_mean = mean(unlist(tphomo), na.rm = TRUE),
    tphomo_sd = sd(unlist(tphomo), na.rm = TRUE),
    tphet_mean = mean(unlist(tphet), na.rm = TRUE),
    tphet_sd = sd(unlist(tphet), na.rm = TRUE)
  )

combined_genotype_data_54 <- combined_genotype_data_54 %>%
  rowwise() %>%
  mutate(
    tphomo_mean = mean(unlist(tphomo), na.rm = TRUE),
    tphomo_sd = sd(unlist(tphomo), na.rm = TRUE),
    tphet_mean = mean(unlist(tphet), na.rm = TRUE),
    tphet_sd = sd(unlist(tphet), na.rm = TRUE)
  )

# Generate Plots
alldata_geno <- plot_genotyping(combined_genotype_data)
alldata_geno_125 <- plot_genotyping(combined_genotype_data_125)
alldata_geno_54 <- plot_genotyping(combined_genotype_data_54)

combined_plot <- plot_combined(combined_genotype_data)
combined_plot_125 <- plot_combined(combined_genotype_data_125)
combined_plot_54 <- plot_combined(combined_genotype_data_54)

# Extract Legends
# Extract the legend specifically from panel I for each coverage
combined_plot_legend_125 <- ggpubr::get_legend(
  combined_plot_125 +
    theme(legend.position = "bottom",
          legend.title = element_blank())
)

# Build Final Figure
genotyping_figure <- ggpubr::ggarrange(
  alldata_geno + theme(legend.position = "none"),
  alldata_geno_125 + theme(legend.position = "none"),
  alldata_geno_54 + theme(legend.position = "none"),
  combined_plot + theme(legend.position = "none"), 
  combined_plot_125 + theme(legend.position = "none"),
  combined_plot_54 + theme(legend.position = "none"),
  
  labels = c("A","B","C","D","E","F"),
  ncol = 3, nrow = 2,
  legend = "bottom",
  legend.grob = combined_plot_legend_125
)

finalgeno <- ggarrange(
  titles_row, genotyping_figure,
  ncol = 1, heights = c(0.06, 1)  # tweak the header height as you like
)
# Save
ggsave("/nas/longleaf/home/adaigle/TEforest/revised_plots/genotyping_ref.svg", finalgeno, dpi=300, width = 8.5, height = 8.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/revised_plots/genotyping_ref.jpg", finalgeno, dpi=300, width = 8.5, height = 8.5)