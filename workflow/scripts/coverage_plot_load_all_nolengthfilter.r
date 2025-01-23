library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
library(viridis)
library(scales)
library(cowplot)

# coverage_values <- c(5, 10, 20, 30, 40, 50)

basepath_nolengthfilter <- "/nas/longleaf/home/adaigle/work/test_TEforest/nolengthfilter_"
basepath_basenorm_feats <- "/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_"

read_df <- function(basepath, covg, genome) {
    path <- paste0(basepath, covg, "X/2L_2R_plots/performance_plots/")
    df <- readRDS(paste0(path, genome, ".rds"))
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
                TRUE ~ caller # Keep the original value for other cases
            )
        ) %>%
        filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))

    return(df)
}

read_df_trim <- function(basepath, covg, genome) {
    path <- paste0(basepath, "/2L_2R_plots/performance_plots/")
    print(path)
    df <- readRDS(paste0(path, genome, ".rds"))
    df$covg <- covg
    print(df$caller)
    df <- df %>%
        mutate(
            caller = case_when(
                caller == "TEforest_trimmed" ~ "TEforest",
                caller == "temp" ~ "TEMP",
                caller == "temp2" ~ "TEMP2",
                caller == "teflon" ~ "TEFLoN",
                caller == "retroseq" ~ "RetroSeq",
                caller == "popoolationte" ~ "PopoolationTE",
                caller == "popoolationte2" ~ "PopoolationTE2",
                TRUE ~ caller # Keep the original value for other cases
            )
        ) %>%
        filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))

    return(df)
}

read_df_ref <- function(basepath, covg, genome) {
    path <- paste0(basepath, covg, "X/2L_2R_plots/performance_plots/")
    df <- readRDS(paste0(path, genome, "_ref.rds"))
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
                TRUE ~ caller # Keep the original value for other cases
            )
        ) %>%
        filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))

    return(df)
}

#coverage_values <- c(5, 10, 20, 30, 40, 50)
## Use lapply to read the data frames for each coverage and combine them into a single data frame
#nolengthfilter_54bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_nolengthfilter, covg, "A2_A3")))
#nolengthfilter_54bp_plot <- nolengthfilter_54bp %>% select(caller, covg, Metric, Score)
#
#nolengthfilter_125bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_nolengthfilter, covg, "AKA-017_GIM-024")))
#nolengthfilter_125bp_plot <- nolengthfilter_125bp %>% select(caller, covg, Metric, Score)
#
#nolengthfilter_150bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_nolengthfilter, covg, "JUT-008_MUN-009")))
#nolengthfilter_150bp_plot <- nolengthfilter_150bp %>% select(caller, covg, Metric, Score)

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
    # Reorder the levels of the caller factor to match the order in colors
    df$caller <- factor(df$caller, levels = names(colors))
    
    # Create a list of unique metrics
    metrics <- unique(df$Metric)
    
    # Initialize an empty list to store the plots
    plot_list <- list()
    
    # Loop over each metric and generate separate plots
    for (metric in metrics) {
        # Filter the data for the current metric
        df_metric <- df[df$Metric == metric, ]
        
        # Create a plot for the current metric
        p <- ggplot(df_metric, aes(x = covg, y = Score, color = caller)) +
            geom_line(size = 1.5, alpha = 0.5) +
            geom_point(size = 3.5) +
            scale_y_continuous(limits = c(0, 1)) +
            scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50), labels = c(5, 10, 20, 30, 40, 50)) +
            labs(
                title = "",
                x = "Coverage",
                y = metric  # Use the current Metric as the y-axis label
            ) +
            theme_minimal() +
            scale_color_manual(values = colors) +
            theme(
                axis.title.y = element_text(size = 12), # Add the y-axis label
                legend.title = element_blank(),
                legend.text = element_text(size = 16),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 12),
                strip.text = element_text(size = 18)
            )
        
        # Add the plot to the list, with the metric name as the key
        plot_list[[metric]] <- p
    }
    
    # Return the list of plots
    return(plot_list)
}

#plot_coverage(nolengthfilter_54bp_plot)
#plot_coverage(nolengthfilter_125bp_plot)
#plot_coverage(nolengthfilter_150bp_plot)

coverage_values <- c(5,10,20, 30, 40, 50)
# Use lapply to read the data frames for each coverage and combine them into a single data frame
basenorm_feats_54bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_basenorm_feats, covg, "A2_A3")))
basenorm_feats_54bp_plot <- basenorm_feats_54bp %>% select(caller, covg, Metric, Score)
#teforest_nolength_54 <- nolengthfilter_54bp_plot %>% filter(caller=="TEforest") %>%
#  mutate(caller = ifelse(caller == "TEforest", "TEforest_nolengthfilter", caller))
#basenorm_feats_54bp_plot2 <- rbind(teforest_nolength_54, basenorm_feats_54bp_plot)

basenorm_feats_125bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_basenorm_feats, covg, "AKA-017_GIM-024")))
basenorm_feats_125bp_plot <- basenorm_feats_125bp %>% select(caller, covg, Metric, Score)
#teforest_nolength_125 <- nolengthfilter_125bp_plot %>% filter(caller=="TEforest") %>%
#  mutate(caller = ifelse(caller == "TEforest", "TEforest_nolengthfilter", caller))
#basenorm_feats_125bp_plot2 <- rbind(teforest_nolength_125, basenorm_feats_125bp_plot)

basenorm_feats_150bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_basenorm_feats, covg, "JUT-008_MUN-009")))
basenorm_feats_150bp_plot <- basenorm_feats_150bp %>% select(caller, covg, Metric, Score)
#teforest_nolength_150 <- nolengthfilter_150bp_plot %>% filter(caller=="TEforest") %>%
#  mutate(caller = ifelse(caller == "TEforest", "TEforest_nolengthfilter", caller))
#basenorm_feats_150bp_plot2 <- rbind(teforest_nolength_150, basenorm_feats_150bp_plot)

#plot_coverage(basenorm_feats_54bp_plot2)
#plot_coverage(basenorm_feats_125bp_plot2)
#plot_coverage(basenorm_feats_150bp_plot2)

coverage_values <- c(5, 10, 20, 30, 40, 50)
basepath_homozygote_ref <- "/nas/longleaf/home/adaigle/work/test_TEforest/homozygote_ref_allbp_"
# Use lapply to read the data frames for each coverage and combine them into a single data frame
homozygote_ref_54bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df_ref(basepath_homozygote_ref, covg, "A2_A2")))
homozygote_ref_54bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df_ref(basepath_homozygote_ref, covg, "A4_A4")))

homozygote_ref_54bp_plot <- homozygote_ref_54bp %>% select(caller, covg, Metric, Score)

homozygote_ref_125bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df_ref(basepath_homozygote_ref, covg, "AKA-017_AKA-017")))
homozygote_ref_125bp_plot <- homozygote_ref_125bp %>% select(caller, covg, Metric, Score)

homozygote_ref_150bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df_ref(basepath_homozygote_ref, covg, "JUT-008_JUT-008")))
homozygote_ref_150bp_plot <- homozygote_ref_150bp %>% select(caller, covg, Metric, Score)

ref_54 <- plot_coverage(homozygote_ref_54bp_plot %>% filter(!(caller %in% c("RetroSeq", "TEMP2", "PopoolationTE"))))
ref_125 <-plot_coverage(homozygote_ref_125bp_plot %>% filter(!(caller %in% c("RetroSeq", "TEMP2", "PopoolationTE"))))
ref_150 <-plot_coverage(homozygote_ref_150bp_plot %>% filter(!(caller %in% c("RetroSeq", "TEMP2", "PopoolationTE"))))

reference_figure <- ggarrange(
    ref_54, ref_125, ref_150, 
    labels = c("A", "B", "C"),
    ncol = 1, nrow = 3, common.legend = T, legend="bottom"
  )

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/reference_figure.png", reference_figure, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/reference_figure.svg", reference_figure, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/reference_figure.tiff", reference_figure, dpi=300, width = 8.5, height = 10.5)

#norm_feats_125bp %>% select(caller, covg, Metric, Score)

basenorm_feats_150bp <- do.call(rbind, lapply(coverage_values, function(covg) read_df(basepath_basenorm_feats, covg, "JUT-008_MUN-009")))
basenorm_feats_150bp_plot <- basenorm_feats_150bp %>% select(caller, covg, Metric, Score)
plot_coverage(basenorm_feats_54bp_plot)
plot_coverage(basenorm_feats_125bp_plot)
plot_coverage(basenorm_feats_150bp_plot)

#load('/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X/2L_2R_plots/JUT-008_MUN-009.RData')


read_commoncaller_csv <- function(basepath, covg, genome, num_callers) {
    path <- paste0(basepath, covg, "X/2L_2R_plots/")
    df <- read.csv(paste0(path, genome, "_", num_callers, "_common_callers.txt"))
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
                TRUE ~ caller # Keep the original value for other cases
            )
        ) %>%
        filter(caller %in% c("TEforest", "TEMP", "TEMP2", "TEFLoN", "RetroSeq", "PopoolationTE", "PopoolationTE2"))

    return(df)
}

plot_bp <- function(df) {
    # Reorder the levels of the caller factor to match the order in colors
    df$caller <- factor(df$caller, levels = names(colors))
    cov <- ggplot(df, aes(x = covg, y = distance_mean, color = caller)) +
        geom_line(size = 1.5, alpha = 0.5) +
        geom_point(size = 3.5) +
        # Adding error bars for distance_sd
        #geom_errorbar(
        #  aes(ymin = distance_mean - distance_sd, ymax = distance_mean + distance_sd),
        #  width = 0.2,  # Adjust width as needed for clarity
        #  size = 0.7,   # Adjust line thickness of error bars
        #  alpha = 0.7   # Adjust transparency
        #) +
        #facet_wrap(~Metric) +
        #scale_y_continuous(limits = c(0, 70)) +
        scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50), labels = c(5, 10, 20, 30, 40, 50)) +
        labs(
            title = "",
            x = "Coverage",
            y = "Distance from breakpoint (bp)"
        ) +
        theme_minimal() +
        scale_color_manual(values = colors) +
        theme(
            #axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 18)
        )
    return(cov)
}

coverage_values <- c(20, 30, 40, 50)

basenorm_feats_common_calls_54bp_three <- do.call(rbind, lapply(coverage_values, function(covg) read_commoncaller_csv(basepath_basenorm_feats, covg, "A2_A3", "three")))
basenorm_feats_common_calls_125bp_three <- do.call(rbind, lapply(coverage_values, function(covg) read_commoncaller_csv(basepath_basenorm_feats, covg, "AKA-017_GIM-024", "three")))
basenorm_feats_common_calls_154bp_three <- do.call(rbind, lapply(coverage_values, function(covg) read_commoncaller_csv(basepath_basenorm_feats, covg, "JUT-008_MUN-009", "three")))

plot_bp(basenorm_feats_common_calls_54bp_three)
plot_bp(basenorm_feats_common_calls_125bp_three)
plot_bp(basenorm_feats_common_calls_154bp_three)

coverage_values <- c(20, 30, 40, 50)

#basenorm_feats_common_calls_54bp_six <- do.call(rbind, lapply(coverage_values, function(covg) read_commoncaller_csv(basepath_basenorm_feats, covg, "A2_A3", "six")))
#basenorm_feats_common_calls_125bp_six <- do.call(rbind, lapply(coverage_values, function(covg) read_commoncaller_csv(basepath_basenorm_feats, covg, "AKA-017_GIM-024", "six")))
#basenorm_feats_common_calls_154bp_six <- do.call(rbind, lapply(coverage_values, function(covg) read_commoncaller_csv(basepath_basenorm_feats, covg, "JUT-008_MUN-009", "six")))

#plot_bp(basenorm_feats_common_calls_54bp_six)
#plot_bp(basenorm_feats_common_calls_125bp_six)
#plot_bp(basenorm_feats_common_calls_154bp_six)

plot_bp(basenorm_feats_54bp)
plot_bp(basenorm_feats_125bp)
plot_bp(basenorm_feats_150bp)

basenorm_feats_150bp_top3 <- basenorm_feats_150bp %>% filter(caller %in% c("TEMP2", "TEforest", "RetroSeq"))
plot_bp(basenorm_feats_150bp_top3)


plot_genotyping <- function(df) {
    # Reorder the levels of the caller factor to match the order in colors
    df$caller <- factor(df$caller, levels = names(colors))
    df <- df[order(df$caller, decreasing = TRUE), ]
    cov <- ggplot(df, aes(x = covg, y = macro_f1_score, color = caller)) +
        geom_line(size = 1.5, alpha = 0.5) +
        geom_point(size = 3.5, alpha=0.8) +
        # Adding error bars for distance_sd
        #geom_errorbar(
        #  aes(ymin = distance_mean - distance_sd, ymax = distance_mean + distance_sd),
        #  width = 0.2,  # Adjust width as needed for clarity
        #  size = 0.7,   # Adjust line thickness of error bars
        #  alpha = 0.7   # Adjust transparency
        #) +
        #facet_wrap(~Metric) +
        #scale_y_continuous(limits = c(0, 70)) +
        scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50), labels = c(5, 10, 20, 30, 40, 50)) +
        labs(
            title = "",
            x = "Coverage",
            y = "f1_score"
        ) +
        theme_minimal() +
        scale_color_manual(values = colors) +
        theme(
            #axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 18)
        )
    return(cov)
}

plot_genotyping(basenorm_feats_common_calls_54bp_three)
plot_genotyping(basenorm_feats_common_calls_125bp_three)
plot_genotyping(basenorm_feats_common_calls_154bp_three)
#plot_genotyping(basenorm_feats_common_calls_54bp_six)
#plot_genotyping(basenorm_feats_common_calls_125bp_six)
#plot_genotyping(basenorm_feats_common_calls_154bp_six)


performance_150 <- plot_coverage_threeplot(basenorm_feats_150bp_plot)
breakpoints_generic_150 <- plot_bp(basenorm_feats_150bp)
#breakpoints_6commoncallers_150 <- plot_bp(basenorm_feats_common_calls_154bp_six)
breakpoints_3commoncallers_150 <- plot_bp(basenorm_feats_common_calls_154bp_three)
#geno_6commoncallers_150 <- plot_genotyping(basenorm_feats_common_calls_154bp_six)
geno_3commoncallers_150 <- plot_genotyping(basenorm_feats_common_calls_154bp_three)



fn_data <- data.frame(
  Comparison = c("54", "125", "151"),
  FN_candidate = c(72, 67, 60),   # FN in candidate regions stage
  FN_class = c(32, 51, 27),       # FN in classification stage
  TP = c(317, 252, 300)           # TP in classification stage
)

# Transform the data to long format for ggplot
fn_data_long <- fn_data %>%
  pivot_longer(cols = c(FN_candidate, FN_class, TP), 
               names_to = "Category", 
               values_to = "Count") %>%
  group_by(Comparison) %>%
  mutate(Proportion = Count / sum(Count)) %>%  # Calculate the proportions
  mutate(Category = recode(Category, 
                           "FN_candidate" = "False negative (No candidate region identified)", 
                           "FN_class" = "False negative (misclassified)", 
                           "TP" = "True positive"))
fn_data_long$Comparison <- factor(fn_data_long$Comparison, levels = c("54", "125", "151"))

stacked_bar <- ggplot(fn_data_long, aes(x = Comparison, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +  # Thinner bars with black outlines
  geom_text(aes(label = percent(Proportion)), 
            position = position_stack(vjust = 0.5),  # Position labels in the middle of each segment
            size = 3, 
            color = "black") +  # Adjust text size and color as needed
  labs(title = "", 
       x = "Read length", 
       y = "Proportion of calls") +
  scale_fill_manual(values = c("False negative (No candidate region identified)" = "lightblue", 
                               "False negative (misclassified)" = "lightcoral", 
                               "True positive" = "lightgreen")) +
  theme_minimal() +
  theme(legend.title = element_blank())

figure <- ggarrange(
    performance_150$f1_score, performance_150$precision, performance_150$recall, stacked_bar, breakpoints_generic_150, breakpoints_3commoncallers_150, 
    labels = c("A", "B", "C", "D", "E", "F"),
    ncol = 3, nrow = 2, common.legend = T, legend="bottom"
  )


ggsave("/nas/longleaf/home/adaigle/TEforest/plots/performance_150.svg", figure, dpi=300, width = 8.5, height = 8.5)

performance_125 <- plot_coverage_threeplot(basenorm_feats_125bp_plot)
breakpoints_generic_125 <- plot_bp(basenorm_feats_125bp)
#breakpoints_6commoncallers_125 <- plot_bp(basenorm_feats_common_calls_125bp_six)
breakpoints_3commoncallers_125 <- plot_bp(basenorm_feats_common_calls_125bp_three)
#geno_6commoncallers_125 <- plot_genotyping(basenorm_feats_common_calls_125bp_six)
geno_3commoncallers_125 <- plot_genotyping(basenorm_feats_common_calls_125bp_three)


figure_125 <- ggarrange(
    performance_125$f1_score, performance_125$precision, performance_125$recall, stacked_bar, breakpoints_generic_125, breakpoints_3commoncallers_125,
    labels = c("A", "B", "C", "D", "E", "F"),
    ncol = 3, nrow = 2, common.legend = T, legend="bottom"
  )

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/performance_125.svg", figure_125, dpi=300, width = 8.5, height = 8.5)

performance_54 <- plot_coverage_threeplot(basenorm_feats_54bp_plot)
breakpoints_generic_54 <- plot_bp(basenorm_feats_54bp)
#breakpoints_6commoncallers_54 <- plot_bp(basenorm_feats_common_calls_54bp_six)
breakpoints_3commoncallers_54 <- plot_bp(basenorm_feats_common_calls_54bp_three)
#geno_6commoncallers_54 <- plot_genotyping(basenorm_feats_common_calls_54bp_six)
geno_3commoncallers_54 <- plot_genotyping(basenorm_feats_common_calls_54bp_three)


figure_54 <- ggarrange(
    performance_54$f1_score, performance_54$precision, performance_54$recall, stacked_bar, breakpoints_generic_54, breakpoints_3commoncallers_54,
    labels = c("A", "B", "C", "D", "E", "F"),
    ncol = 3, nrow = 2, common.legend = T, legend="bottom"
  )

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/performance_50.svg", figure_54, dpi=300, width = 8.5, height = 8.5)



read_genotype_data <- function(basepath, covg, genome) {
    # Initialize an empty list to store all loaded objects
    loaded_data <- list()
    
    # Iterate over the covg values (c(5, 30))
    for (coverage in covg) {
        # Construct the path based on basepath, coverage, and genome names
        path <- paste0(basepath, coverage, "X/2L_2R_plots/")
        
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


covg <- c(5, 10,20,30, 40, 50)
all_data <- read_genotype_data(basepath_basenorm_feats, covg, "JUT-008_MUN-009")
all_data_125 <- read_genotype_data(basepath_basenorm_feats, covg, "AKA-017_GIM-024")
all_data_54 <- read_genotype_data(basepath_basenorm_feats, covg, "A2_A3")

combine_genotype_data <- function(all_data, covg) {
    combined_genotype_data <- bind_rows(
        lapply(covg, function(c) all_data[[paste0("genotype_data_", c, "X")]])
    )
    return(combined_genotype_data)
}

combined_genotype_data <- combine_genotype_data(all_data, covg)
combined_genotype_data_125 <- combine_genotype_data(all_data_125, covg)
combined_genotype_data_54 <- combine_genotype_data(all_data_54, covg)

plot_genotyping(combined_genotype_data)

combined_genotype_data <- combined_genotype_data %>%
  rowwise() %>%
  mutate(
    tphomo_mean = mean(unlist(tphomo), na.rm = TRUE),
    tphomo_sd = sd(unlist(tphomo), na.rm = TRUE),
    tphomo_se = sd(unlist(tphomo), na.rm = TRUE) / sqrt(sum(!is.na(unlist(tphomo)))),
    tphet_mean = mean(unlist(tphet), na.rm = TRUE),
    tphet_sd = sd(unlist(tphet), na.rm = TRUE),
    tphet_se = sd(unlist(tphet), na.rm = TRUE) / sqrt(sum(!is.na(unlist(tphet))))
  )
combined_genotype_data_125 <- combined_genotype_data_125 %>%
  rowwise() %>%
  mutate(
    tphomo_mean = mean(unlist(tphomo), na.rm = TRUE),
    tphomo_sd = sd(unlist(tphomo), na.rm = TRUE),
    tphomo_se = sd(unlist(tphomo), na.rm = TRUE) / sqrt(sum(!is.na(unlist(tphomo)))),
    tphet_mean = mean(unlist(tphet), na.rm = TRUE),
    tphet_sd = sd(unlist(tphet), na.rm = TRUE),
    tphet_se = sd(unlist(tphet), na.rm = TRUE) / sqrt(sum(!is.na(unlist(tphet))))
  )

combined_genotype_data_54 <- combined_genotype_data_54 %>%
  rowwise() %>%
  mutate(
    tphomo_mean = mean(unlist(tphomo), na.rm = TRUE),
    tphomo_sd = sd(unlist(tphomo), na.rm = TRUE),
    tphomo_se = sd(unlist(tphomo), na.rm = TRUE) / sqrt(sum(!is.na(unlist(tphomo)))),
    tphet_mean = mean(unlist(tphet), na.rm = TRUE),
    tphet_sd = sd(unlist(tphet), na.rm = TRUE),
    tphet_se = sd(unlist(tphet), na.rm = TRUE) / sqrt(sum(!is.na(unlist(tphet))))
  )

plot_combined <- function(df) {
    set.seed(123)  # For reproducibility
    df <- df %>%
        mutate(jittered_covg = jitter(covg, amount = 1))  # Jitter the x-axis (covg) values
    df$caller <- factor(df$caller, levels = names(colors))
    df <- df[order(df$caller, decreasing = TRUE), ]
    ggplot(df) +
        # Plot for Tphomo
        geom_line(aes(x = jittered_covg, y = tphomo_mean, color = caller), size = 1, alpha = 0.8) +
        geom_point(aes(x = jittered_covg, y = tphomo_mean, color = caller), size = 3) +
        #geom_errorbar(aes(x = jittered_covg, ymin = tphomo_mean - tphomo_sd, ymax = tphomo_mean + tphomo_sd, color = caller),
        #              width = 0.5, size = 1, alpha = 0.5) +
        # Plot for Tphet with lighter color
        geom_line(aes(x = jittered_covg, y = tphet_mean, color = caller), size = 1, alpha=0.5) +
        geom_point(aes(x = jittered_covg, y = tphet_mean, color = caller), size = 3, alpha=0.5, shape = 17) +
        #geom_errorbar(aes(x = covg, ymin = tphet_mean - tphet_sd, ymax = tphet_mean + tphet_sd, color = caller),
        #              width = 0.5, size = 1, alpha = 0.5) +
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
basenorm_feats_common_calls_154bp_three
# Generate the combined plot
#combined_genotype_data %>% filter(caller %in% c("TEforest", "RetroSeq", "TEMP2"))
combined_plot <- plot_combined(combined_genotype_data)
combined_plot_125 <- plot_combined(combined_genotype_data_125)
combined_plot_54 <- plot_combined(combined_genotype_data_54)

combined_genotype_data %>% filter(covg==30) %>% select(caller, tphomo_mean, tphet_mean)
combined_genotype_data_125 %>% filter(covg==30) %>% select(caller, tphomo_mean, tphet_mean)
combined_genotype_data_54 %>% filter(covg==30) %>% select(caller, tphomo_mean, tphet_mean)

figure <- ggarrange(
    all_data$TEforest_conf_matrix_40X, performance_150$precision, performance_150$recall, stacked_bar, breakpoints_generic_150, breakpoints_3commoncallers_150, 
    labels = c("A", "B", "C", "D", "E", "F"),
    ncol = 3, nrow = 2, common.legend = T, legend="bottom"
  )


freqplot_df <- all_data$freqplot_data_30X %>% 
    filter(caller %in% c("TEforest_classifier_filter_bps", "retroseq", "temp2")) %>% 
    mutate(caller = recode(caller, 
                           "TEforest_classifier_filter_bps" = "TEforest",
                           "retroseq" = "RetroSeq",
                           "temp2" = "TEMP2"))


freqplot_df_teforest <- freqplot_df %>% 
    filter(caller %in% c("TEforest"))

freqplot_df_125 <- all_data_125$freqplot_data_30X %>% 
    filter(caller %in% c("TEforest_classifier_filter_bps", "retroseq", "temp2")) %>% 
    mutate(caller = recode(caller, 
                           "TEforest_classifier_filter_bps" = "TEforest",
                           "retroseq" = "RetroSeq",
                           "temp2" = "TEMP2"))


freqplot_df_125_teforest <- freqplot_df_125 %>% 
    filter(caller %in% c("TEforest"))

freqplot_df_54 <- all_data_54$freqplot_data_30X %>% 
    filter(caller %in% c("TEforest_classifier_filter_bps", "retroseq", "temp2")) %>% 
    mutate(caller = recode(caller, 
                           "TEforest_classifier_filter_bps" = "TEforest",
                           "retroseq" = "RetroSeq",
                           "temp2" = "TEMP2"))


freqplot_df_54_teforest <- freqplot_df_54 %>% 
    filter(caller %in% c("TEforest"))


frequency_plot <- function(df) {
plt <- ggplot(df, aes(x = truth, y = as.numeric(heterozygosity), color = caller)) +
  geom_jitter(alpha = 0.2, width = 0.1, height = 0) + 
  geom_violin(aes(group = truth, fill = caller), alpha = 0.3) + # Light fill for the violins
  stat_summary(fun = median, geom = "crossbar", aes(group = truth, linetype = "Median"), 
               color = "black", width = 0.05, position = position_dodge(0.9), show.legend = F) + # Crossbar for median with legend
  facet_wrap(~caller) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + # Match fill colors to color values
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("absent", "heterozygous", "homozygous")) + # Custom x-axis labels
  scale_linetype_manual(name = "Statistics", values = c(Median = "solid")) + # Add linetype to legend
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14),
    legend.position = "none",
    axis.text.x = element_text(size = 7.5)
  ) +
  labs(
    title = "",
    x = "True Genotype",
    y = "Predicted prevalence"
  )
return(plt)
}

freqplot_df_teforest <- freqplot_df %>% 
    filter(caller %in% c("TEforest"))

freqplot_df_retroseq <- freqplot_df %>% 
    filter(caller %in% c("RetroSeq"))

freqplot_df_temp2 <- freqplot_df %>% 
    filter(caller %in% c("TEMP2"))

frqplt_teforest <- frequency_plot(freqplot_df_teforest)
frqplt_retroseq <- frequency_plot(freqplot_df_retroseq)
frqplt_temp2 <- frequency_plot(freqplot_df_temp2)

freqplot_df_125_teforest <- freqplot_df_125 %>% 
    filter(caller %in% c("TEforest"))

freqplot_df_125_retroseq <- freqplot_df_125 %>% 
    filter(caller %in% c("RetroSeq"))

freqplot_df_125_temp2 <- freqplot_df_125 %>% 
    filter(caller %in% c("TEMP2"))

frqplt_teforest_125 <- frequency_plot(freqplot_df_125_teforest)
frqplt_retroseq_125 <- frequency_plot(freqplot_df_125_retroseq)
frqplt_temp2_125 <- frequency_plot(freqplot_df_125_temp2)

freqplot_df_54_teforest <- freqplot_df_54 %>% 
    filter(caller %in% c("TEforest"))

freqplot_df_54_retroseq <- freqplot_df_54 %>% 
    filter(caller %in% c("RetroSeq"))

freqplot_df_54_temp2 <- freqplot_df_54 %>% 
    filter(caller %in% c("TEMP2"))

frqplt_teforest_54 <- frequency_plot(freqplot_df_54_teforest)
frqplt_retroseq_54 <- frequency_plot(freqplot_df_54_retroseq)
frqplt_temp2_54 <- frequency_plot(freqplot_df_54_temp2)

alldata_geno <- plot_genotyping(combined_genotype_data)
alldata_geno_125 <-plot_genotyping(combined_genotype_data_125)
alldata_geno_54 <-plot_genotyping(combined_genotype_data_54)
# Extract the legend from combined_plot
combined_plot_legend <- get_legend(combined_plot + theme(legend.position = "bottom"))

confmat <- all_data$TEforest_conf_matrix_30X

confmat +
  labs(x = "Predicted Genotype", y = "True Genotype") # Manually add the axis labels


genotyping_figure <- ggarrange(
  frqplt_teforest, frqplt_retroseq, frqplt_temp2, 
  all_data$TEforest_conf_matrix_30X, all_data$retroseq_conf_matrix_30X, 
  all_data$temp2_conf_matrix_30X, alldata_geno, geno_3commoncallers_150,
  combined_plot,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  ncol = 3, nrow = 3,
  common.legend = TRUE, legend="none", font.label = list(size = 12)
)

# Combine the genotyping figure and the legend manually
final_plot <- plot_grid(
  genotyping_figure,
  combined_plot_legend,
  ncol = 1,
  rel_heights = c(1, 0.1) # Adjust the relative height to control the space for the legend
)

genotyping_figure_125 <- ggarrange(
  frqplt_teforest_125, frqplt_retroseq_125, frqplt_temp2_125, 
  all_data_125$TEforest_conf_matrix_30X, all_data_125$retroseq_conf_matrix_30X, 
  all_data_125$temp2_conf_matrix_30X, alldata_geno_125, geno_3commoncallers_125,
  combined_plot_125,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  ncol = 3, nrow = 3,
  common.legend = TRUE, legend="none", font.label = list(size = 12)
)

# Combine the genotyping figure and the legend manually
final_plot_125 <- plot_grid(
  genotyping_figure_125,
  combined_plot_legend,
  ncol = 1,
  rel_heights = c(1, 0.1) # Adjust the relative height to control the space for the legend
)

genotyping_figure_54 <- ggarrange(
  frqplt_teforest_54, frqplt_retroseq_54, frqplt_temp2_54, 
  all_data_54$TEforest_conf_matrix_30X, all_data_54$retroseq_conf_matrix_30X, 
  all_data_54$temp2_conf_matrix_30X, alldata_geno_54, geno_3commoncallers_54,
  combined_plot_54,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  ncol = 3, nrow = 3,
  common.legend = TRUE, legend="none", font.label = list(size = 12)
)

# Combine the genotyping figure and the legend manually
final_plot_54 <- plot_grid(
  genotyping_figure_54,
  combined_plot_legend,
  ncol = 1,
  rel_heights = c(1, 0.1) # Adjust the relative height to control the space for the legend
)


ggsave("/nas/longleaf/home/adaigle/TEforest/plots/genotyping_150.svg", final_plot, dpi=300, width = 8.5, height = 8.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/genotyping_150.jpg", final_plot, dpi=300, width = 8.5, height = 8.5)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/genotyping_125.svg", final_plot_125, dpi=300, width = 8.5, height = 8.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/genotyping_125.jpg", final_plot_125, dpi=300, width = 8.5, height = 8.5)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/genotyping_54.svg", final_plot_54, dpi=300, width = 8.5, height = 8.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/genotyping_54.jpg", final_plot_54, dpi=300, width = 8.5, height = 8.5)


# investigate temp2 for 125 bp

combined_genotype_data_125_temp2 <- combined_genotype_data_125 %>% filter(caller=="TEMP2")

combined_genotype_data_125_temp2 %>% pull(tphomo_mean)

#load('/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_40X/2L_2R_plots/AKA-017_GIM-024.RData')

#benchmark_mapping_results_genotypes %>% filter(caller=="temp2")

basepath_trimmed_read_experiment <- "/nas/longleaf/home/adaigle/work/test_TEforest/inference_trimreads"


trimmed_reads_125bp <- do.call(rbind, lapply(30, function(covg) read_df_trim(basepath_trimmed_read_experiment, covg, "AKA-017_GIM-024")))
basenorm_feats_125bp_30 <- do.call(rbind, lapply(30, function(covg) read_df(basepath_basenorm_feats, covg, "AKA-017_GIM-024")))
trimmed_reads_125bp$bp <- "100 bp (trimmed)"
basenorm_feats_125bp_30$bp <- "125 bp"
trimmed_reads_125bp <- trimmed_reads_125bp %>%
  add_row(caller = "TEforest", 
          bp = "100 bp (trimmed)", 
          Metric = "f1_score", 
          Score = 0.77, 
          covg = 30) %>%
  add_row(caller = "TEforest", 
          bp = "100 bp (trimmed)", 
          Metric = "recall", 
          Score = 0.65, 
          covg = 30) %>%
  add_row(caller = "TEforest", 
          bp = "100 bp (trimmed)", 
          Metric = "precision", 
          Score = 0.96, 
          covg = 30)
trimplot_125 <- rbind(trimmed_reads_125bp, basenorm_feats_125bp_30)
trimplot_125$bp <- factor(trimplot_125$bp, levels = c("125 bp", "100 bp (trimmed)"))
trimplot_125_F1 <- trimplot_125 %>% filter(Metric=="f1_score")
trimplot_125_recall <- trimplot_125 %>% filter(Metric=="recall")
trimplot_125_precision <- trimplot_125 %>% filter(Metric=="precision")


trimmed_reads_150bp <- do.call(rbind, lapply(30, function(covg) read_df_trim(basepath_trimmed_read_experiment, covg, "JUT-008_MUN-009")))
basenorm_feats_150bp_30 <- do.call(rbind, lapply(30, function(covg) read_df(basepath_basenorm_feats, covg, "JUT-008_MUN-009")))
trimmed_reads_150bp$bp <- "100 bp (trimmed)"
basenorm_feats_150bp_30$bp <- "151 bp"
trimmed_reads_150bp <- trimmed_reads_150bp %>%
  add_row(caller = "TEforest", 
          bp = "100 bp (trimmed)", 
          Metric = "f1_score", 
          Score = 0.84, 
          covg = 30) %>%
  add_row(caller = "TEforest", 
          bp = "100 bp (trimmed)", 
          Metric = "recall", 
          Score = 0.76, 
          covg = 30) %>%
  add_row(caller = "TEforest", 
          bp = "100 bp (trimmed)", 
          Metric = "precision", 
          Score = 0.94, 
          covg = 30)
trimplot_150 <- rbind(trimmed_reads_150bp, basenorm_feats_150bp_30)
trimplot_150$bp <- factor(trimplot_150$bp, levels = c("151 bp", "100 bp (trimmed)"))
trimplot_150_F1 <- trimplot_150 %>% filter(Metric=="f1_score")
trimplot_150_recall <- trimplot_150 %>% filter(Metric=="recall")
trimplot_150_precision <- trimplot_150 %>% filter(Metric=="precision")


plot_f1_125 <- ggplot(trimplot_125_F1, aes(x = caller, y = Score, fill = bp)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "F1 Score", y = "Score", x = "", fill = "Read length") +  
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12), 
        legend.spacing.y = unit(1, "cm"),      
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = bp),
            position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

plot_recall_125 <- ggplot(trimplot_125_recall, aes(x = caller, y = Score, fill = bp)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Recall", y = "Score", x = "", fill = "Read length") +  
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(1, "cm"), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = bp),
            position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

plot_precision_125 <- ggplot(trimplot_125_precision, aes(x = caller, y = Score, fill = bp)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Precision", y = "Score", x = "Caller", fill = "Read length") +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(1, "cm"),     
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = bp),
            position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

# Combine plots with adjusted heights
trimmed_reads_125_plot <- ggarrange(
  plot_f1_125, 
  plot_recall_125, 
  plot_precision_125, 
  ncol = 1, 
  nrow = 3, 
  labels = c("A", "B", "C"),
  common.legend = TRUE, 
  legend = "bottom",
  heights = c(1, 1, 1.4)
)


plot_f1_150 <- ggplot(trimplot_150_F1, aes(x = caller, y = Score, fill = bp)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "F1 Score", y = "Score", x = "", fill = "Read length") +  
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12), 
        legend.spacing.y = unit(1, "cm"),      
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = bp),
            position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

plot_recall_150 <- ggplot(trimplot_150_recall, aes(x = caller, y = Score, fill = bp)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Recall", y = "Score", x = "", fill = "Read length") +  
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(1, "cm"), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = bp),
            position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

plot_precision_150 <- ggplot(trimplot_150_precision, aes(x = caller, y = Score, fill = bp)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Precision", y = "Score", x = "Caller", fill = "Read length") +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(1, "cm"),     
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f", Score), y = Score, group = bp),
            position = position_dodge(width = 0.8), size = 2.75, vjust = -0.5)

# Combine plots with adjusted heights
trimmed_reads_150_plot <- ggarrange(
  plot_f1_150, 
  plot_recall_150, 
  plot_precision_150, 
  ncol = 1, 
  nrow = 3, 
  labels = c("A", "B", "C"),
  common.legend = TRUE, 
  legend = "bottom",
  heights = c(1, 1, 1.4)
)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/trimmed_reads_150_plot.png", trimmed_reads_150_plot, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/trimmed_reads_150_plot.svg", trimmed_reads_150_plot, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/trimmed_reads_150_plot.tiff", trimmed_reads_150_plot, dpi=300, width = 8.5, height = 10.5)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/trimmed_reads_125_plot.png", trimmed_reads_125_plot, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/trimmed_reads_125_plot.svg", trimmed_reads_125_plot, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/trimmed_reads_125_plot.tiff", trimmed_reads_125_plot, dpi=300, width = 8.5, height = 10.5)



make_bp_plot <- function(caller_name, benchmark_mapping_results_summary) {
#if (caller_name == "TEforest_classifier_filter_bps") {
#  caller_name <- "TEforest"
#}

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

num <- as.numeric(row.names(benchmark_mapping_results_summary)[benchmark_mapping_results_summary$caller == caller_name])
# Extract the vector directly

original_distance_data <- benchmark_mapping_results_summary %>% 
    filter(caller==caller_name) %>% filter(covg==50) %>% filter(Metric=="f1_score") %>%
    pull(distance_vectors_df)

original_distance_data_caller <- original_distance_data[[1]]$distance
# Calculate mean and median from the original data
original_mean <- mean(original_distance_data_caller)
original_median <- median(original_distance_data_caller)

# Modify the data: replace values greater than 200 with 200 for the final bin
modified_distance_data <- ifelse(original_distance_data_caller > 25, 25, original_distance_data_caller)

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
  theme_minimal() +
  theme(axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), 
        plot.title=element_text(size=12), 
        legend.title=element_blank(), 
        legend.text=element_text(size=12)) +
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
  ) +
  guides(fill = "none")

  return(plt)
}

temp2_bp_54 <- make_bp_plot("TEMP2", basenorm_feats_54bp)
teforest_bp_54 <- make_bp_plot("TEforest", basenorm_feats_54bp)
retroseq_bp_54 <- make_bp_plot("RetroSeq", basenorm_feats_54bp)
teflon_bp_54 <- make_bp_plot("TEFLoN", basenorm_feats_54bp)

temp2_bp_125 <- make_bp_plot("TEMP2", basenorm_feats_125bp)
teforest_bp_125 <- make_bp_plot("TEforest", basenorm_feats_125bp)
retroseq_bp_125 <- make_bp_plot("RetroSeq", basenorm_feats_125bp)
teflon_bp_125 <- make_bp_plot("TEFLoN", basenorm_feats_125bp)

temp2_bp_150 <- make_bp_plot("TEMP2", basenorm_feats_150bp)
teforest_bp_150 <- make_bp_plot("TEforest", basenorm_feats_150bp)
retroseq_bp_150 <- make_bp_plot("RetroSeq", basenorm_feats_150bp)
teflon_bp_150 <- make_bp_plot("TEFLoN", basenorm_feats_150bp)

bp_distribution_figure <- ggarrange(
    teforest_bp_150, temp2_bp_150, retroseq_bp_150, teflon_bp_150,
    labels = c("A", "B", "C", "D"),
    ncol = 1, nrow = 4, common.legend = T, legend="right"
  )

bp_distribution_figure_125 <- ggarrange(
    teforest_bp_125, temp2_bp_125, retroseq_bp_125, teflon_bp_125,
    labels = c("A", "B", "C", "D"),
    ncol = 1, nrow = 4, common.legend = T, legend="right"
  )

bp_distribution_figure_54 <- ggarrange(
    teforest_bp_54, temp2_bp_54, retroseq_bp_54, teflon_bp_54,
    labels = c("A", "B", "C", "D"),
    ncol = 1, nrow = 4, common.legend = T, legend="right"
  )

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_150_plot.png", bp_distribution_figure, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_150_plot.svg", bp_distribution_figure, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_150_plot.tiff", bp_distribution_figure, dpi=300, width = 8.5, height = 10.5)


ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_125_plot.png", bp_distribution_figure_125, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_125_plot.svg", bp_distribution_figure_125, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_125_plot.tiff", bp_distribution_figure_125, dpi=300, width = 8.5, height = 10.5)

ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_54_plot.png", bp_distribution_figure_54, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_54_plot.svg", bp_distribution_figure_54, dpi=300, width = 8.5, height = 10.5)
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/bp_distribution_figure_54_plot.tiff", bp_distribution_figure_54, dpi=300, width = 8.5, height = 10.5)

