# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggh4x)

path <- "imaging/profiles/"
fixed.colors <- c(
  "Mesp2" = "#92F2AB",
  "Tbx18" = "#C200C2",
  "Uncx" = "#9ABDE4",
  "Hoxd4" = "#C200C2",
  "Hoxb9" = "#ECE819",
  "FOXC1" = "#C200C2",
  "SOX2" = "#9ABDE4"
)
input.file <- "HCR.Mesp.Uncx.Tbx18.KO.csv"
# input.file <- "HCR.Hox.KO.csv"
# input.file <- "IF.FOXC1.SOX2.Mosaic.csv"
# Load the data
data <- read.csv(paste0(path, input.file))
# The data for each pixel need to be binned
# Set the number of bins to split each gastruloid to the same number of bins
nbins <- 200
# # In another analysis, bins are of fixed size:
# abs.bin.length <- 10
# Define plotting directory
path_for_plots <- file.path("output.files", path)

my_theme <-
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank()
    )

# Create the output directory
if (!dir.exists(path_for_plots)) {
  dir.create(path_for_plots, recursive = TRUE)
}

# Breaks as proportion are computed
breaks <- seq(from = 0, to = 1, length.out = nbins)
# As well as the coordinate of the middle of each bin
mid.bins <- breaks[-nbins] + 1 / (2 * (nbins - 1))
# # Breaks abs are computed
# breaks.abs <- seq(from = 0, to = max(data$Distance + abs.bin.length), by = abs.bin.length)
# # As well as the coordinate of the middle of each bin
# mid.bins.abs <- breaks.abs[-length(breaks.abs)] + 1 / 2 * (abs.bin.length)
# Normalize Distance and Fluorescence for each Marker, Replicate, and Condition between 0 and 1
data_long <- data %>%
  mutate(Replicate = as.factor(Replicate)) %>% # Ensure Replicate is a factor
  group_by(Marker, Condition, Replicate, Time) %>%
  mutate(
    # Normalize Distance within each Marker/Condition/Replicate group
    Normalized_Distance = (Distance - min(Distance, na.rm = TRUE)) /
      (max(Distance, na.rm = TRUE) - min(Distance, na.rm = TRUE)),
    Bin = cut(Normalized_Distance, breaks = breaks, include.lowest = TRUE, labels = mid.bins),
    # Bin.abs = cut(Distance, breaks = breaks.abs, include.lowest = TRUE, labels = mid.bins.abs)
  ) %>%
  ungroup()

# Calculate Average of the normalized fluorescence for each normalized distance, Condition, and Marker
average_data_replicate <- data_long %>%
  group_by(Marker, Bin, Condition, Replicate, Time) %>%
  summarise(
    Average_single_Fluorescence = mean(Fluorescence, na.rm = TRUE)
  ) %>%
  group_by(Marker, Condition, Replicate, Time) %>%
    mutate(
    # Normalize Fluorescence within each Marker/Condition/Replicate group
    Normalized_Fluorescence = (Average_single_Fluorescence - min(Average_single_Fluorescence, na.rm = TRUE)) /
      (max(Average_single_Fluorescence, na.rm = TRUE) - min(Average_single_Fluorescence, na.rm = TRUE))
  )
# Average/Compute SD between replicates and time
average_data <- average_data_replicate %>%
  group_by(Marker, Bin, Condition) %>%
  summarise(
    Average_NormFluorescence = mean(Normalized_Fluorescence, na.rm = TRUE),
    SD_NormFluorescence = sd(Normalized_Fluorescence, na.rm = TRUE)
  )
# # Same for bins of absolute length
# average_data_replicate.abs <- data_long %>%
#   group_by(Marker, Bin.abs, Condition, Replicate, Time) %>%
#   summarise(
#     Average_single_Fluorescence = mean(Fluorescence, na.rm = TRUE)
#   ) %>%
#   group_by(Marker, Condition, Replicate, Time) %>%
#   mutate(
#     # Normalize Fluorescence within each Marker/Condition/Replicate group
#     Normalized_Fluorescence = (Average_single_Fluorescence - min(Average_single_Fluorescence, na.rm = TRUE)) /
#       (max(Average_single_Fluorescence, na.rm = TRUE) - min(Average_single_Fluorescence, na.rm = TRUE))
#   )
# Define the span of the loess smooth
my.span <- 0.05

# Individual plots
g <- ggplot(data = average_data_replicate, aes(x = as.numeric(as.character(Bin)))) +
  # Plot the smoothed average lines
  geom_smooth(aes(y = Normalized_Fluorescence, color = Marker),
    method = "loess", span = my.span, linewidth = 0.5, se = FALSE
  ) +
  facet_nested(Time + Condition ~ Replicate, axes = "all") +
  # Add labels and themes
  labs(x = "A-P axis", y = "Fluorescence (A.U.)") +
  my_theme +
  # Set custom colors
  scale_color_manual(values = fixed.colors) +
  xlim(0.1, 1) # Set x-axis range

nb_rows <- nrow(unique(average_data_replicate[, c("Time", "Condition")]))
nb_cols <- length(unique(average_data_replicate$Replicate))

ggsave(
  file.path(path_for_plots, paste0("individual", gsub(".csv", "", input.file), ".pdf")),
  width = 2 + 2 * nb_cols, height = 2 + 2 * nb_rows
)

# Summary plots
g <- ggplot(data = average_data, aes(x = as.numeric(as.character(Bin)))) +
  # Plot the smoothed average lines
  geom_line(aes(y = Average_NormFluorescence, color = Marker),
    linewidth = 0.5
  ) +
  # Add ribbon
  geom_ribbon(
    aes(
      ymin = Average_NormFluorescence - SD_NormFluorescence,
      ymax = Average_NormFluorescence + SD_NormFluorescence,
      fill = Marker
    ),
    alpha = 0.2
  ) +
  facet_nested(. ~ Condition, axes = "all") +
  # Add labels and themes
  labs(x = "A-P axis", y = "Fluorescence (A.U.)") +
  my_theme +
  # Set custom colors
  scale_color_manual(values = fixed.colors) +
  scale_fill_manual(values = fixed.colors) +
  xlim(0.1, 1) # Set x-axis range

ggsave(
  file.path(path_for_plots, paste0("average_not_smoothed_time_pooled_", gsub(".csv", "", input.file), ".pdf")),
  width = 14, height = 2.5
)

if (length(unique(average_data_replicate$Time)) > 1) {
  # Average/Compute SD between replicates and Time
  average_data <- average_data_replicate %>%
    group_by(Marker, Bin, Condition, Time) %>%
    summarise(
      Average_NormFluorescence = mean(Normalized_Fluorescence, na.rm = TRUE),
      SD_NormFluorescence = sd(Normalized_Fluorescence, na.rm = TRUE)
    )
  # Summary plots
  g <- ggplot(data = average_data, aes(x = as.numeric(as.character(Bin)))) +
    # Plot the smoothed average lines
    geom_line(aes(y = Average_NormFluorescence, color = Marker),
      linewidth = 0.5
    ) +
    # Add ribbon
    geom_ribbon(
      aes(
        ymin = Average_NormFluorescence - SD_NormFluorescence,
        ymax = Average_NormFluorescence + SD_NormFluorescence,
        fill = Marker
      ),
      alpha = 0.2
    ) +
    facet_nested(. ~ Time + Condition, axes = "all") +
    # Add labels and themes
    labs(x = "A-P axis", y = "Fluorescence (A.U.)") +
    my_theme +
    # Set custom colors
    scale_color_manual(values = fixed.colors) +
    scale_fill_manual(values = fixed.colors) +
    xlim(0.1, 1) # Set x-axis range

  ggsave(
    file.path(path_for_plots, paste0("average_not_smoothed_", gsub(".csv", "", input.file), ".pdf")),
    width = 14, height = 2.5
  )
}
if ("Uncx" %in% average_data_replicate$Marker) {
  # One figure main: conditions as colums 1 row
  # KO: 1
  # WT: 1
  # TLS: 1
  g <- ggplot(
    data = subset(average_data_replicate, Replicate == 1 & Condition != "Gastruloid-Mosaic"),
    aes(x = as.numeric(as.character(Bin)))
  ) +
    # Plot the smoothed average lines
    geom_smooth(aes(y = Normalized_Fluorescence, color = Marker),
      method = "loess", span = my.span, linewidth = 1, se = FALSE
    ) +
    facet_nested(. ~ Time + Condition, axes = "all") +
    # Add labels and themes
    labs(x = "A-P axis", y = "Fluorescence (A.U.)") +
    my_theme +
    # Set custom colors
    scale_color_manual(values = fixed.colors) +
    xlim(0, 1) # Set x-axis range
  ggsave(
    file.path(path_for_plots, paste0("first_rep_", gsub(".csv", "", input.file), ".pdf")),
    width = 10, height = 2
  )
  # One figure sup (S10): replicate as row and conditions as columns
  # Mosaic: 1, 2, 6
  # KO et WT: replicate 4:6
  selected <- data.frame(
    Condition = rep(paste0("Gastruloid-", c("WT", "Mosaic", "KO")), each = 3),
    Replicate = c(4:6, 1, 2, 6, 4:6)
  )
  average_data_replicate_selected <- merge(average_data_replicate, selected)
  average_data_replicate_selected$Condition <- factor(
    average_data_replicate_selected$Condition,
    levels = unique(selected$Condition)
  )
  # Rename replicates for visualization purposes
  average_data_replicate_selected <- average_data_replicate_selected %>%
    mutate(fakeReplicate = case_when(
      Replicate == 1 ~ "4",
      Replicate == 2 ~ "5",
      .default = Replicate
    ))
  g <- ggplot(data = average_data_replicate_selected, aes(x = as.numeric(as.character(Bin)))) +
    # Plot the smoothed average lines
    geom_smooth(aes(y = Normalized_Fluorescence, color = Marker),
      method = "loess", span = my.span, linewidth = 0.5, se = FALSE
    ) +
    facet_nested(fakeReplicate ~ Time + Condition, axes = "all") +
    # Add labels and themes
    labs(x = "A-P axis", y = "Fluorescence (A.U.)") +
    my_theme +
    # Set custom colors
    scale_color_manual(values = fixed.colors) +
    xlim(0.1, 1) # Set x-axis range

  ggsave(
    file.path(path_for_plots, paste0("selected1_", gsub(".csv", "", input.file), ".pdf")),
    width = 37 / 4, height = 15 / 4
  )
  # Second figure sup (S9): nested genotype - replicate 1 row
  # KO et WT: replicate 2:3
  # TLS: replicate 2:3
  g <- ggplot(
    data = subset(average_data_replicate, Replicate %in% 2:3 & Condition != "Gastruloid-Mosaic"),
    aes(x = as.numeric(as.character(Bin)))
  ) +
    # Plot the smoothed average lines
    geom_smooth(aes(y = Normalized_Fluorescence, color = Marker),
      method = "loess", span = my.span, linewidth = 0.7, se = FALSE
    ) +
    facet_nested(. ~ Time + Condition + Replicate, axes = "all", nest_line = element_line()) +
    # Add labels and themes
    labs(x = "A-P axis", y = "Fluorescence (A.U.)") +
    my_theme +
    # Set custom colors
    scale_color_manual(values = fixed.colors) +
    xlim(0.1, 1) # Set x-axis range
  ggsave(
    file.path(path_for_plots, paste0("second_third_rep_", gsub(".csv", "", input.file), ".pdf")),
    width = 16, height = 2.5
  )
}
