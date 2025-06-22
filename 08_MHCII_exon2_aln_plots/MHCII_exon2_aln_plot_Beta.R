library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggseqlogo)

# Define the alignment file path
alignment_file <- "MHCII_exon2_herring_and_comparative_set_aln_curated_beta.fa"

# Read the alignment data
alignment <- readAAStringSet(alignment_file)
alignment_matrix <- as.matrix(alignment)

# Convert alignment to a data frame
alignment_df <- as.data.frame(alignment_matrix)
alignment_df$sequence <- rownames(alignment_df)

# Generate y axis coordinates with 0.5 gap between clades
clade_sizes <- c(14, 8, 17, 5, 4, 5, 3, 9)
y_coords <- c()
offset <- 0

for (size in clade_sizes) {
  y_coords <- c(y_coords, seq_len(size) + offset)
  offset <- max(y_coords) + 0.5  # move starting point for next clade
}

alignment_df$index <- y_coords


# Convert to long format
alignment_long <- alignment_df %>%
  pivot_longer(-c(sequence, index), names_to = "position", values_to = "residue") %>%
  mutate(position = as.numeric(sub("V", "", position)))


# Alignment Panel
alignment_plot <- ggplot(alignment_long, aes(x = position, y = index)) +
  annotate("rect", xmin = 14 - 0.5, xmax = 14 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "green2", alpha = 0.6) +
  annotate("rect", xmin = 88 - 0.5, xmax = 88 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "green2", alpha = 0.6) +
  annotate("rect", xmin = 90 - 0.5, xmax = 90 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "tan1", alpha = 0.6) +
  annotate("rect", xmin = 91 - 0.5, xmax = 91 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "tan1", alpha = 0.6) +

  geom_hline(yintercept = c(14.75, 23.25, 40.75, 46.25, 50.75, 56.25, 59.75), linetype = "solid", color = "black") +
  
  geom_text(aes(
    label = residue,
    color = case_when(
      residue %in% c("G", "S", "T", "Y", "C") ~ "#009900",                      # Polar
      residue %in% c("D", "E") ~ "#CC0000",                                   # Acidic
      residue %in% c("K", "R", "H") ~ "#0073E6",                             # Basic
      residue %in% c("Q", "N") ~ "purple",                           # Neutral
      residue %in% c("A", "V", "L", "I", "P", "W", "F", "M") ~ "black",     # Hydrophobic
      TRUE ~ "black"  # fallback for unknowns (e.g., gaps or X)
    )
  ), size = 3.5) + #fontface = "bold"
  scale_color_identity() +
  
  scale_y_continuous(breaks = alignment_df$index, labels = alignment_df$sequence, trans = "reverse", expand = c(0, 0)) +
  
  scale_x_continuous(breaks = c(1, seq(10, max(alignment_long$position), by = 10)), labels = c(21, seq(30, max(alignment_long$position) + 20, by = 10)), limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0), position = "top") +

  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(size = 12, vjust = 1),
        axis.ticks.x.bottom = element_line(color = "black"),
        axis.text.x.top = element_text(size = 12, vjust = 0),
        axis.ticks.x.top = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 5, 0, 0) # Reduced margins for alignment plot
        )


color_strip_df <- data.frame(
  y = rev(y_coords),
  labels = rev(alignment_df$sequence),
  colors = rev(c(rep("#FFB199", 14), rep("#80CCFF", 8), rep("#B6EDB0", 17), rep("#CCCCCC", 5), rep("white", 4), rep("#FFE68C", 5), rep("#C7B3FF", 3), rep("#FFA6C9", 9))))

color_strip_plot <- ggplot(color_strip_df, aes(x = 1, y = y, fill = colors)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_identity() +  # tell ggplot to use actual color codes
  geom_text(aes(label = labels), color = "black", size = 4, hjust = 0, nudge_x = -0.5) +
  scale_y_continuous(breaks = y_coords, expand = c(0, 0), trans = "reverse") +
  coord_cartesian(xlim = c(0.5, 1.5)) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 5)
  )


# Annotation Panel
annotation_plot <- ggplot() +
  annotate("point", x = c(8, 10, 12, 28, 30, 32, 37, 38, 47, 65, 69, 70, 74, 77, 79, 80, 83, 87, 90, 91, 95, 96, 98, 99), y = 2, shape = 15, size = 2.5, color = "black") +
  annotate("point", x = c(8, 9, 10, 12, 26, 28, 30, 31, 32, 37, 47, 66, 69, 70, 76, 79, 80, 83, 86, 87, 90, 91, 95, 96), y = 1, shape = 15, size = 2.5, color = "red3") +
  
  scale_y_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("", "")) +
  scale_x_continuous(limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0)) +
  
  theme_void() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 5, 5, 0))  # Reduced margins for annotation plot


# Create a blank plot
blank_plot <- ggplot() + theme_void()


# Combine all panels
combined_plot <- plot_grid(
  plot_grid(color_strip_plot, alignment_plot, ncol = 2, rel_widths = c(1, 6.3), align = "h"),
  plot_grid(blank_plot, annotation_plot, ncol = 2, rel_widths = c(1, 6.3), align = "h"),
  ncol = 1,
  rel_heights = c(0.98, 0.02)
)

ggsave("beta_alignment_plot.pdf", plot = combined_plot, device = "pdf", width = 15, height = 17)
