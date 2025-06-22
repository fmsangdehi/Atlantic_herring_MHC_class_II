library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggseqlogo)

# Define the alignment file path
alignment_file <- "MHCII_exon2_herring_and_comparative_set_aln_curated_alpha.fa"

# Read the alignment data
alignment <- readAAStringSet(alignment_file)
alignment_matrix <- as.matrix(alignment)

# Convert alignment to a data frame
alignment_df <- as.data.frame(alignment_matrix)
alignment_df$sequence <- rownames(alignment_df)

# Generate y axis coordinates with 0.5 gap between clades
clade_sizes <- c(13, 6, 17, 5, 4, 5, 3, 9)
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
  annotate("rect", xmin = 12 - 0.5, xmax = 12 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "green2", alpha = 0.6) +
  annotate("rect", xmin = 91 - 0.5, xmax = 91 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "green2", alpha = 0.6) +
  annotate("rect", xmin = 79 - 0.5, xmax = 79 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "tan1", alpha = 0.6) +
  annotate("rect", xmin = 94 - 0.5, xmax = 94 + 0.5, ymin = 0.5, ymax = max(y_coords) + 0.5, fill = "tan1", alpha = 0.6) +
  annotate("rect", xmin = 15 - 0.5, xmax = 18 + 0.5, ymin = 5 - 0.5, ymax = 5 + 0.5, fill = "#FFEF00", alpha = 1) +
  
  geom_hline(yintercept = c(13.75, 20.25, 37.75, 43.25, 47.75, 53.25, 56.75), linetype = "solid", color = "black") +
  
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
  
  scale_x_continuous(breaks = seq(1, max(alignment_long$position), by = 10), labels = seq(20, max(alignment_long$position) + 19, by = 10), limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0), position = "top") +

  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(size = 12, vjust = 1),
        axis.ticks.x.bottom = element_line(color = "black"),
        axis.text.x.top = element_text(size = 12, vjust = 0),
        axis.ticks.x.top = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 5, 0, 0)
        )


color_strip_df <- data.frame(
  y = rev(y_coords),
  labels = rev(alignment_df$sequence),
  colors = rev(c(rep("#FFB199", 13), rep("#80CCFF", 6), rep("#B6EDB0", 17), rep("#CCCCCC", 5), rep("white", 4), rep("#FFE68C", 5), rep("#C7B3FF", 3), rep("#FFA6C9", 9))))

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
  annotate("point", x = c(5, 8, 10, 25, 27, 29, 34, 35, 53, 61, 63, 64, 65, 75, 79, 90, 93, 94, 97, 101), y = 2, shape = 15, size = 2.5, color = "grey10") +
  annotate("point", x = c(5, 8, 10, 25, 27, 34, 35, 53, 61, 62, 63, 64, 65, 75, 78, 79, 90, 91, 93, 94, 97, 98, 101), y = 1, shape = 15, size = 2.5, color = "red3") +
  
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
  rel_heights = c(0.978, 0.022)
)

ggsave("alpha_alignment_plot.pdf", plot = combined_plot, device = "pdf", width = 15, height = 15)
