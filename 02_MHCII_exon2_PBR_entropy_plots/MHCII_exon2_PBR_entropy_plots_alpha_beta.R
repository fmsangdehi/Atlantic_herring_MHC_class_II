library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggseqlogo)

#------------------------------
# Generate plot for alpha chain
#------------------------------

# Define the alignment file path
alignment_file_dir <- "path/to/aln/files/dir"
alignment_file <- "AA_aln_DAA.fa"

# Read the alignment data and trim it
alignment <- readAAStringSet(paste0(alignment_file_dir, alignment_file))
alignment <- subseq(alignment, start = 20, end = 107)
alignment_matrix <- as.matrix(alignment)

# Convert alignment to a data frame
alignment_df <- as.data.frame(alignment_matrix)
alignment_df$sequence <- rownames(alignment_df)
alignment_df$index <- seq_len(nrow(alignment_df))

# Convert to long format
alignment_long <- alignment_df %>%
  pivot_longer(-c(sequence, index), names_to = "position", values_to = "residue") %>%
  mutate(position = as.numeric(sub("V", "", position)))


# Annotation Panel
annotation_plot <- ggplot() +
  annotate("point", x = c(c(5, 7, 9, 22, 24, 26, 31, 32, 43, 51, 53, 54, 55, 61, 65, 68, 71, 72, 75, 79)), y = 2, shape = 15, size = 2.5, color = "grey10") +
  annotate("point", x = c(5, 7, 9, 22, 24, 31, 32, 43, 51, 52, 53, 54, 55, 61, 64, 65, 68, 69, 71, 72, 75, 76, 79), y = 1, shape = 15, size = 2.5, color = "red3") +
  
  scale_y_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("", "")) +
  scale_x_continuous(limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0)) +
  
  theme_void() +
  theme(axis.text.y = element_text(size = 10, hjust = 0),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"))


# Function to calculate Shannon entropy for each position
calculate_entropy <- function(column) {
  column <- column[column != "-"]
  if (length(column) == 0) return(0)
  
  # Calculate frequency of each amino acid in the column
  aa_freq <- table(column) / length(column)
  
  # Compute Shannon entropy
  -sum(aa_freq * log2(aa_freq))
}

# Calculate entropy for each column (position) in the alignment matrix
Herring_consensus_seq_count <- 195
entropy_values <- apply(alignment_matrix[1:Herring_consensus_seq_count,], 2, calculate_entropy)
normalized_entropy <- entropy_values/log2(20)

# Convert to a data frame
entropy_df <- data.frame(Position = seq_along(normalized_entropy), ShannonEntropy = normalized_entropy)

asterisk_data <- data.frame(
  x_pos = c(14, 15, 16, 17),
  asterisk_label = "*"
)

# Create entropy plot
entropy_plot <- ggplot(entropy_df, aes(x = Position, y = ShannonEntropy, fill = ShannonEntropy)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.4, limits = c(0, 0.8)) +
  labs(y = "Normalized\nentropy") +
  scale_x_continuous(breaks = seq(1, max(alignment_long$position), by = 10), labels = seq(20, max(alignment_long$position) + 19, by = 10), limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0), position = "bottom") +
  scale_y_continuous(breaks = c(0, 0.4, 0.8), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.8), expand = c(0, 0)) +
  geom_text(data = asterisk_data, mapping = aes(x = x_pos, label = asterisk_label), y = -0.09, size = 10, fontface = "bold", color = alpha("#FFEF00", 1), inherit.aes = FALSE) +
  coord_cartesian(clip = "off") + # Critical for visibility of yellow segment
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.ticks.y = element_blank(),
    axis.ticks.x.bottom = element_line(color = "black")
  )


# Function to generate a smooth coil
generate_shifted_gradient_wave_coil <- function(start_x = 0, end_x = 36, y_center = 1, coils = 10, amplitude = 0.3, resolution = 5000) {
  x_positions <- seq(start_x, end_x, length.out = resolution)
  
  y_positions <- y_center + amplitude * sin(seq(0, coils * 2 * pi, length.out = resolution))
  
  gradient <- sin(seq(0, coils * 2 * pi, length.out = resolution) - pi / 2)
  
  data.frame(x = x_positions, y = y_positions, gradient = gradient)
}

# Generate data for the shifted gradient coil
coil_data <- generate_shifted_gradient_wave_coil(
  start_x = 0, end_x = 36, y_center = 1, 
  coils = 10, amplitude = 0.3, resolution = 5000
)

# Generate secondary structure elements plot
alpha_helix <- data.frame(
  start = c(58),
  end = c(80),
  y = c(1),
  label = c("A1")) %>%
  mutate(length = end - start)

beta_strand <- data.frame(
  start = c(2, 21, 29, 39),
  end = c(11, 26, 35, 43),
  y = c(1, 1, 1, 1),
  label = c("B1", "B2", "B3", "B4")
)

arrow_length <- 0.3

secondary_structure_plot <- ggplot() +
  geom_segment(data = beta_strand, aes(x = start - 0.5, xend = end - 2*arrow_length + 0.5, y = y, yend = y), arrow = arrow(length = unit(arrow_length, "cm"), type = "closed", angle = 25), color = "#ECCD00", linewidth = 3, linejoin = "mitre") +
  
  geom_path(data = coil_data[coil_data$x <= alpha_helix$length[1] + 1,], aes(x = x + alpha_helix$start[1] - 0.5, y = y, color = gradient), linewidth = 1.7, linejoin = "round", lineend = "round") +
  scale_color_gradient(low = "#87d0f5", high = "#085997") +

  geom_text(data = alpha_helix, aes(x = (start + end) / 2, y = y - 0.7, label = label), size = 7, color = "black") +
  geom_text(data = beta_strand, aes(x = (start + end) / 2, y = y - 0.7, label = label), size = 7, color = "black") +
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


# Combine plots
combined_plot_DAA <- ggdraw() +
  draw_label("a", x = 0, y = 1, hjust = 0, vjust = 0.8, fontface = "bold", size = 32) +
  draw_plot(plot_grid(annotation_plot, entropy_plot, secondary_structure_plot, ncol = 1, align = "v", rel_heights = c(0.4, 2.5, 1)),
            x = 0.015, y = 0, width = 0.985, height = 0.9) +
  theme(plot.margin = margin(20, 10, 10, 10))


#------------------------------
# Generate plot for beta chain
#------------------------------

# Define the alignment file path
alignment_file_dir <- "path/to/aln/files/dir"
alignment_file <- "AA_aln_DAB.fa"

# Read the alignment data and trim it
alignment <- readAAStringSet(paste0(alignment_file_dir, alignment_file))
alignment <- subseq(alignment, start = 21, end = 117)
alignment_matrix <- as.matrix(alignment)

# Convert alignment to a data frame for ggplot2
alignment_df <- as.data.frame(alignment_matrix)
alignment_df$sequence <- rownames(alignment_df)
alignment_df$index <- seq_len(nrow(alignment_df))

# Convert to long format
alignment_long <- alignment_df %>%
  pivot_longer(-c(sequence, index), names_to = "position", values_to = "residue") %>%
  mutate(position = as.numeric(sub("V", "", position)))


# Annotation Panel
annotation_plot <- ggplot() +
  annotate("point", x = c(8, 10, 12, 30, 32, 34, 39, 40, 49, 58, 62, 63, 67, 70, 72, 73, 76, 80, 83, 84, 87, 88, 90, 91), y = 2, shape = 15, size = 2.5, color = "black") +
  annotate("point", x = c(8, 9, 10, 12, 28, 30, 32, 33, 34, 39, 49, 59, 62, 63, 69, 72, 73, 76, 79, 80, 83, 84, 87, 88), y = 1, shape = 15, size = 2.5, color = "red3") +
  
  scale_y_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("", "")) +
  scale_x_continuous(limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0)) +
  
  theme_void() +
  theme(axis.text.y = element_text(size = 10, hjust = 0),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"))


# Function to calculate Shannon entropy for each position
calculate_entropy <- function(column) {
  column <- column[column != "-"]
  if (length(column) == 0) return(0)
  
  # Calculate frequency of each amino acid in the column
  aa_freq <- table(column) / length(column)
  
  # Compute Shannon entropy
  -sum(aa_freq * log2(aa_freq))
}

# Calculate entropy for each column (position) in the alignment matrix
Herring_consensus_seq_count <- 220
entropy_values <- apply(alignment_matrix[1:Herring_consensus_seq_count,], 2, calculate_entropy)
normalized_entropy <- entropy_values/log2(20)

# Convert to a data frame
entropy_df <- data.frame(Position = seq_along(normalized_entropy), ShannonEntropy = normalized_entropy)

# Create entropy plot
entropy_plot <- ggplot(entropy_df, aes(x = Position, y = ShannonEntropy, fill = ShannonEntropy)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.4, limits = c(0, 0.8)) +
  labs(y = "Normalized\nentropy") +
  scale_x_continuous(breaks = c(1, seq(10, max(alignment_long$position), by = 10)), labels = c(21, seq(30, max(alignment_long$position) + 20, by = 10)), limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0), position = "bottom") +
  scale_y_continuous(breaks = c(0, 0.4, 0.8), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.8), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.ticks.y = element_blank(),
    axis.ticks.x.bottom = element_line(color = "black")
  )


# Function to generate a smooth coil
generate_shifted_gradient_wave_coil <- function(start_x = 0, end_x = 36, y_center = 1, coils = 10, amplitude = 0.3, resolution = 5000) {
  x_positions <- seq(start_x, end_x, length.out = resolution)
  
  y_positions <- y_center + amplitude * sin(seq(0, coils * 2 * pi, length.out = resolution))
  
  gradient <- sin(seq(0, coils * 2 * pi, length.out = resolution) - pi / 2)
  
  data.frame(x = x_positions, y = y_positions, gradient = gradient)
}

# Generate data for the shifted gradient coil
coil_data <- generate_shifted_gradient_wave_coil(
  start_x = 0, end_x = 36, y_center = 1, 
  coils = 10, amplitude = 0.3, resolution = 5000
)

# Generate secondary structure elements plot
alpha_helix <- data.frame(
  start = c(54, 67, 83),
  end = c(65, 80, 89),
  y = c(1, 1, 1),
  label = c("A1", "A2", "A3")) %>%
  mutate(length = end - start)

beta_strand <- data.frame(
  start = c(7, 25, 37, 47),
  end = c(18, 34, 43, 51),
  y = c(1, 1, 1, 1),
  label = c("B1", "B2", "B3", "B4")
)

arrow_length <- 0.3

secondary_structure_plot <- ggplot() +
  geom_segment(data = beta_strand, aes(x = start - 0.5, xend = end - 2*arrow_length + 0.5, y = y, yend = y), arrow = arrow(length = unit(arrow_length, "cm"), type = "closed", angle = 25), color = "#ECCD00", linewidth = 3, linejoin = "mitre") +
  
  geom_path(data = coil_data[coil_data$x <= alpha_helix$length[1] + 1,], aes(x = x + alpha_helix$start[1] - 0.5, y = y, color = gradient), linewidth = 1.7, linejoin = "round", lineend = "round") +
  geom_path(data = coil_data[coil_data$x <= alpha_helix$length[2] + 1,], aes(x = x + alpha_helix$start[2] - 0.5, y = y, color = gradient), linewidth = 1.7, linejoin = "round", lineend = "round") +
  geom_path(data = coil_data[coil_data$x <= alpha_helix$length[3] + 1,], aes(x = x + alpha_helix$start[3] - 0.5, y = y, color = gradient), linewidth = 1.7, linejoin = "round", lineend = "round") +
  scale_color_gradient(low = "#87d0f5", high = "#085997") +

  geom_text(data = alpha_helix, aes(x = (start + end) / 2, y = y - 0.7, label = label), size = 7, color = "black") +
  geom_text(data = beta_strand, aes(x = (start + end) / 2, y = y - 0.7, label = label), size = 7, color = "black") +
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, max(alignment_long$position) + 1), expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


# Combine plots
combined_plot_DAB <- ggdraw() +
  draw_label("b", x = 0, y = 1, hjust = 0, vjust = 0.8, fontface = "bold", size = 32) +
  draw_plot(plot_grid(annotation_plot, entropy_plot, secondary_structure_plot, ncol = 1, align = "v", rel_heights = c(0.4, 2.5, 1)),
            x = 0.015, y = 0, width = 0.985, height = 0.9) +
  theme(plot.margin = margin(20, 10, 10, 10))


#-----------------------------
# Combine alpha and beta plots
#-----------------------------

combined_plot <- plot_grid(combined_plot_DAA, combined_plot_DAB, ncol = 1, align = "v", rel_heights = c(1, 1))

ggsave("PBR_plots_combined.pdf", plot = combined_plot, device = "pdf", width = 12, height = 7.5)
