library(ape)
library(dendextend)

# Load the Newick trees
alpha_tree <- read.tree("L4_DAA4_cluster.nwk")
beta_tree <- read.tree("L4_DAB4_cluster.nwk")

alpha_tree$tip.label <- gsub("DAA", "DAB", alpha_tree$tip.label)
    
# Convert trees to dendrogram format
alpha_dend <- as.dendrogram(as.hclust(alpha_tree))
beta_dend <- as.dendrogram(as.hclust(beta_tree))

# Extract labels
alpha_labels <- labels(alpha_dend)
beta_labels <- labels(beta_dend)

# Define colors in alpha tree
alpha_col <- c(rep("#CE16EB", 1), rep("#ff6900", 19), rep("#4169e1", 11))

# Define a Mapping Between beta Labels and the colors
alpha_label_color_map <- data.frame(
  alpha_label = alpha_labels,
  colors  = alpha_col
)

beta_label <- as.data.frame(beta_labels)
beta_label_color_map <- merge(beta_label, alpha_label_color_map, by.x = "beta_labels", by.y = "alpha_label", sort = F)

# Apply colors to the dendrograms
alpha_dend <- alpha_dend %>%
  set("labels_colors", adjustcolor(alpha_col, alpha = 0.5)) %>%  # Paler alpha labels
  set("labels_cex", 1) %>%
  set("branches_lwd", 1.5)

beta_dend <- beta_dend %>%
  set("labels_colors", beta_label_color_map$colors) %>%
  set("labels_cex", 1) %>%
  set("branches_lwd", 1.5)


dend_list <- dendlist(alpha_dend, beta_dend)
entanglement(dend_list)

dend_list_untangled <- untangle(dend_list, method = "step2side")
entanglement(dend_list_untangled)

# Assign colors to the connecting lines
line_colors <- labels_colors(dend_list_untangled[[1]]) %>%
  adjustcolor(alpha = 2)

# Save tanglegram as PDF
pdf("tanglegram_L4.pdf", width = 6, height = 8)
tanglegram(dend_list_untangled,
           lwd = 2,
           edge.lwd = 2,
           columns_width = c(5, 2, 5),
           margin_inner = 10,
           intersecting = T,
           axes = T,
           type = "r",
           lab.cex = 1,
           main_left = "    Alpha chains",
           main_right = "Beta chains  ",
           cex_main = 2,
           sub = "entanglement = 0.05",
           cex_sub = 1,
           match_order_by_labels = T,
           highlight_distinct_edges = F,
           common_subtrees_color_lines = F,
           color_lines = line_colors)
dev.off()
