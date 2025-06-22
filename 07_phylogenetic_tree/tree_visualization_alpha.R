library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(cowplot)

## Reroot the tree
unrooted_tree <- read.tree("ML_tree_alpha.treefile")
rerooted_tree <- root(unrooted_tree, outgroup = "Nurse_shark_MHCIIA", resolve.root = TRUE)
write.tree(rerooted_tree, file = "ML_tree_alpha_rerooted.nwk")


## Read in the rooted tree, and visualize with ggtree.
alpha_tree <- read.newick("ML_tree_alpha_rerooted.nwk", node.label='support')
alpha_tree@phylo$tip.label <- gsub("'", "", alpha_tree@phylo$tip.label)
alpha_tree@phylo$tip.label <- paste0("  ", alpha_tree@phylo$tip.label)


## Make the metadata and then ggtree.
metadata <- data.frame(tip_labels = alpha_tree@phylo$tip.label) %>%
  mutate(category = ifelse(grepl("Herring", tip_labels), "Herring", "Non_herring"))

clades_df <- data.frame(clade = c("Tetrapods classical", "Tetrapods DO", "Tetrapods DM", "Teleosts DE", "Teleosts DA/DB", "Clupeid DA", "Clupeid DB", "Atlantic salmon DEA"),node = c(89, 93, 85, 81, 58, 60, 65, 30))

p1 <- ggtree(alpha_tree, linetype=NA) %<+% metadata +
  geom_highlight(data=clades_df, aes(node=node, fill=clade), alpha=1, align="right", extend=0.03, type = "gradient", gradient.direction = 'tr') +
  geom_tree(linewidth = 1) +
    geom_tiplab(aes(color = category), align = T, size=3.5, linetype = "dotted", linesize = 0.3, hjust = 0, offset=0) +
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("Clupeid DA"="#FFB199", "Clupeid DB"="#80CCFF", "Teleosts DA/DB"="#B6EDB0", "Teleosts DE"="#CCCCCC", "Atlantic salmon DEA"="#CCCCCC", "Tetrapods classical"="#FFA6C9", "Tetrapods DO"="#C7B3FF", "Tetrapods DM"="#FFE68C")) +
  geom_nodelab(aes(label = support), hjust=-.25, size = 2.5) +
  geom_treescale(y = -0.2, linesize = 1, width = NULL, offset = 0.3) +
  xlim(0, 5) +
  theme(legend.position = "None",
        plot.margin = margin(2, 0, 2, 0, "pt")
  )
p2 <- rotate(p1, 84)
p3 <- flip(p2, 92, 95)

ggsave("ML_tree_alpha_rooted_on_shark.pdf", plot = p3, device = "pdf", width = 8, height = 10)
