library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(cowplot)

## Reroot the tree
unrooted_tree <- read.tree("ML_tree_beta.treefile")
rerooted_tree <- root(unrooted_tree, outgroup = "Nurse_shark_MHCIIB", resolve.root = TRUE)
write.tree(rerooted_tree, file = "ML_tree_beta_rerooted.nwk")


## Read in the rooted tree, and visualize with ggtree.
beta_tree <- read.newick("ML_tree_beta_rerooted.nwk", node.label='support')
beta_tree@phylo$tip.label <- gsub("'", "", beta_tree@phylo$tip.label)
beta_tree@phylo$tip.label <- paste0("  ", beta_tree@phylo$tip.label)


## Make the metadata and then ggtree.
metadata <- data.frame(tip_labels = beta_tree@phylo$tip.label) %>%
  mutate(category = ifelse(grepl("Herring", tip_labels), "Herring", "Non_herring"))

clades_df <- data.frame(clade = c("Tetrapods classical", "Tetrapods DO", "Tetrapods DM", "Teleosts DE", "Teleosts DA/DB", "Clupeid DA", "Clupeid DB", "Atlantic salmon DEA"),node = c(95, 100, 91, 77, 60, 62, 68, 35))

p1 <- ggtree(beta_tree, linetype=NA) %<+% metadata +
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
p2 <- flip(p1, 104, 98)

ggsave("ML_tree_beta_rooted_on_shark.pdf", plot = p2, device = "pdf", width = 8, height = 10)
