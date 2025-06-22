library(tidyverse)
library(gggenes)
require(RColorBrewer)
library(Cairo)
library(cowplot)

genes_in <- read.table("Locus4_supertypes.txt", sep = "\t", header=T, check.names = FALSE)

genes <- genes_in %>% 
  pivot_longer(cols = "0":"17", names_to = "x_axis", values_to = "genes") %>%
  filter(genes != "") %>%
  mutate(x_axis = as.numeric(x_axis)) %>%
  separate(col = genes, into = c("gene_name", "gene_group", "strand", "null"), sep = "_", remove = F) %>%
  mutate(start = ifelse(strand == "1", x_axis - 0.9, x_axis + 0.9)) %>%
  mutate(end = ifelse(strand == "1", x_axis + 0.9, x_axis - 0.9)) %>%
  mutate(status = ifelse(!is.na(null), "Non-functional", gene_group))

genes <- genes %>%
  mutate(null = ifelse(null == "pg", "Ψ", null)) %>%
  mutate(gene_name = ifelse(!is.na(null), paste0(gene_name, null), gene_name)) %>%
  mutate(gene_name_color = ifelse(grepl("DAB|DBB", gene_name) & !grepl("N", gene_name), "white", "black"))

genes$color <- factor(genes$gene_group, levels = c("4aA", "4aB", "4bA", "4bB", "4cA", "4cB", "DBA", "DBB"))
genes$fill <- factor(genes$status, levels = c("4aA", "4aB", "4bA", "4bB", "4cA", "4cB", "DBA", "DBB", "Non-functional"),
                     labels = c("4a", "4a ", "4b", "4b ", "4c", "4c ", "DB", "DB ", "Non-functional"))


outline_color_brewer <- c("#99aeef", "#4169e1", "#FFB684", "#ff6900", "#EA94F7", "#CE16EB", "#c4c1c1", "#A6A6A6")
fill_color_brewer <- c("#99aeef", "#4169e1", "#FFB684", "#ff6900", "#EA94F7", "#CE16EB", "#c4c1c1", "#A6A6A6", "#FFFFFF")


p <- ggplot(as.data.frame(genes), aes(xmin = start, xmax = end, y = y_axis, color = color, fill = fill, label = I(gene_name))) +
  geom_segment(mapping = aes(x = x_min, y = y_axis, xend = x_max, yend = y_axis), linewidth=0.5, color = "gray50") +
  geom_gene_arrow(size = 0.8) +
  geom_gene_label(color = genes$gene_name_color, align = "centre", fontface = "bold", family = "Arial", grow = F, reflow = F, min.size = 1, size = 6) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = 1:nrow(genes_in), labels = rev(genes_in$SH)) +
  scale_color_manual(values = outline_color_brewer) +
  scale_fill_manual(values = fill_color_brewer) +
  labs(fill = "Alpha chain \nBeta chain" ) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 10, face = "bold", color = "gray20", hjust = 0, family = "Arial"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 15, 0, 0),
        legend.title = element_text(size = 12, face = "bold", color = "gray20", hjust = 0, lineheight = 1.4, family = "Arial"),
        legend.direction = "horizontal",
        legend.byrow = F,
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14, face = "bold", color = "gray20", hjust = 0, family = "Arial"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  guides(color = "none", fill = guide_legend(nrow = 2))


final_plot <- ggdraw() +
  draw_plot(p, x = 0.01, y = 0, width = 0.99, height = 0.97) +
  draw_label("b", x = 0.02, y = 1, hjust = 0, vjust = 1.07, fontfamily = "Arial", fontface = "bold", size = 22)

CairoPDF("structural_haplotypes_L4.pdf", family = "Arial", width = 7.7, height = 8)
print(final_plot)
dev.off()
