library(Biostrings)
library(ape)
require(stringdist)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(grid)

aln_dir <- "path/to/aln/dir"
out_dir <- "path/to/output/dir"

aln_file <- c("DAA", "DAB", "DBA", "DBB")

for(i in seq_along(aln_file)){
  aa_aln <- Biostrings::readAAStringSet(paste0(aln_dir, "AA_aln_", aln_file[i], "_Mafft_orderAligned.fa"), format="fasta")
  seq_length <- unique(width(aa_aln))
  
  if(length(aa_aln) > 3){
    aa_vec <- as.character(aa_aln)
    AA_dist <- stringdistmatrix(aa_vec, method = "hamming", useNames = "names") * 100 / seq_length
    dist_mat <- as.matrix(AA_dist)
    
    row_order <- column_order <- hclust(AA_dist, method = "average")$order
#    dist_mat <- dist_mat[row_order,column_order]


    ## Make a heatmap
    dend <- as.dendrogram(hclust(AA_dist, method = "average"))

    Loci_annot <- data.frame(Loci = paste(" ", substr(gsub("^.*_", "", row.names(dist_mat)), 4, 4), sep = ""))
    annot_col <- list(Loci=c("1"="#E41A1C", "2"="#377EB8", "3"="#4DAF4A", "4"="#984EA3", "5"="#FF7F00", "6"="#FFFF33", "7"="#A65628", "8" = "#F781BF", "9" = "#999999"))
    names(annot_col$Loci) <- paste0(" ", names(annot_col$Loci))
    
    complex_heatmap <- Heatmap(dist_mat, name = "% amino acid differences",
                               cluster_rows = dend, cluster_columns = dend,
                               row_order = NULL, column_order = NULL,
                               col = colorRamp2(seq(min(dist_mat), max(dist_mat), length = 2), c("white", "blue"), space = "LAB"),
                               show_row_dend = F, show_column_dend = T,
                               column_dend_height = unit(4, "cm"), column_dend_side = "top",
                               show_heatmap_legend=T,
                               heatmap_legend_param = list(at = c(round(min(dist_mat), 0), round(max(dist_mat), 0)), labels = paste(" ", c(round(min(dist_mat), 0), round(max(dist_mat), 0)), sep=""), legend_height = unit(4, "cm"), grid_width = unit(0.5, "cm"), title = gt_render("% amino acid<br>differences", r = unit(2, "pt"), padding = unit(c(0, 0, 1, 0), "mm")), title_position = "topcenter", title_gap = unit(2.5, "mm"), title_gp = gpar(fontsize = 16, fontface = "bold", margin = margin(b = 3)), labels_gp = gpar(fontsize = 16)),
                               show_row_names = F, show_column_names = F,
                               row_names_side = "left", column_names_side = "top",
                               row_names_max_width = unit(1.7, "cm"), column_names_max_height = unit(1.7, "cm"),
                               top_annotation = HeatmapAnnotation(df=Loci_annot, name="", col=annot_col, show_legend=T, gap=unit(1,"mm"), show_annotation_name = T, annotation_name_side = "right", annotation_name_offset = unit(0.1, "cm"), annotation_name_gp = gpar(fontsize = 13, fontface = "bold", col = "white"), simple_anno_size = unit(0.5, "cm"),
                                                                annotation_legend_param = list(grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"), title = gt_render("Loci", r = unit(2, "pt"), padding = unit(c(0, 0, 0.5, 0), "mm")), title_position = "topcenter", title_gap = unit(2.5, "mm"), title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 16))),
                               left_annotation = HeatmapAnnotation(df=Loci_annot, name="", col=annot_col, show_legend=F, gap=unit(1,"mm"), show_annotation_name = F, which="row", simple_anno_size = unit(0.5, "cm")),
    )
    
    pdf(paste0(out_dir, gsub(".fa", '', aln_file[i]), "_heatmap", ".pdf"), width = 8.5, height = 8)
    draw(complex_heatmap, padding = unit(c(3, 3, 4, 14), "mm"), adjust_annotation_extension = T, legend_gap = unit(10, "mm"), merge_legend = TRUE)
    dev.off()
  }
}
