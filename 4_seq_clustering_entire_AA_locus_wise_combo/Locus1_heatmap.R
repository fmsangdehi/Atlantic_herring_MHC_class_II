library(Biostrings)
library(ape)
require(stringdist)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)

aln_dir <- "path/to/locus-wise/aln/dir"
output_dir <- "path/to/output/dir"

aln_file_alpha <- "AA_aln_DAA1.fa"; row_title <- "Alleles of DAA1.1, DAA1.2 and DAA1.3 genes"
aln_file_beta <- "AA_aln_DAB1.fa"; col_title <- "Alleles of DAB1.1, DAB1.2 and DAB1.3 genes"
Locus <- "L1"

## Import the alignment
aa_aln_alpha <- Biostrings::readAAStringSet(paste(aln_dir, aln_file_alpha, sep = ""), format="fasta")
seq_length_alpha <- unique(width(aa_aln_alpha))
aa_aln_beta <- Biostrings::readAAStringSet(paste(aln_dir, aln_file_beta, sep = ""), format="fasta")
seq_length_beta <- unique(width(aa_aln_beta))

# Reorder DAB sequences based on DAA sequence order, to ensure DAB sequences are exactly in the same order as DAA
DAA_seq_names <- names(aa_aln_alpha)
expected_DAB_seq_order <- gsub("DAA", 'DAB', DAA_seq_names)
aa_aln_beta <- aa_aln_beta[expected_DAB_seq_order]

## Calculate and order edit distance matrix triangle for beta
    aa_vec_beta <- as.character(aa_aln_beta)
    AA_dist_beta <- stringdistmatrix(aa_vec_beta, method = "hamming", useNames = "names") * 100 / seq_length_beta
    order_beta = hclust(AA_dist_beta, method = "average")$order
    dist_mat_beta <- as.matrix(AA_dist_beta)
    dist_mat_beta <- dist_mat_beta[order_beta,order_beta]
    
## Calculate and order edit distance matrix triangle for alpha
    aa_vec_alpha <- as.character(aa_aln_alpha)
    AA_dist_alpha <- stringdistmatrix(aa_vec_alpha, method = "hamming", useNames = "names") * 100 / seq_length_alpha
    dist_mat_alpha <- as.matrix(AA_dist_alpha)
    dist_mat_alpha <- dist_mat_alpha[order_beta,order_beta]

## Combine alpha and beta edit distance matrices
    dist_mat <- dist_mat_alpha
    dist_mat[upper.tri(dist_mat, diag = F)] <- dist_mat_beta[upper.tri(dist_mat_beta, diag = F)]
    colnames(dist_mat) <- gsub("DAA", "DAB", colnames(dist_mat))

## Extract the raw matrix
    write.csv(dist_mat, paste(output_dir, Locus, "_combined_mtrx", ".csv", sep = ""))
    
dend <- as.dendrogram(hclust(AA_dist_alpha, method = "average"))

## Annotation for L2
Locus_annot <- data.frame(Group = rep("1a",63))
annot_col_1 <- list(Group=c("1a"="#7139c4"))


## Make a heatmap, with annotation
    complex_heatmap <- Heatmap(dist_mat, name = "Number of amino acid changes",
        row_title = gt_render(row_title, r = unit(2, "pt"), padding = unit(c(0, 0, 1, 0), "mm")), row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title = gt_render(col_title, r = unit(2, "pt"), padding = unit(c(0, 0, 1, 0), "mm")), column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        cluster_rows = F, cluster_columns = F,
        row_order = NULL, column_order = NULL,
        col = colorRamp2(seq(min(dist_mat), max(dist_mat), length = 2), c("white", "blue"), space = "LAB"),
        show_row_dend = T, show_column_dend = F,
        border = T, border_gp = gpar(col = "black", lwd = 2),
        show_heatmap_legend=T,
        heatmap_legend_param = list(at = c(round(min(dist_mat), 0), round(max(dist_mat), 0)), legend_height = unit(4, "cm"), grid_width = unit(0.5, "cm"), title = gt_render("% amino acid<br>differences", r = unit(2, "pt"), padding = unit(c(0, 0, 2, 0), "mm")), title_position = "topcenter", title_gap = unit(2.5, "mm"), title_gp = gpar(fontsize = 12, fontface = "bold", margin = margin(b = 3)), labels_gp = gpar(fontsize = 12)),
        show_row_names = T, show_column_names = T,
        row_names_side = "left", column_names_side = "top",
        row_names_max_width = unit(1.7, "cm"), column_names_max_height = unit(1.7, "cm"),
        column_names_gp = gpar(fontsize = 4+50/nrow(dist_mat)),
        row_names_gp = gpar(fontsize = 4+50/nrow(dist_mat)),
        top_annotation=HeatmapAnnotation(df=Locus_annot, name="", col=annot_col_1, show_legend=F, gap=unit(1,"mm"), show_annotation_name = T, annotation_name_side = "right", annotation_name_offset = unit(0.2, "cm"), annotation_name_gp = gpar(fontsize = 5, fontface = "bold", col = "white"), simple_anno_size = unit(0, "cm"),
                                         annotation_legend_param = list(grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"), title = gt_render("Beta", r = unit(2, "pt"), padding = unit(c(0, 0, 1.2, 0), "mm")), title_position = "topcenter", title_gap = unit(2.5, "mm"), title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12))),
    )
        
    
    pdf(paste(output_dir, Locus, "_combined_heatmap", ".pdf", sep = ""), width = 8.5, height = 7)
    draw(complex_heatmap, padding = unit(c(4, 4, 4, 12), "mm"), adjust_annotation_extension = T, legend_gap = unit(10, "mm"), merge_legend = F)
    decorate_heatmap_body("Number of amino acid changes", {
      grid.lines(c(0, 1), c(1, 0), gp = gpar(col = "black", lwd = 2))
    })
    grid.text("c", x = unit(3, "mm"), y = unit(1, "npc") - unit(0.5, "mm"),
              just = c("left", "top"), gp = gpar(fontsize = 32, fontface = "bold"))
    dev.off()
