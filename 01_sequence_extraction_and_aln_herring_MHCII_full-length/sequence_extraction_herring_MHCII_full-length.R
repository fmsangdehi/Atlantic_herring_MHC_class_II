library(Biostrings)
library(rtracklayer)
library(Rsamtools)
library(stringr)

# Define gene groups
gene_groups <- c("DAA", "DAB", "DBA", "DBB")

# List of GTF file paths and assembly sequences
gtf_files <- list.files(path = "path/to/gtf/files", pattern = "\\.gtf$", full.names = TRUE)
assembly_files <- list.files(path = "path/to/assemblies", pattern = "\\.fa$", full.names = TRUE)

# Set up output directory
output_dir <- "herring_MHCII_full-length_CDS_and_AA"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define output files for CDS and AA sequences for each gene group
output_files_cds <- setNames(file.path(output_dir, paste0("CDS_", gene_groups, ".fa")), gene_groups)
output_files_aa <- setNames(file.path(output_dir, paste0("AA_", gene_groups, ".fa")), gene_groups)

for (i in seq_along(gtf_files)) {
  assembly <- FaFile(assembly_files[i])
  assembly_name <- tools::file_path_sans_ext(basename(assembly_files[i]))
  
  gtf <- rtracklayer::import(gtf_files[i])
  
  gtf_filtered <- subset(gtf, type == "CDS" & !str_detect(mcols(gtf)$gene_id, "NullG$|PseudoG$|ψG"))
  
  for (gene_group in gene_groups) {
    gene_specific <- subset(gtf_filtered, str_starts(mcols(gtf_filtered)$gene_id, gene_group))
    gene_specific <- gene_specific[order(mcols(gene_specific)$gene_id, mcols(gene_specific)$exon_number)]
    
    lapply(split(gene_specific, mcols(gene_specific)$gene_id), function(x) {
      seq <- getSeq(assembly, x)
      full_cds <- DNAStringSet(paste(as.character(seq), collapse = ""))
      
      gene_id <- paste(assembly_name, sub("G$", "", unique(mcols(x)$gene_id)), sep = "_")
      names(full_cds) <- gene_id
      
      writeXStringSet(full_cds, filepath = output_files_cds[gene_group], append = TRUE, width = 10000)
      
      full_aa <- translate(full_cds)
      names(full_aa) <- gene_id  # Set the same name for AA
      writeXStringSet(full_aa, filepath = output_files_aa[gene_group], append = TRUE, width = 10000)
    })
  }
}

