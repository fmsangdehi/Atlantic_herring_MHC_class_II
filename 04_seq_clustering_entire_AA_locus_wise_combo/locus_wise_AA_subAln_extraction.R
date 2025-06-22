library(Biostrings)
library(stringr)

# Define input/output directories
input_dir <- "path/to/aln/dir"
output_dir <- "path/to/output/dir"

# Define input files
aln_files <- c("AA_aln_DAA.fa", "AA_aln_DAB.fa")

# Define patterns to extract and process
DAA_patterns <- c("DAA1", "DAA2", "DAA4")
DAB_patterns <- c("DAB1", "DAB2", "DAB4")

for (file in aln_files) {
  alignment <- readAAStringSet(paste0(input_dir, file))
  seq_names <- names(alignment)
  
  if (grepl("DAA", file)) {
    patterns <- DAA_patterns
  } else if (grepl("DAB", file)) {
    patterns <- DAB_patterns
  }
  
  for (pattern in patterns) {
    matches <- grepl(pattern, seq_names)
    
    filtered_sequences <- alignment[matches]
      
    fasta_filename <- paste0(output_dir, "AA_aln_", pattern, ".fa")
      
    writeXStringSet(filtered_sequences, fasta_filename)
  }
}
