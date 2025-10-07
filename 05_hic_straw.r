## converting .hic matrices to sparse upper trianlge format
# .txt files for each chromosome

#load packages
suppressPackageStartupMessages({
  library(strawr)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript hic_compare_all_chromosomes.R <Sample1> <Resolution>")
}
sample1 <- args[1]
resolution <- as.integer(args[2])




chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

hic.list <- lapply(chromosomes, function(chr) {
  message("Processing: ", chr)
  mat1 <- straw("NONE", sample1, chr, chr, "BP", resolution)  
  create.hic.table(mat1, mat2, chr = chr, exclude.regions = exclude, exclude.overlap = 0.2)
})
names(hic.list) <- chromosomes