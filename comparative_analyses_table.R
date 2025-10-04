# author: SXL

#******************
#* Load libraries
#******************
suppressPackageStartupMessages({
  library(strawr)
  library(HiCcompare)
  library(ggplot2)
  library(Cairo)
})

#******************
#* Parse arguments
#******************
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript hic_compare_all_chromosomes.R <Sample1> <Sample2> <Resolution>")
}
sample1 <- args[1]
sample2 <- args[2]
resolution <- as.integer(args[3])

#******************
#* Set file paths
#******************
base_path <- "/gpfs/commons/home/sliu/ALS_MicroC_analyses/output/04_hic_per_sample"
hic_file1 <- file.path(base_path, paste0("Sample_", sample1), paste0("Sample_", sample1, ".hic"))
hic_file2 <- file.path(base_path, paste0("Sample_", sample2), paste0("Sample_", sample2, ".hic"))

#******************
#* Load blacklist
#******************
data("hg38_blacklist")
exclude <- hg38_blacklist

#******************
#* Define chromosomes
#******************
chromosomes <- c(paste0("chr",1:22), "chrX", "chrY")

#******************
#* Create hic.tables for all chromosomes
#******************
hic.list <- lapply(chromosomes, function(chr) {
  message("Processing: ", chr)
  mat1 <- straw("NONE", hic_file1, chr, chr, "BP", resolution)
  mat2 <- straw("NONE", hic_file2, chr, chr, "BP", resolution)
  create.hic.table(mat1, mat2, chr = chr, exclude.regions = exclude, exclude.overlap = 0.2)
})
names(hic.list) <- chromosomes

#******************
#* Apply LOESS normalization in parallel
#******************
hic.list <- hic_loess(hic.list, parallel = TRUE)

#******************
#* Perform comparison and save MD plots
#******************

for (chr in chromosomes) {
  message("Comparing: ", chr)
  plot_file <- paste0("hic_compare_MD_plot_", sample1, "_vs_", sample2, "_", chr, ".png")
  CairoPNG(filename = plot_file, width = 1000, height = 600, res = 150)
  hic.list2 <-  hic_compare(hic.list[[chr]], A.min = 15, adjust.dist = TRUE, p.method = "fdr", Plot = TRUE)
 # hic.list <- do.call(rbind, hic.list)
#  hic.list <- hic.list[hic.list$p.adj<0.05,]
  table_name  <- paste0("hic_compare_", sample1, "_vs_", sample2, "_", chr, ".csv")
  write.table(hic.list2, table_name, sep = ",")
  dev.off()
}

# note SXL: if saving plots don't work with png then next is to try save pdf from each loop as hic_compare naturally saves a pdf.
# or save the hic_compare tables  (try to understand what the information can do , if it tells yuou which region of the chromosome it is)

