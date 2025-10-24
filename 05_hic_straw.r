## converting .hic matrices to sparse upper trianlge format
# .txt files for each chromosome

#load packages
library(strawr)

file_hic = file.path('/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output/Tests/als_hic/Sample_GWF2162/Sample_GWF2162.hic')

mat <- straw("NONE", file_hic, "chr1", "chr1", "BP", 100000)
write.table(mat, "/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output/Tests/als_hic/Sample_GWF2162/Sample_GWF216.NONE.chr1.100000.txt")


