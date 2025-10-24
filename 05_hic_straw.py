#convert .hic to usable sparse matrix for downstrream analysis using straw
# doesn't work because of operating system on cluster not compatible with straw
import hicstraw

SAMPLE1 = '/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output/Tests/als_hic/Sample_GWF2162/Sample_GWF2162.hic'

result = hicstraw.straw('NONE', SAMPLE1, 'X', 'X', 'BP', 1000000)
# the values returned are in x / y / counts
for i in range(len(result[0])):
   print("{0}\t{1}\t{2}".format(result[0][i], result[1][i], result[2][i]))

