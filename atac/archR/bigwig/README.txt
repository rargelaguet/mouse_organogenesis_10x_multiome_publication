############
## Signac ##
############

https://timoast.github.io/sinto/basic_usage.html#filter-cell-barcodes-from-bam-file
https://github.com/timoast/sinto


Hi, it is not currently possible to create a bigwig for different groups of cells in Signac. I'd suggest writing the cell names to a file and then splitting the bam file by cell using the sinto package (https://github.com/timoast/sinto), and then creating bigwig tracks using deeptools (https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)


###########
## ArchR ##
###########

getGroupBW()