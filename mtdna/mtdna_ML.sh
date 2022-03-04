NAME=bf186

# Create consensus sequences for NGS data
angsd -dofasta 3 -doCounts 10 -minQ 30 -i /steveData/cindy/bams/${NAME}.mtdna.sort.bam -out ${NAME}_mtDNA

# Concatenate fasta per sample in single file
# See Manuscript for aligning and GBlock cleanup

# IQTree
iqtree -s aligned_ebd_rogers.gb.phy -m MFP -bb 10000
