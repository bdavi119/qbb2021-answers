# Week 6 Homework
______________________________________

# Setup/Getting Data

My working directory is `06-week-hw`, and inside that I placed the data in a `data` directory.

$`conda create -n hic cooltools cooler matplotlib numpy bedtools ucsc-bigWigToBedGraph`
$`conda activate hic`  
$`bigWigToBedGraph data/K562_hg19_H3K27me3_chr3.bw data/K562_hg19_H3K27me3_chr3.bg`

# Compartment Analysis

$`pip install jupyter` (Needed so I can use jupyter notebook from hic environment)

The rest of this was completed in `compartment_analysis.ipynb`.

# Expression vs. Repression
$`awk 'BEGIN{OFS="\t"}{$5=($3-$2)*$4; print $1,$2,$3,$4,$5}' data/K562_hg19_H3K27me3_chr3.bg > data/K562_hg19_H3K27me3_chr3_norm.bg`  
$`bedtools map -a data/K562_hg19_FPKM_chr3.bed -b data/K562_hg19_H3K27me3_chr3_norm.bg -c 5 -o sum > data/K562_hg19_FPKM_chr3_map.bed`  
$`awk 'BEGIN{OFS="\t"}{$7=$7 / ($3-$2); print $1,$2,$3,$4,$5,$6,$7}' data/K562_hg19_FPKM_chr3_map.bed > data/K562_hg19_FPKM_chr3_mapnorm.bed`  
$`bedtools map -a data/K562_hg19_FPKM_chr3_mapnorm.bed -b compartments.bed -c 5 -o distinct -f 0.5 > FPKM_with_compartments.bed`

THe final graphing was complete in `expression_analysis.ipynb`.
