# Week 5 Homework
_________________________________

# Part 1: ChIP-seq

## Setup
I copied the .tar data file to a `data` directory in `05-week-hw` and then unpacked the zip file within `data`. Then I copied the mm10 index file to `05-week-hw`, keeping the directory name `mm10`.

I started in the `base` environment and will switch to `chipseq` later.

## Mapping
$`mkdir mapped`  
$`bowtie2 -x mm10/mm10 -U data/CTCF_ER4.fastq > mapped/CTCF_ER4mapped.sam`  
$`bowtie2 -x mm10/mm10 -U data/CTCF_G1E.fastq > mapped/CTCF_G1Emapped.sam`  
$`bowtie2 -x mm10/mm10 -U data/input_ER4.fastq > mapped/input_ER4mapped.sam`  
$`bowtie2 -x mm10/mm10 -U data/input_G1E.fastq > mapped/input_G1Emapped.sam`

## Calling peaks
$`mkdir ER4-peaks G1E-peaks`  
$`conda activate chipseq`  
$`macs2 callpeak -t mapped/CTCF_ER4mapped.sam -c mapped/input_ER4mapped.sam --name=ER4 --gsize=61420004 --bdg --outdir ER4-peaks`  
$`macs2 callpeak -t mapped/CTCF_G1Emapped.sam -c mapped/input_G1Emapped.sam --name=G1E --gsize=61420004 --bdg --outdir G1E-peaks`  
$`conda deactivate`  
$`bedtools intersect -a ER4-peaks/*.narrowPeak -b all_chr19.bed  > ER4-peaks/filtered.narrowPeak`  
$`bedtools intersect -a G1E-peaks/*.narrowPeak -b all_chr19.bed  > G1E-peaks/filtered.narrowPeak`  

I looked online to find the size of mouse chromosome 19 is 61,420,004 (used in macs2 calls).  
Additionally, for the last two commands, I made a file called `all_chr19.bed` that has a single entry for the entire length of chromosome 19. This allows me to filter by chromosome, because I noticed an aberrant entry for another chromosome in one fo the .narrowPeak files.


## Differential binding
Get binding lost from G1E->ER4 differentiation:  
$`bedtools intersect -v -bed -a G1E-peaks/filtered.narrowPeak -b ER4-peaks/filtered.narrowPeak > diff/lost_binding.narrowPeak`

Get binding gained from G1E->ER4 differentiation:  
$`bedtools intersect -v -bed -a ER4-peaks/filtered.narrowPeak -b G1E-peaks/filtered.narrowPeak > diff/gained_binding.narrowPeak`


## Feature overlapping
$`wget https://raw.githubusercontent.com/bxlab/qbb2020/master/week5/Mus_musculus.GRCm38.94_features.bed`  
$`bedtools intersect -c -a Mus_musculus.GRCm38.94_features.bed -b G1E-peaks/filtered.narrowPeak > feature-count/G1E-overlap.bed`  
$`bedtools intersect -c -a Mus_musculus.GRCm38.94_features.bed -b ER4-peaks/filtered.narrowPeak > feature-count/ER4-overlap.bed`  

## Plotting
Performed in jupyter notebook.

# Part 2: Motif Discovery
## Data
$`mkdir motif-data`  
$`cd motif-data/`  
$`wget http://jaspar2018.genereg.net/download/CORE/JASPAR2018_CORE_non-redundant_pfms_meme.zip`  
$`unzip JASPAR2018_CORE_non-redundant_pfms_meme.zip`  
$`cd ..`  

## Motif finding
$`sort -nrk 5 ER4-peaks/filtered.narrowPeak | head -n 100 > enriched.narrowPeak`  
$`bedtools getfasta -fi /Users/cmdb/data/genomes/mm10.fa -bed enriched.narrowPeak -fo seqs_for_motif_finding.fasta`  
$`conda activate chipseq`  
$`meme-chip -oc motifs-discovered -maxw 20 seqs_for_motif_finding.fasta`  
$`tomtom -png -oc tomtom-out motifs-discovered/meme_out/meme.txt motif-data/MA*`  
$`cp tomtom-out/align_RCCACYAGRKGGCRS_134_+MA0139.1.png .`  
$`mv tomtom-out/align_RCCACYAGRKGGCRS_134_+MA0139.1.png best_match.png`  

Looking at the three matches found, MA0139.1 gave the best p-value (6.93e-26), so I saved it's logo comparison in the main working directory as `best_match.png`.

## Submitting

Later I performed the following commands to get copies of the narrowPeak files for turning in:  
$`mkdir bed-files-for-submit`  
$`cp G1E-peaks/filtered.narrowPeak bed-files-for-submit/G1E_peaks.narrowPeak`  
$`cp ER4-peaks/filtered.narrowPeak bed-files-for-submit/ER4_peaks.narrowPeak`  
$`cp diff/* bed-files-for-submit/`  



