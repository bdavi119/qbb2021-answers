# Week 4 Homework

## Data setup
Placed data files in my `04-week-hw` directory with `chr6.fa.gz` in an `index` subdirectory. I also made an empty `fq-files` and `fq-output` for later.  

Then I ran:  
$`tar xzf methylation_fastq.tar`  
$`mv ./*.fastq fq-files`  
$`rm *.tar`  
$`fastqc fq-files/* -o fq-outputs/`  
$`open fq-outputs/SRR3083926_1.chr6_fastqc.html`

This opened the html document for the first set of reads. I thought it was most interesting that the GC content was so low, at only 21%. This means that the area that was read was heavily enriched in A's and T's.

## Bisulfite mapping
### Initial mapping
$`bismark_genome_preparation --genomic_composition --parallel 7 index/`  
$`mkdir bismark-outputs`  
$`bismark -B condn926 -o bismark-outputs/ index/ -1 fq-files/SRR3083926_1.chr6.fastq -2 fq-files/SRR3083926_2.chr6.fastq`  
$`bismark -B condn929 -o bismark-outputs/ index/ -1 fq-files/SRR3083929_1.chr6.fastq -2 fq-files/SRR3083929_2.chr6.fastq`  
$`mkdir bismark-outputs`  

This created two .bam files, one for each condition, and they were put in the `bismark-outputs` directory (for organization).

### Deduplicate
$`mkdir deduplicated`  
$`deduplicate_bismark --output_dir deduplicated/ bismark-outputs/*.bam`

### Sorting
$`mkdir sorted-reads`  
$`samtools sort -o sorted-reads/sorted926.bam deduplicated/condn926_pe.deduplicated.bam`  
$`samtools sort -o sorted-reads/sorted929.bam deduplicated/condn929_pe.deduplicated.bam`  
$`samtools index sorted-reads/sorted926.bam`  
$`samtools index sorted-reads/sorted929.bam`  

### Methylation extraction 
$`bismark_methylation_extractor --bedgraph --comprehensive -o extracted deduplicated/*.bam`

### IGV visualization 
$`igv sorted-reads/*.bam extracted/*bedGraph*`


## Analysis
Command from assignment:  
$`awk 'BEGIN{OFS="\t"}{if ($4 == "+") print $3,$5 - 2000,$5,$13,$12,$4; else print $3,$6,$6 + 2000,$13,$12,$4;}' mm10_refseq_genes_chr6_50M_60M.bed | grep -v Rik | uniq -f 3 | sort -k2,2n > promoters.bed`  

Then I used `bedtools map` to get the total methylation score for each promoter.

$`bedtools map -a promoters.bed -b extracted/condn926_pe.deduplicated.bedGraph.gz -c 4 -o sum > mapped/mapped926.tsv`  
$`bedtools map -a promoters.bed -b extracted/condn929_pe.deduplicated.bedGraph.gz -c 4 -o sum > mapped/mapped929.tsv`

The rest of the analysis was completed in a jupyter notebook.
