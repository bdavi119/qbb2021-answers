#!/bin/bash

# dir="./test9/*.fastq"
# Above directory used for testing with one file

dir="./reads/*.fastq"
ref="ref/sacCer3.fa"

# If file is run with -b flag, it generates bam files (steps 2 and 3),
# which takes a long time. If no flag is given, it assumes the bam files
# already exist, and the program can run much more quickly.
get_bams=0
while getopts 'b' flag
do
	case "${flag}" in
		b) get_bams=1
	esac
done

if [ $get_bams -eq 1 ]
then
	for file in $dir
	do
		# Step 2
		NUM=$(echo $file | cut -c13-14)
		sam_file="02step/aln_$NUM.sam"
		bwa mem -t 4 -R "@RG\tID:Sample$NUM\tSM:Sample$NUM" ref/sacCer3.fa $file > $sam_file
	
		# Step 3
		bam_sorted="03step/sorted_$NUM.bam"
		samtools sort -O bam -o $bam_sorted $sam_file
		samtools index $bam_sorted
	done
fi
	
#Step 4
step04="04step/04.vcf"
freebayes -f $ref -p 1 --genotype-qualities 03step/*.bam > $step04

#Step 5
step05="05step/05.vcf"
vcffilter -f "QUAL > 20" $step04 > $step05

# Step 6
step06="06step/06.vcf"
vcfallelicprimitives -k -g $step05 > $step06

# Step 7
final_all="final_all.vcf"
snpeff ann R64-1-1.99 $step06 > $final_all

# Make shortened vcf file to turn in
final_1000="final_1000.vcf"
head -n 1000 $final_all > $final_1000