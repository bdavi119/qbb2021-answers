# Brendon Davis - Week 1 HW

Setup:
With asm.tgz in my week-01-hw directory, I entered `tar -xvf asm.tgz`. This created a folder called 'asm' in my directory with all the data.

# Question 1
## A
`cd asm` Note: nearly all of the rest of this is performed with asm as the working directory.   
`samtools faidx ref.fa -n 1`

The length of the reference sequence is 233,806 bases, which I found after opening the ref.fa.fai file that was created in asm/. The second column shows the length of the sequence.

## B 
`fastqc *.fq -extract`  
`less frag180.1_fastqc/fastqc_data.txt`  
`less frag180.2_fastqc/fastqc_data.txt`  
`less jump2k.1_fastqc/fastqc_data.txt`  
`less jump2k.2_fastqc/fastqc_data.txt`

From reading the beginning of each file, we see that frag180.1 and frag180.2 each have 35178 reads of length 100, while jump2k.1 and jump2k.2 each have 70355 reads of length 50. This confirms what we know from the README.

## C
I can expect around 60x coverage. Taking the total number of bases read divided by the length of the reference sequence, we get

![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B2%28100%29%2835178%29&plus;2%28150%29%2870355%29%7D%7B233806%7D%3D60.1828)

## D
`cat *.fq > merged.fq`  
`fastqc merged.fq -extract`  
`open merged_fastqc.html`  

I took a screenshot of the proper graph in the html file, displayed here: [Average Read Quality Graph](./avg_read_quality.png)

# Question 2
## A
`jellyfish count -m 21 -C -s 1000000 merged.fq`  
`jellyfish histo mer_counts.jf > reads.histo`  

There are 1030 21-mers that occur exactly 50 times, which I found by opening the reads.histo file and looking at the row with 50 in the first column.

## B
`jellyfish dump -c mer_counts.jf > dump.fa`  
`sort -nr -k 2 dump.fa | head -n 10`  

This printed the top 10 most frequent kmers and their counts:

GCCCACTAATTAGTGGGCGCC 104  
CGCCCACTAATTAGTGGGCGC 104  
CCCACTAATTAGTGGGCGCCG 104  
ACGGCGCCCACTAATTAGTGG 102  
AACAGGCCAGCTTATAAGCTG 100  
ACAGGCCAGCTTATAAGCTGG 99  
CAGGCCAGCTTATAAGCTGGC 96  
AGGCCAGCTTATAAGCTGGCC 94  
AGCATCGCCCACATGTGGGCG 82  
GCATCGCCCACATGTGGGCGA 80

## C
After uploading reads.histo to GenomeScope, it says the Genome Haploid Length is between 233,510 and 233,799 bp's. 

## D
The estimated genome size is very close to the length of the reference genome. The reference genome (233,806 bp's) is 7 bp's longer than the maximum estimated genome size (233,799 bp's) from GenomeScope.

# Problem 3
## A
`spades.py --pe1-1 frag180.1.fq --pe1-2 frag180.2.fq --mp1-1 jump2k.1.fq --mp1-2 jump2k.2.fq -o spades-output -t 4 -k 31`  
`grep -c '>' spades-output/contigs.fasta`

This printed '4', so 4 contigs were produced by Spades.

## B
`samtools faidx spades-output/contigs.fasta`  
`total=0; for i in $( cut -f 2 spades-output/contigs.fasta.fai ) ; do total=$( expr $total + $i ); done; echo $total`  

The above script printed 234,467, so that is the total length of the 4 contigs.

## C
`cut -f 2 spades-output/contigs.fasta.fai | sort -nr`

The largest number printed was 105,830, so this is the largest contig.

## D
`../calc_n50.py spades-output/contigs.fasta.fai 233806`

I wrote a short script (calc_n50.py) to calculate this. Not really necessary for just 4 contigs, but could be useful in the future. The script printed 47,860, so 47,860 bp's (or 47.86 kb) is the N50.

# Problem 4
## A
`mkdir dnadiff_output`  
`cd dnadiff_output`  
`dnadiff ../ref.fa ../spades-output/contigs.fasta`

I then opened the out.report file generated, and under '[Alignments]' it says that AvgIdentity is 100.00 for the reference.

## B
`show-coords dnadiff_output/out.delta`

This displayed the alignments made by dnadiff, and the first one is longest with 105,830 bp's.

## C
By opening dnadiff/out.report, I can see my assembly has 1 insertion. There are also 5 deletions, as indicated by the listing that there are 5 insertions in the reference compared to the query.

# Problem 5
## A
`show-coords dnadiff_output/out.delta`

The reference aligns up until position 26789, after which there is an insertion in the query.

## B
Looking at the alignments, there is a gap in the alignable query between 26787 and 27500 bp's exclusive on the NODE_3 sequence. This corresponds to a 712 bp insertion.

## C
`samtools faidx spades-output/contigs.fasta NODE_3_length_41351_cov_20.528098:26788-27499 > ../insertion.fasta`

The DNA sequence this produced was   CGCCCATGCGTAGGGGCTTCTTTAATTACTTGATTGACGCATGCCCCTCGTTCTACATGTCTAGCTTCGTAACTGCCCCGATTTATACAGGAGCATATGCGTTTCGTAGTGCCGGGAATGCATACCAAAGGGCTCACGGCGGGTACGCCACAATGGCTCAAGTCGAAAATGAATCGAAGACAACAAGGAATACCGTACCCAATTACTCAAGGACCTCATACACCATCCCATGCTACTTATCTACAGACATACACGCCAGCACCCAGCAACCAAAGCACACCGACGATAAGACTACAATCGCGATAAGCACAACTTACATTAGGAGGCCCGGCAAATCTTGACGGCGTTAAGTCCGACACGAATACCCCCCGACAAAAGCCTCGTATTCCGAGAGTACGAGAGTGCACAAAGCACCAAGGCGGGGCTTCGGTACATCCACCAGTAGTCCCGTCGTGGCGGATTTTCGTCGCGGATGATCCGAGGATTTCCTGCCTTGCCGAACACCTTACGTCATTCGGGGATGTCATAAAGCCAAACTTAGGCAAGTAGAAGATGGAGCACGGTCTAAAGGATTAAAGTCCTCGAATAACAAAGGACTGGAGTGCCTCAGGCATCTCTGCCGATCTGATTGCAAGAAAAAATGACAATATTAGTAAATTAGCCTATGAATAGCGGCTTTAAGTTAATGCCGAGGTCAATATTGACATCGGTA

## D
`./dna-decode.py --decode --input ../insertion.fasta > ../secret_message.txt`

The message is "Congratulations to the 2021 CMDB @ JHU class!  Keep on looking for little green aliens..."

Thanks for the nice message! I'm keeping my eyes peeled for the little green aliens.
