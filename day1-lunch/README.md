1c)
cp ~/qbb2021/data/*.bed ~/qbb2021/data/dm6.chrom.sizes .

2a)
wc -l *.bed > feature_count.txt

2c) Observations:
Of the three methylation types, K27me3 and K4me3 combined occur more frequently than K9me3, which may suggest that, in general, more genes are being activated than repressed.q
The total number of methylation sites is about half of the number of genes, indicating at least half of all genes (and likely more) have no associated methylation sites of these types.

3a)
sort fbgenes.bed | cut -f 1 | uniq -c > fbgenes.info

3c) Observations:
Chromosome Y has a remarkably small number of genes, especially compared to chromosome X.
The vast majority of genes are held on chromosomes 2, 3, and X. Chromome 4 has ~5% the number of genes as any of these bigger cromosomes.

4a)
bedtools intersect -u -a fbgenes.bed -b K9me3.bed | sort | cut -f 1 | uniq -c > chr-with-fbgenes-k9.txt

4d) Observations:
~90% of the genes on Chromosome 4 have K9me3 methylation, suggesting Chr4 is mostly under repressive control.
About 1/5th of genes on each of ChrX and ChrY are under K9me3 control, so even though they have stark size differences, they are proportionally under the same amount of control.

