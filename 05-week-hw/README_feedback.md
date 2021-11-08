Hi Brendon, 

--Your motif sequence doesn't look quite right (-1). I think you used the correct file (the narrowPeak file, enriched for just the top 100), 
but review the rest of the steps, i.e., 
$bedtools getfasta -fi /Users/cmdb/data/genomes/mm10.fa -bed enriched.narrowPeak -fo seqs_for_motif_finding.fasta
$conda activate chipseq
$meme-chip -oc motifs-discovered -maxw 20 seqs_for_motif_finding.fasta
$tomtom -png -oc tomtom-out motifs-discovered/meme_out/meme.txt motif-data/MA*

5/6
