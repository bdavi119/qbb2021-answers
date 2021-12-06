Regrade 12/6/21 - SC 
Great job, Brendon, thanks! 4/4

This is really good! One small note: for your expression histogram, it looks like you're plotting H3K27me3 presence rather than gene expression, and you're choosing your cutoff based on that. You should be plotting a histogram of gene expression (i.e. column 4 in your data frame), not H3K27me3 presence, and choosing your cutoff based on that (-0.25)

3.75/4



From Brendon: I fixed the issue, so now the histogram shows gene expression. I also changed the expression cutoff and used that in the subsequent code blocks, and this yields results for the violin plot that make more sense than what I had previously (ie. majority of B genes are off, majority of A genes are on).
