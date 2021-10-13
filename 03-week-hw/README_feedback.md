- QQ plot looks good - just wanted you to know you can also use log10pvalues for the axes instead of quantiles 
- Your MAF histogram looks like the y-axis values are doubled... try again? (-0.5) 
- Manhattan correctly has significant p-values tagged in red, but there are additional red dots along the bottom of the chart... try again (-0.5) 
- 6/7


From Brendon:

- I removed the bins parameter from the histogram so now it is the default number of bins, which means the y-values have changed. I am not sure if this fixes the problem but my histogram now looks like what other students got.
- I have added a aparmater to the manhattan plots to remove the extra dots. The plot was originally trying to mark other SNPs that may be linked to significant SNPs, but now it only marks significant SNPs.
