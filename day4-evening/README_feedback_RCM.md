Excellent work--especially on the encoding of the t-test as a linear model and the natural extension to Poisson GLM.

Very minor critique of the histogram: in this case, holding the number of bins constant between the mothers and fathers leads to a situation where the width of the bins differs between the groups, which could be visually misleading. Instead of fixing the number of bins at 30, can you fix the width of the bins to be equal?

Brendon: Fixed the bin sizes. I used an algorithm from site below to set a bin number for each dataset that creates the same width for every bin.
https://www.kite.com/python/answers/how-to-set-the-bin-size-of-a-matplotlib-histogram-in-python
