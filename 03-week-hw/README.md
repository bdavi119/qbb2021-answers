# Week 3 homework

Below commands are those performed on the command line. Most analysis was done in Jupyter Notebook.
____________________________________________
## Problem 1
`plink --pca --vcf genotypes.vcf`

`plink.eigenvec` used for graphing principle components.

## Problem 2
`plink --freq --vcf genotypes.vcf`

`plink.frq` used for graphing allele frequencies.

## Problem 3
`plink --vcf genotypes.vcf --assoc --pheno CB1908_IC50.txt --covar plink.eigenvec --covar-number 1-10 --allow-no-sex`

`mv plink.qassoc CB1908.qassoc`

`plink --vcf genotypes.vcf --assoc --pheno GS451_IC50.txt --covar plink.eigenvec --covar-number 1-10 --allow-no-sex`

`mv plink.qassoc GS451.qassoc`

The .qassoc files were renamed so they would be recognizeable and so the second one would not overwrite the first. These files were used for subsequent graphing and analysis.

## Problem 4
`pip install qmplot`

`qmplot` was installed to use for making the Manhattan plots.

## Problem 5
Entirely on Jupyter Notebook.