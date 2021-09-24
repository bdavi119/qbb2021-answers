# File structure
02-week-hw/  
|  
|---files for solution  
|  
|---reads/  
|   |  
|   |---the 10 fastq files with reads  
|  
|---ref/  
|   |  
|   |---reference genome, including indexing files  

Steps below are done in `02-week-hw/` unless otherwise specified. Also, all the necessary directories (02step, 03step, 04step, 05step, 06step) already exist in this directory as well.

# Step 1
`bwa index ref/sacCer3.fa`

# Step 2-7
`snpeff download R64-1-1.99` Required for step 7  
`./var_caller.sh`

I wrote a bash script to process all ten of the .fasta files. This generates output
for each step in directories called ##step (## replaced with the step number). At the end of running, it produced `final_all.vcf` and `final_1000.vcf` in the working directory.

# Step 8
All of the data wrangling and graphing for step 8 was done in the jupyter notebook `var_grapher.ipynb`. The final graph is saved as `graphs.jpg`.