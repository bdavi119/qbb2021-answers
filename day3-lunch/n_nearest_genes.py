#!/usr/bin/env python3
"""
Brendon Davis
CMDB Bootcamp Day 3 Lunch - Advanced Exercises Part 2

This program requires 4 additional arguments from the command line. Extra command is the
number of nearest genes that you want.

This program does not yet work perfectly. It's close but doesn't print the nearest genes
in order. Further debugging required.
"""


import sys
import search
    

def main():
    """
    Parses .gtf file into ordered list of gene names and positions and then uses
    that to find genes that are closest to a postion of interest on a particular
    chromosome.
    """
    # Read in arguments
    gtf_file = sys.argv[1]
    mut_chrom = sys.argv[2]
    mut_pos = int(sys.argv[3])
    num_genes = int(sys.argv[4])
    f = open(gtf_file)

    genes = []
    for line in f:
        if line.startswith("#"):
            # Ignores header lines
            continue
        # Split line into a list and then define start and end
        fields = line.strip("\r\n").split("\t")
        start = int(fields[3])
        end = int(fields[4])
        
        # Use line only if it is gene on our chromosome of interest
        if (fields[0] == mut_chrom) and (fields[2] == "gene") and ('gene_biotype "protein_coding"' in line):
            subfields = fields[-1].split(';') # Take last field
            for field in subfields:
                if "gene_name" in field:
                    # Find 'gene name' and then take the string after
                    gene_name = field.split()[1]
            genes.append((gene_name, start, end)) # Builds list of tuples
    
    for i in range(num_genes):
        closest = search.binary_search(genes, mut_pos)
        idx = genes.index(closest[:3])
        genes.pop(idx)
        print("The number {} nearest gene is {}, which runs from positions {} to {}. "
               "It is {} base pairs from your position. This took {} "
               "iterations.\n".format(i+1, *closest))
        print()

    
    
if __name__ == '__main__':
    main()