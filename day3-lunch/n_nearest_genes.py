#!/usr/bin/env python3
"""
Brendon Davis
CMDB Bootcamp Day 3 Lunch - Advanced Exercises Part 2

This program requires 4 additional arguments from the command line: gtf file, chromosome of interest,
position of interest, and the number of nearest genes that you want.
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

    genes = search.gtf_parser(gtf_file, mut_chrom)
    
    for i in range(num_genes):
        # Performs same thing as in main() of search.py except does it multiple times
        # once for each gene we want to find
        closest = search.binary_search(genes, mut_pos)
        idx = genes.index(closest[:3])
        genes.pop(idx)
        print("The number {} nearest gene is {}, which runs from positions {} to {}. "
               "It is {} base pairs from your position. This took {} "
               "iterations.".format(i+1, *closest))
    
    
if __name__ == '__main__':
    main()