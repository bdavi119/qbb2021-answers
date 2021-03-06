#!/usr/bin/env python3
"""
Brendon Davis
CMDB Bootcamp Day 3 Lunch

This program requires 3 additional arguments from the command line: gtf file, chromosome of interest,
position of interest, and the number of nearest genes that you want.
"""


import sys

def distance_to_gene(gene, pos):
    """
    Takes a list of genes as tuples and a position of interest as an integer. Returns the
    distance from the position of interest to the nearest end of the gene.
    """
    if (pos >= gene[1]) and (pos <= gene[2]):
        # If position is in gene
        return 0
    elif pos < gene[1]:
        return gene[1] - pos
    else:
        return pos - gene[2]


def encloses(gene1, gene2):
    """
    Takes two tuples representing genes and returns True if the first gene encloses 
    the second gene and False if it does not.
    """
    a, gene1_start, gene1_end = gene1
    b, gene2_start, gene2_end = gene2
    
    if (gene1_start <= gene2_start) and (gene1_end >= gene2_end):
        return True
    else:
        return False
    

def binary_search(genes, pos, iteration=0):
    """
    Takes three args: a list of genes as tuples with (name, start position, end position),
    a integer that is the position of interest, and an optional iteration number (leave as
    default if running from main). Recursively finds the closest gene in the list of genes
    to the position given, and returns the gene name, start pos, end pos, distance from
    the nearest ent of the gen to the position of interest, and the number of iterations
    it took.
    """
    iteration += 1
    check_idx = len(genes) // 2
    gene_name, gene_start, gene_end = genes[check_idx]
    dist = distance_to_gene(genes[check_idx], pos)
    
    # Woohoo recursion is fun
    if dist == 0:
            return (gene_name, gene_start, gene_end, dist, iteration)
        
    elif pos < gene_start:
        if check_idx == 0:
            return (gene_name, gene_start, gene_end, dist, iteration)
        
        elif dist < distance_to_gene(genes[check_idx - 1], pos):
            return (gene_name, gene_start, gene_end, dist, iteration)
            
        else:
            return binary_search(genes[:check_idx], pos, iteration)
    
    # Same as previous elif except for reverse case
    elif pos > gene_end:
        if check_idx == (len(genes) - 1):
            return (gene_name, gene_start, gene_end, dist, iteration)
        
        elif dist < distance_to_gene(genes[check_idx + 1], pos):
            return (gene_name, gene_start, gene_end, dist, iteration)
            
        else:
            return binary_search(genes[check_idx:], pos, iteration)
    

def gtf_parser(gtf_filename, mut_chrom):
    """
    Takes name of gtf file and chromosome of interest, and parses gtf into a list of
    tuples (gene_name, start, end) that includes only protein coding genes on the
    chromosome of interest. Excludes duplicates or genes enclosed by other genes.
    """
    f = open(gtf_filename)
    prev_gene = ('', 0, 0)

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
        # I modified this so it won't add 'duplicate' genes (genes that start and end at the same postion) or enclosed genes (genes enclosed by other genes)
        if (fields[0] == mut_chrom) and (fields[2] == "gene") and ('gene_biotype "protein_coding"' in line):
            subfields = fields[-1].split(';') # Take last field
            for field in subfields:
                if "gene_name" in field:
                    # Find 'gene name' and then take the string after
                    gene_name = field.split()[1]
        
                    current_gene = (gene_name, start, end)
                    if not encloses(prev_gene, current_gene):
                        # Checks to make sure this gene isn't enclosed by previous one
                        genes.append(current_gene) # Builds list of tuples
                        prev_gene = current_gene
    
    return genes
    
    

def main():
    """
    Parses .gtf file into ordered list of gene names and positions and then uses
    that to find which gene is closest to a postion of interest on a particular
    chromosome.
    """
    # Read in arguments
    gtf_file = sys.argv[1]
    mut_chrom = sys.argv[2]
    mut_pos = int(sys.argv[3])

    genes = gtf_parser(gtf_file, mut_chrom)

    closest = binary_search(genes, mut_pos) # No iteration param so that default is used
    print("The nearest gene is {}, which runs from positions {} to {}. "
           "It is {} base pairs from your position. This took {} "
           "iterations.".format(*closest))

    
    
if __name__ == '__main__':
    main()