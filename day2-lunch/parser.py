#!/usr/bin/env python3

import sys

def main():
    """
    Takes one file name as an argument. FIle should be "messy" flybase data.
    For each gene entry, prints gene_id and protein_id seperated by tabs.
    """
    args = list(sys.argv)   
    og_file = open(args[1])
    
    for line in og_file:
        if "DROME" not in line:
            # If this is one of the header lines
            continue
        
        prot_id = line[40:46] # Counted out the characters for these ID's
        gene_id = line[52:63]
        print(gene_id + '\t' + prot_id)
    
    og_file.close()
            

if __name__ == '__main__':
    main()