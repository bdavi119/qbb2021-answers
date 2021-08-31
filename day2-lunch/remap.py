#!/usr/bin/env python3

import sys

def main():
    """
    Takes two arguments (file names) and optional third argument. First file
    must be gene_id to protein_id mapping file. Second file must be ctab with
    gene_id's to be replaced. Third argument can be a word to use in replacement
    of the protein_id if the gene_id of the corresponding entry is not in the
    mapping file. Entry will be skipped if no third line given. Prints results
    that can be redirected to output file.
    """
    args = list(sys.argv)
    
    # First create dictionary with mapping file
    map_file = open(args[1])
    
    mapping = dict()
    for line in map_file:
        key, value = line.strip().split('\t')
        mapping[key] = value
    
    map_file.close()
    
    # Next go through ctab and change gene_id to protein_id
    ctab = open(args[2])
    header = True
    for gene in ctab:
        parsed_gene = gene.split('\t')
        
        # This edits the header line instead of looking for gene_id
        if header:
            parsed_gene[8] = 'protein_id'
            print('\t'.join(parsed_gene))
            header = False
            continue
        
        gene_id = parsed_gene[8]
        
        if gene_id in mapping:
            parsed_gene[8] = mapping[gene_id]
        elif len(args) == 3:
            # If no third arg given
            continue
        else:
            # If third arg given
            parsed_gene[8] = args[3]
        
        print('\t'.join(parsed_gene))
    
    ctab.close()
            

if __name__ == '__main__':
    main()