#!/usr/bin/env python3
"""
Running this file requires exactly 3 additional arguments given from command line.
"""
import sys
from fasta_reader import FASTAReader

def get_target_kmers(target, k):
    """
    Takes parsed fasta file and integer k and builds a dictionary of all k-mers
    that appear in the file's sequences, and these k-mers are mapped to a list of
    tuples where each tuple is a sequence name and location where that k-mer was found.
    """
    
    kmer_locs = {}
    for name, seq in target:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:(i+k)]
            kmer_locs.setdefault(kmer, [])
            kmer_locs[kmer].append((name, i+1)) # Start counting bp's at 1
    
    return kmer_locs


def main():
    """
    Takes three arguments from command line: a target fasta file, a query fasta file, 
    and a number to define how long of kmers you want. Prints a list of all matches 
    between kmers in the query and target. Lines printed are in the format 
    ----target_sequence_name    target_start    query_start k-mer----
    """
    
    args = list(sys.argv)
    assert len(args) == 4
    target_file, query_file, k = args[1:]
    k = int(k)
    
    # Get database of kmers from target fasta file
    target = FASTAReader(open(target_file))
    database = get_target_kmers(target, k)
    
    # Get the query sequence from FASATReader
    query = FASTAReader(open(query_file))
    query_seq = query[0][1]
    
    for i in range(len(query_seq) - k + 1):
        kmer = query_seq[i:(i+k)]
        if kmer in database:
            for match in database[kmer]:
                print('{}\t{}\t{}\t{}'.format(match[0], match[1], i+1, kmer))
    

if __name__ == '__main__':
    main()
