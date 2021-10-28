#!/usr/bin/env python3
"""
Brendon Davis

This program contains a small suite of functions and a main method
for performing a global alignment on two AA/DNA sequences.

Running this program requires at four additional argument:
    - The file path for the FASTA file with the two sequences to be aligned
    - The file path for the scoring matrix
    - An integer to be used as the gap penalty
    - A filepath to write the alignment output to
"""

import numpy as np
from fasta_reader import FASTAReader
import sys


def read_sigma(filename):
    """
    filename = String the with filepath to a txt file containing the
               the scoring matrix to be used. File is expected to be
               whitespace delimited with single-character residues listed
               along the first line and as the first element of each row,
               creating the "coordinate system" of the matrix.
    
    Returns the sigma array as a 2-D numpy array and the residue order
    (as it is written left-to-right or top-to-bottom) as a list.
    """
    f = open(filename, 'r')
    lines = f.readlines()
    residue_order = lines[0].strip().split()
    
    size = len(residue_order)
    sig = np.zeros((size, size), dtype=int)
    line_num = -1
    for line in lines[1:]:
        line_num += 1
        line = line.strip().split()
        sig[line_num, :] = line[1:]
    
    return sig, residue_order

def needleman_array(seq1, seq2):
    """
    seq1 = String with DNA sequence #1
    seq2 = String with DNA sequence #2
    
    Returns two numpy arrays and an integer. First array is
    the Needleman-Wunsh (N-W) array with sequence #1 on the rows and
    sequence #2 on the columns. Second array documents which
    entry(ies) produced the similarity score for the corresponding
    entry in the N-W array (v=vertical, h=horizontal, d=diagonal).
    Returned integer is the final similarity score for the two
    sequences.
    """
    m = len(seq1) + 1
    n = len(seq2) + 1
    
    F = np.zeros((m, n), dtype=int)
    
    # Note: history has to be dtype object instead of string
    # to allow for strings of different length
    history = np.empty((m, n), dtype=object)
    
    # First two for-loops define the top right and
    # left column of N-W array
    for i in range(1, m):
        F[i, 0] = -GAP * i
        history[i][0] = 'v'
    for j in range(1, n):
        F[0, j] = -GAP * j
        history[0][j] = 'h'
        
    # Fill in the rest of the N-W array and keep track
    # of what entry(ies) produce next entry in history array
    for i in range(1, m):
        for j in range(1, n):
            base1 = seq1[i-1]
            base2 = seq2[j-1]
            
            v = F[i-1, j] - GAP
            h = F[i, j-1] - GAP
            d = F[i-1, j-1] + SIGMA[ORDER.index(base1), ORDER.index(base2)]
                                # ORDER is used to properly index into the correct entry of SIGMA
            F[i,j] = max(v, h, d)
                
            origin = ''
            if (v >= h) and (v >= d):
                origin += 'v'         
            if (h >= v) and (h >= d):
                origin = origin + 'h'
            if (d >= v) and (d >= h):
                origin += 'd'
            # origin can have more than one character if there is
            # a tie. This means there is more than one optimal alignment
            
            history[i, j] = origin
    
    return (F, history, F[len(seq1), len(seq2)])


def reverse_traverse(hist):
    """
    hist = Numpy array with character(s) in each entry
    corresponding to entries in N-W array that describe
    how that N-W entry was produced
    
    Auxiliary function to aligner() that uses recursion to
    produce a Strng with characters describing the optimal
    path through the N-W array. Returns a list with one String.
    
    Note 2: This function could probably be more easily solved
    with a while loop if we didn't care about multiple optimal
    alignments. However, recursion allows us to more easily
    retrieve multiple optimal alignments if desired.
    -----------------------------------
    """
    directions = ''
    while hist.shape != (1, 1):
        char = hist[-1, -1]
        if len(char) > 1:
            if 'd' in char:
                char = 'd'
            else:
                char = 'v'
            # This prioritzes d over v over h, which is kind of arbitrary,
            # but by prioritzing d it should generally shorten the alignment
        
        if char == 'v':
            hist = hist[:-1, :]
        elif char == 'h':
            hist = hist[:, :-1]
        elif char == 'd':
            hist = hist[:-1, :-1]
        
        directions = char + directions # Builds final string in forward order 
                               # (as directions from top left of array)
    return (directions)
        


def aligner(seq1, seq2, hist, similarity, to_write):
    """
    seq1 = String with DNA sequence #1
    seq2 = String with DNA sequence #2
    hist = Numpy array with character(s) in each entry
    corresponding to entries in N-W array that describe
    how that N-W entry was produced
    
    Returns integer with number of optimal alignments that
    were performed. Also prints out the alignment in lines
    of 100 bp's to the terminal.
    """
    a_seq1 = ''
    a_seq2 = ''
    same = '' # Typographical representation of matching bases in alignment
              # * represents mismatch and | represents match
    directions = reverse_traverse(hist)
        
    loop_num = 0
    current_seq1, current_seq2 = (0, 0) # current_seq vars (representing indices) allow us to move seperately though seq1 and seq2
    gap1_count, gap2_count, mismatch_count = (0, 0, 0)
    
    f = open(to_write, 'w')
    for char in directions:
        loop_num += 1
        
        if char == 'h':
            a_seq1 = a_seq1 + '-'
            a_seq2 = a_seq2 + seq2[current_seq2]
            same = same + '*'
            current_seq2 += 1
            gap1_count += 1
            
        if char == 'v':
            a_seq1 = a_seq1 + seq1[current_seq1]
            a_seq2 = a_seq2 + '-'
            same = same + '*'
            current_seq1 += 1
            gap2_count += 1
            
        if char == 'd':
            a_seq1 = a_seq1 + seq1[current_seq1]
            a_seq2 = a_seq2 + seq2[current_seq2]
            if seq1[current_seq1] == seq2[current_seq2]:
                same = same + '|'
            else:
                same = same + '*'
                mismatch_count += 1
            current_seq1 += 1
            current_seq2 += 1
        
        if len(same) % 100 == 0:
            # Once strings get long enough we print and move to new line
            f.write(a_seq1 + '\n')
            f.write(same + str(loop_num) + '\n')
            f.write(a_seq2 + '\n\n')
            a_seq1, a_seq2, same = ('', '', '')
        elif loop_num == len(directions):
            # If we're at the end of the alignment
            f.write(a_seq1 + '\n')
            f.write(same + str(loop_num) + '\n')
            f.write(a_seq2 + '\n')
            a_seq1, a_seq2, same = ('', '', '')
        
    print('The sequence similarity is {}.'.format(similarity))
    print('Sequence 1 has {} gaps.'.format(gap1_count))
    print('Sequence 2 has {} gaps.'.format(gap2_count))   
    print('There are {} direct mismatches.'.format(mismatch_count))

    return(True)
    
    
def main():
    """
    Reads in FASTA file (argument 1 from command line) and
    writes alignment to specified file (argument 4).
    """
    filename = sys.argv[1]
    to_write = sys.argv[4]

    seqs = FASTAReader(open(filename)) # FASTAReader() comes from bootcamp
    assert len(seqs) == 2 # Only two sequences allowed
    seq1 = seqs[0][1]
    seq2 = seqs[1][1]
    
    # arr (the N-W array) isn't needed for hw assignment
    # but feels right for needleman_array() to return it
    arr, hist, similarity = needleman_array(seq1, seq2)
    print(arr)

    aligner(seq1, seq2, hist, similarity, to_write) 
    # aligner() automatically writes alignments so we don't need to do anything else from main()


if __name__ == '__main__':
    SIGMA, ORDER = read_sigma(sys.argv[2])
    GAP = int(sys.argv[3])
    main()