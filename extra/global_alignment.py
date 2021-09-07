#!/usr/bin/env python3
"""
Brendon Davis

This program contains a small suite of functions and a main method
for performing a global alignment on two DNA sequences with the
HOXD70 scoring matrix.

Running this program requires at least one additional argument:
the file path for the FASTA file with the two sequences to be
aligned. Optional second argument is the word 'all', which
will make the program produce all optimal alignments but may
drastically increase runtime for long sequences.

test.fasta and test_mult_optimal.fasta have been provided as well
as simple SIGMA and GAP values below for testing. test.fasta contains
the two sequences from the example in the assignment and
test_mult_optimal.fasta contain two sequences that have multiple optimal
alignments when using the simple SIGMA and GAP values (use this to test
'all' optional command line argument). The other provided file,
to_align.fasta, contains the longer sequences from the assignment to align.
"""

import numpy as np
from fasta_reader import FASTAReader
import sys

base_order = ['A', 'C', 'G', 'T'] # Needed for interpreting sigma correctly

# HOXD70 scoring matrix
SIGMA = np.array([[91, -114, -31, -123],
                 [-114, 100, -125, -31],
                 [-31, -125, 100, -114],
                 [123, -31, -114, 91]])

GAP = 300 # gap penalty

"""
# These sigma and gap values are for testing with
# the example sequences in the assignment
SIGMA = np.array([[1,-1,-1,-1],
                  [-1,1,-1,-1],
                  [-1,-1,1,-1],
                  [-1,-1,-1,1]])

GAP = 1
"""


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
            d = F[i-1, j-1] + SIGMA[base_order.index(base1), base_order.index(base2)]
                                # base_order is used to properly index into the correct entry of SIGMA
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


def reverse_traverse(hist, all_alignments):
    """
    hist = Numpy array with character(s) in each entry
    corresponding to entries in N-W array that describe
    how that N-W entry was produced
    all_alignments = Boolean describing whether all
    optimal alignments are desired or just one
    
    Auxiliary function to aligner() that uses recursion to
    produce a Strng with characters describing the optimal
    path through the N-W array. Returns a list with one String
    all_alignments is set to False or potentially multiple
    Strings if all_alignments is set to True.
    
    Note 1: Has a long runtime for mildly long (>100 bp's) sequences
    if all_alignments is set to True. Not suggested for long sequences.
    
    Note 2: This function could probably be more easily solved
    with a while loop if we didn't care about multiple optimal
    alignments. However, recursion allows us to more easily
    retrieve multiple optimal alignments if desired.
    """
    if hist.shape == (1, 1):
        # If we're done reversing through the array
        return ['']
    
    current_dir = hist[-1, -1]
    if (not all_alignments) and (len(current_dir) > 1):
        # If we only want one alignment instead of all optimal alignments
        if 'd' in current_dir:
            current_dir = 'd'
        else:
            current_dir = 'v'
        # This prioritzes d over v over h, which is kind of arbitrary,
        # but by prioritzing d it should generally shorten the alignment
    
    directions = []
    for char in current_dir:
        # Cuts down hist array into new array for next recursive call
        if char == 'v':
            new_hist = hist[:-1, :]
        elif char == 'h':
            new_hist = hist[:, :-1]
        elif char == 'd':
            new_hist = hist[:-1, :-1]
        
        # Recursion is fun, don't @ me
        dir_to_here = reverse_traverse(new_hist, all_alignments)
        for string in dir_to_here:
            string = string + char # Builds final string in forward order 
                                   # (as directions from top left of array)
            directions.append(string) # Using a list allows for multiple optimal alignments
   
    return directions
        


def aligner(seq1, seq2, hist, all_alignments):
    """
    seq1 = String with DNA sequence #1
    seq2 = String with DNA sequence #2
    hist = Numpy array with character(s) in each entry
    corresponding to entries in N-W array that describe
    how that N-W entry was produced
    all_alignments = Boolean describing whether all
    optimal alignments are desired or just one
    
    Returns integer with number of optimal alignments that
    were performed. Will always be 1 if all_alignments is
    set to False. Also prints out the alignment(s) in lines
    of 100 bp's to the terminal.
    """
    a_seq1 = ''
    a_seq2 = ''
    same = '' # Typographical representation of matching bases in alignment
              # * represents mismatch and | represents match
    directions = reverse_traverse(hist, all_alignments)
    
    for i in range(len(directions)):
        print('Alignment {}:'.format(i+1))
        
        loop_num = 0
        current_seq1 = 0 # current_seq vars (representing indices) allow us to move seperately though seq1 and seq2
        current_seq2 = 0
        for char in directions[i]:
            loop_num += 1
            
            if char == 'h':
                a_seq1 = a_seq1 + '-'
                a_seq2 = a_seq2 + seq2[current_seq2]
                same = same + '*'
                current_seq2 += 1
                
            if char == 'v':
                a_seq1 = a_seq1 + seq1[current_seq1]
                a_seq2 = a_seq2 + '-'
                same = same + '*'
                current_seq1 += 1
                
            if char == 'd':
                a_seq1 = a_seq1 + seq1[current_seq1]
                a_seq2 = a_seq2 + seq2[current_seq2]
                if seq1[current_seq1] == seq2[current_seq2]:
                    same = same + '|'
                else:
                    same = same + '*'
                current_seq1 += 1
                current_seq2 += 1
            
            if len(same) % 100 == 0:
                # Once strings get long enough we print and move to new line
                print(a_seq1)
                print(same + str(loop_num))
                print(a_seq2 + '\n')
                a_seq1, a_seq2, same = ('', '', '')
            elif loop_num == len(directions[i]):
                # If we're at the end of the alignment
                print(a_seq1)
                print(same + str(loop_num))
                print(a_seq2 + '\n')
                a_seq1, a_seq2, same = ('', '', '')
        
        if i+1 != len(directions):
            # Skipping line to seperate multiple alignments
            print('')
            
    return(len(directions)) # Not necessary for this assignment but could be useful as a return
    


def main():
    """
    Reads in FASTA file (argument 1 from command line) and
    prints out one optimal alignment (all optimal alignments
    if 'all' is given as second argument) for the two sequences.
    """
    filename = sys.argv[1]
    all_alignments = False # default if no second additional argument given

    if len(sys.argv) > 2: 
        if sys.argv[2] == 'all':
            all_alignments = True

    seqs = FASTAReader(open(filename)) # FASTAReader() comes from bootcamp
    assert len(seqs) == 2 # Only two sequences allowed
    seq1 = seqs[0][1]
    seq2 = seqs[1][1]
    
    # arr (the N-W array) isn't needed for this assignment
    # but feels right for needleman_array() to return it
    arr, hist, similarity = needleman_array(seq1, seq2)

    print("The sequence similarity is {}.\n".format(similarity))
    aligner(seq1, seq2, hist, all_alignments) 
    # aligner() automatically prints alignments so we don't need to do anything else from main()


if __name__ == '__main__':
    main()