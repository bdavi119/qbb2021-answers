def main():
    """
    Reads in .sam file in same directory as this program. Performs various
    summary calculations and writes results to a file named 'day1-evening.out'.
    """
    alignment_count = 0
    perfect_count = 0
    total_mapq = 0
    in_range = 0
    uniq_rds = set()
    # duplicates = 0
    # prev_pos = None
    # prev_cigars = set()
    
    f = open('SRR072893.sam', 'r')
    for line in f:
        if line.startswith('@'):
            continue # If this is a header line
        
        # Part 1
        alignment_count += 1
        
        # Used in future parts
        read = line.strip().split('\t')
    
        # Part 2
        # Makes sure that the '40M' check later works
        if len(read[9]) != 40:
            print("Error, there's a read that's not 40 bases long")

        if read[5] == '40M':
            perfect_count += 1
    
        # Part 3
        # Makes sure there's not an empty mapq entry
        if read[4] == '':
            print("Error, at least one MAPQ value missing")
            
        total_mapq += int(read[4])
    
        # Part 4
        pos = int(read[3])
        if (read[2] == '2L') & (pos >= 10000) & (pos <= 20000):
            in_range += 1
    
        # Part 5
        uniq_rds.add((pos, read[5]))
        
        """
        Other method for Part 5. It appears that while positions are in order,
        CIGAR strings are not necessarily ordered. This catches CIGAR strings
        of the same position that are not next to each other. To use, also
        uncomment three lines near beginning and write statement at end

        if (pos == prev_pos) & (read[5] in prev_cigars):
            duplicates += 1
        elif (pos == prev_pos) & (read[5] not in prev_cigars):
            prev_cigars.add(read[5])
        else:
            prev_cigars = set([read[5]])
        prev_pos = pos
        """

    f.close()
    
    # Final calulculations for Part 3 and 5
    avg_mapq = total_mapq / alignment_count
    PCR_dupes = alignment_count - len(uniq_rds)

    # Write to file
    out = open('day1-evening.out', 'w')
    out.write('Answers\n')
    out.write('Part 1: {} total alignments\n'.format(alignment_count))
    out.write('Part 2: {} perfect matching alignments\n'.format(perfect_count))
    out.write('Part 3: {} is the average MAPQ score\n'.format(avg_mapq))
    out.write('Part 4: {} reads on 2L between 10000 and 20000 inclusive\n'.format(in_range))
    out.write('Part 5: {} potential PCR duplicates\n'.format(PCR_dupes))
    # out.write('Part 5 alt: {} potential PCR duplicates\n'.format(duplicates))
    out.close()

if __name__ == '__main__':
    main()