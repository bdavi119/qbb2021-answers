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
    
    f = open('SRR072893.sam', 'r')
    for line in f:
        if not line.startswith('@'):
            alignment_count += 1 # For part 1
            read = line.strip().split('\t')
        
            # Part 2
            # Not required, but makes sure that the '40M' check later works
            if len(read[9]) != 40:
                print("Error, there's a read that's not 40 bases long")
 
            if read[5] == '40M':
                perfect_count += 1
        
            # Part 3
            # Not required 
            if read[4] == None:
                print("Error, at least one MAPQ value missing")
            total_mapq += int(read[4])
        
            # Part 4
            pos = int(read[3])
            if (read[2] == '2L') & (pos >= 10000) & (pos <= 20000):
                in_range += 1
        
            # Part 5
            uniq_rds.add((pos, read[5]))

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
    out.close()

if __name__ == '__main__':
    main()