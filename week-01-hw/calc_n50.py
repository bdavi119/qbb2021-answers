#!/usr/bin/env python3

import pandas as pd
import sys

def main():
    df = pd.read_csv(sys.argv[1], sep='\t', header=None)
    genome_length = int(sys.argv[2])
    df.sort_values(by=1, axis=0)
    
    current_length = 0
    i = -1
    while (current_length / genome_length) < 0.5:
        i += 1
        current_length += df.loc[i, 1]
    
    print(df.loc[i, 1])


if __name__ == '__main__':
    main()