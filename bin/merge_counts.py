#!/usr/bin/env python3

import pandas as pd
import sys
import os

def main():
        # Nextflow passes the file paths as command-line arguments
        # sys.argv[1:] captures all files passed to the script
        files = sys.argv[1:]
        if not files:
            print("No count files provided to merge.")
            sys.exit(1)

        df_list = []

        for f in files:
            # featureCounts has a 1-line header, so we skip row 0
            # We only need 'Geneid' (Col 0) and the actual count (Col 6)
            df = pd.read_csv(f, sep='\t', skiprows=1, usecols=[0, 6])

            # Clean up the BAM column name to extract the sample name
            bam_col = df.columns[1]
            sample_name = os.path.basename(bam_col).replace('_sorted.bam', '')
            df.rename(columns={bam_col: sample_name}, inplace=True)

            df.set_index('Geneid', inplace=True)
            df_list.append(df)

        # Merge all dataframes side-by-side on the Geneid index
        merged_df = pd.concat(df_list, axis=1)

        # Save the output file in the current working directory
        merged_df.to_csv('raw_counts_matrix.csv')
        print("Successfully generated raw_counts_matrix.csv")

if __name__ == '__main__':
        main()
