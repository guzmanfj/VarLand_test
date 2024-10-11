'''
Script to write the variants to files for dbnsfp to process them.
'''

import pandas as pd
import argparse
from pathlib import Path

def parsing(args: list=None) -> argparse.Namespace:
    """
    Creates the argument parser instance and applies it to the command line
    input

    Args:
        args (list, optional): List of the arguments to be parsed (only to be
            used for testing). If none is provided, it is taken from sys.argv.
            Defaults to None.

    Returns:
        argparse.Namespace
    """
    
    def validate_file(d:str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            raise ValueError("The specified file or directory doesn't exist.")
            
        return d
    

    parser = argparse.ArgumentParser(description=('Extract the variants and locations'
            ' and write them in a file for dbNSFP.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('output', help='Output text file', type=Path)
    
    return parser.parse_args(args)


if __name__ == '__main__':
    
    args = parsing()

    variants = pd.read_pickle(args.input)

    if not variants.empty:

        # Write the variants to a file in the format
        # chr pos ref alt res_wt res_mut
        variants_formatted = variants[['#CHROM','POS','REF','ALT','Amino_acids']].copy()
        # Create two columns from the 'Amino_acids' column, splitting it by '/'
        variants_formatted[['res_wt','res_mut']] = variants_formatted['Amino_acids'].str.split('/', expand=True)

        variants_formatted.drop(columns=['Amino_acids'], inplace=True)

        # Remove the 'chr' prefix from the chromosome column
        variants_formatted['#CHROM'] = variants_formatted['#CHROM'].str.replace('chr','')

        # Save to csv without header and index
        variants_formatted.to_csv(args.output, index=False, header=False,
                                sep='\t')
    else:
        args.output.touch()

    print(f'Done.')
