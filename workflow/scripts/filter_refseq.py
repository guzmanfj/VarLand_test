import argparse
import pandas as pd
import numpy as np
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
            raise ValueError("The specified file doesn't exist.")
            
        return d
    

    parser = argparse.ArgumentParser(description=('Filter variants with the '
                'Refseq ID equal to the MANE ID. Saves the resulting DataFrame '
                'into a pickle.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


if __name__=='__main__':

    args = parsing()

    # variants = pd.read_pickle('./04_mane.pkl')
    variants = pd.read_pickle(args.input)

    print(variants.shape)

    # Split the `RefSeq` column and explode to create a row for each transcript
    variants['RefSeq_list'] = variants.RefSeq.str.split('&')
    variants = variants.explode('RefSeq_list')
    variants = variants.reset_index(drop=True)

    # Select the variants where the RefSeq ID is in either the MANE_SELECT or
    # MANE_PLUS_CLINICAL columns
    select_mane = variants['RefSeq_list'] == variants['MANE_SELECT']
    select_mane_plus = variants['RefSeq_list'] == variants['MANE_PLUS_CLINICAL']
    select = np.logical_or(select_mane, select_mane_plus)

    variants = variants[select]

    # Replace `RefSeq` column with column with unique RefSeq ID
    variants['RefSeq'] = variants['RefSeq_list']
    variants = variants.drop(columns=['RefSeq_list'])
    variants = variants.reset_index(drop=True)

    print(variants.shape)

    # variants.to_pickle('./05_variants.pkl')
    variants.to_pickle(args.output)

    print('Done.')

