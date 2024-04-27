'''
Select variants that contain MANE transcripts
'''

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
    

    parser = argparse.ArgumentParser(description=('Filter variants with MANE '
                'transcripts. Saves the resulting DataFrame into a pickle.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


if __name__=='__main__':

    args = parsing()

    # variants = pd.read_pickle('./03_missense.pkl')
    variants = pd.read_pickle(args.input)

    select = np.logical_or(variants.MANE_SELECT.notna(),variants.MANE_PLUS_CLINICAL.notna())

    variants = variants[select]

    variants = variants.reset_index(drop=True)

    # variants.to_pickle('./04_mane.pkl')
    variants.to_pickle(args.output)

    print('Done.')

