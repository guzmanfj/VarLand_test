import argparse
import pandas as pd
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
    
    def validate_dir(d:str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            raise ValueError("The specified directory doesn't exist.")
            
        return d
    

    parser = argparse.ArgumentParser(description=('Takes a directory containing '
                'subdirectories with the .vep files for each chromosome. Reads '
                'all .vep files and concatenates them into a single dataframe.'))

    parser.add_argument('input', help='Input directory', type=validate_dir)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


if __name__=='__main__':
    
    args = parsing()

    # Define the names of the chromosome directories
    chr_dirs = [args.input / f'chr{i}_split' for i in range(1, 23)]
    chr_dirs = chr_dirs + [args.input / 'chrX_split', args.input / 'chrY_split']

    # Read all files ending with .vep for each chromosome directory
    chrom_dfs = []
    for d in chr_dirs:
        print(f'Reading {d}')
        # Data starts at line 38 (1-based)
        df = pd.concat([pd.read_csv(f, sep='\t', header=35) for f in d.glob('*.vep')],
                    ignore_index=True)
        chrom_dfs.append(df)

    # Concatenate all dataframes into one
    df = pd.concat(chrom_dfs, ignore_index=True)

    df.to_pickle(args.output)
