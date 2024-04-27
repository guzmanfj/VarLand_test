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
    

    parser = argparse.ArgumentParser(description=('Combine all the pickles with'
                ' the catalytic residues into a single DataFrame.'))

    parser.add_argument('input', help='Input directory with the DataFrames',
                        type=validate_file)
    
    parser.add_argument('output', help='Output file', type=Path)
    
    return parser.parse_args(args)


if __name__ == '__main__':
    
    args = parsing()

    # njobs = int(sys.argv[1])

    # Read all the catalytic residues dataframes
    catalytic_residues = []
    for file in args.input.glob('12_catalytic_residues_*.pkl'):
        icatalytic_residues = pd.read_pickle(file)
        catalytic_residues.append(icatalytic_residues)

    catalytic_residues = (pd.concat(catalytic_residues).drop_duplicates()
                          .reset_index(drop=True))

    # Save to pickle
    # catalytic_residues.to_pickle('./13_catalytic_residues.pkl')
    catalytic_residues.to_pickle(args.output)