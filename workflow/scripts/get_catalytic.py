'''
Get the catalytic residues for every UniProt ID
'''

import pandas as pd
from structure_features.catalyticSites import parse_data
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
    

    parser = argparse.ArgumentParser(description=('Get the catalytic residues for '
                                                  'every UniProt ID.'))

    parser.add_argument('input', help='Input variants DataFrame', type=validate_file)
    
    parser.add_argument('output', help='Output file', type=Path)
    
    return parser.parse_args(args)


if __name__ == '__main__':
    
    args = parsing()

    # Get the number of job from the arguments.
    # job = int(sys.argv[1])
    # commands_per_job = int(sys.argv[2])

    # variants = pd.read_pickle('./11_variants.pkl')
    variants = pd.read_pickle(args.input)

    # # Create a new column for the residue positions
    # variants['Residue_position'] = (variants['Protein_position'].str.split('/')
    #                                 .str[0].astype(int))

    # # Create a new column for the mutated residue in 3-letter code
    # variants['Residue'] = (variants['Amino_acids'].str.split('/').str[0]
    #                        .map(IUPACData.protein_letters_1to3))

    # Get a list of unique UniProt IDs
    uids = list(variants['UniProt_IDs'])
    uids = [item for sublist in uids for item in sublist]
    uids = [uid.split('-')[0] for uid in uids]
    uids = list(dict.fromkeys(uids)) # Remove duplicates while preserving order

    # Make a list of indices to run in this job
    # indices = np.arange(0, len(uids), commands_per_job)
    # start = indices[job]
    # end = start + commands_per_job

    # Get the catalytic residues for the UniProt IDs in this job
    catalytic_residues = parse_data(uids)

    # Create a dataframe with the catalytic residues
    catalytic_residues = pd.DataFrame(catalytic_residues).drop_duplicates()

    # Save to pickle
    # catalytic_residues.to_pickle(f'./catalytic_residues/12_catalytic_residues_{job}.pkl')
    catalytic_residues.to_pickle(args.output)

    print('Done.')