import pandas as pd
import numpy as np
import pickle
from pathlib import Path
import argparse


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
                ' the variants features into a single DataFrame.'))

    parser.add_argument('input', help='Input directory with the DataFrames',
                        type=validate_file)
    
    parser.add_argument('foldx_energies', help='Input directory with the foldx '
                        'energies dictionaries', type=validate_file)
    
    parser.add_argument('output', help='Output file', type=Path)
    
    return parser.parse_args(args)


if __name__ == '__main__':
    
    args = parsing()

    # njobs = int(sys.argv[1])

    # Read all variants dataframes
    variants = []
    for file in args.input.glob('15_variants_*.pkl'):
        ivariants = pd.read_pickle(file)
        variants.append(ivariants)

    variants = pd.concat(variants)

    # Clean secondary structure column
    variants.secondary_structure.replace(' ', np.nan, inplace=True)

    # Remove from the contacts residues less than 6 residues away in the sequence:
    # Create a dictionary with non-adjacent contacts for every variant.
    # Keep residues that are at least 6 amino acids away in the sequence
    noadj_contacts = {}
    for i, var in variants.iterrows():
        res_pos = var['Residue_position']
        contacts = var['intra_contacts']
        select = np.abs(res_pos - contacts) >= 6
        noadj_contacts[i] = contacts[select]

    variants['intra_contacts'] = pd.Series(noadj_contacts)
    variants['Conservation'] = variants['Conservation'].astype(float)
    variants['BLOSUM62'] = variants['BLOSUM62'].astype(int)


    ## Add foldx energies

    # Read all the dictionaries with the foldx results and merge them
    foldx_energies = {}
    for file in args.foldx_energies.glob('17_foldx_energies_*.pkl'):
        with open(file, 'rb') as f:
            foldx_energies.update(pickle.load(f))

    # Create a dataframe with the foldx results
    foldx = pd.DataFrame(foldx_energies).T

    # Merge with the variants DataFrame by the index
    variants = variants.merge(foldx, left_index=True, right_index=True)

    print(variants.shape)

    # See if there are duplicated rows
    select_columns = [
        '#CHROM', 'POS', 'REF', 'ALT', 'RefSeq', 'PDB_path', 'Amino_acids',
        'Residue_position', 'pLDDT'
    ]

    # Remove the duplicated rows in the full dataframe based on the selected columns
    variants = variants.drop_duplicates(subset=select_columns).reset_index(drop=True)

    # Save to pickle
    # variants.to_pickle('./AM_pathogenic.pkl')
    variants.to_pickle(args.output)

    print(variants.shape)
    print('Done.')
