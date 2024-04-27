from Bio.Data import IUPACData
import pandas as pd
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
    

    parser = argparse.ArgumentParser(description=('Preprocess the variants pickle file.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('dssp_alphafold_database', help='Path to the AlphaFold DSSP database',
                        type=validate_file)
    
    parser.add_argument('dssp_mane_database', help='Path to the AlphaFold MANE DSSP database',
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)



if __name__ == '__main__':

    args = parsing()

    ## Pre-process variants

    # variants = pd.read_pickle('./11_variants.pkl')
    variants = pd.read_pickle(args.input)

    # Remove variants without PDB
    variants = variants[variants.PDB_path.notna()]
    # variants = variants.reset_index(drop=True)

    # Create a new column for the residue positions
    variants['Residue_position'] = (variants['Protein_position'].str.split('/')
                                    .str[0].astype(int))

    # Create a new column for the mutated residue in 3-letter code
    variants['Residue'] = (variants['Amino_acids'].str.split('/').str[0]
                        .map(IUPACData.protein_letters_1to3))

    # Make a column with the DSSP path
    dssp_paths = {}

    # iterate over the variants and choose the correct DSSP path depending on the PDB path
    for i, var in variants.iterrows():
        
        # get the PDB path
        pdb_path = var['PDB_path']
        
        # get the DSSP path
        if 'mane' in str(pdb_path):
            dpath = args.dssp_mane_database / (pdb_path.name+'.dssp')
            dssp_paths[i] = dpath 
        else:
            dpath = args.dssp_alphafold_database / (pdb_path.name+'.dssp')
            dssp_paths[i] = dpath

    # Add the DSSP path to the dataframe
    variants['DSSP_path'] = pd.Series(dssp_paths)

    # Save to pickle
    # variants.to_pickle('./14_variants.pkl')
    variants.to_pickle(args.output)

    print('Done! Variants pre-processed.')