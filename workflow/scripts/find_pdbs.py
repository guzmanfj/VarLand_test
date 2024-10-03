import argparse
import pandas as pd
from pathlib import Path
import biskit as b
import re
import numpy as np
from typing import List


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
    
    def convert_range(s:str) -> List[int]:
        """
        Convert the range of indexes from a string to integers
        """
        print(s)
        l = [int(x) for x in s.split('-')]
        return l

    def validate_dir(d:str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            d.mkdir(parents=True)
            
        return d
    

    parser = argparse.ArgumentParser(description=('Obtain the PDB files for the '
                                                  'variants, if existing.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('variants_ranges', help='String with the start and end '
                        'of the range of variants to process, separated by a dash',
                        type=convert_range)
    
    parser.add_argument('alphafold_database', help='Path to the AlphaFold database',
                        type=validate_file)
    
    parser.add_argument('alphafold_mane_database', help='Path to the AlphaFold MANE database',
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


def check_sequence(pdbfile:Path, seq_length:int, aa:str, aa_pos:int) -> bool:
    """
    Load the pdbfile and check if the sequence length and amino acid match
    """
    m = b.PDBModel(str(pdbfile))
    seq = m.sequence()
    if not len(seq) == seq_length:
        return False

    if not seq[aa_pos-1] == aa:
        return False

    return True


def get_model_number(aapos:float, models:list):
    """
    Determine which model to get for the amino acid position
    """
    # If less than 800, just read the first file
    if aapos <= 800:
        return 0

    # After 800, read files in 200 increments
    else:
        excess = aapos-800
        quotient = int(np.ceil(excess/200))

    # If the amino acid position is higher than the 200-residue limit, it would
    # normally go into the next file. But if there is no other file, return the last file
    if len(models) <= quotient:
        return len(models)-1

    # The amino acid falls into the 200-residue limit of one of the files
    elif len(models) > quotient:
        return quotient


def find_pdbs(variants:pd.DataFrame, af_database:Path) -> dict[int, Path]:
    """
    Iterate over the variants, and for each variant, iterate over the UniProt IDs
    until we find one that has a structure in the `af_database` directory

    Args:
        variants (pd.DataFrame): DataFrame with the variants
        af_database (Path): Path to the AlphaFold database

    Returns:
        dict[int, Path]: Dictionary with the variant index as key, and the path
            to the PDB file as value
    """
    
    pdbs = {}
    for i, row in variants.iterrows():
        pdbs[i] = np.nan
        for uidp in row['UniProt_IDs']:
            uid = uidp.split('-')[0]
            pdb_files = list(af_database.glob(f'AF-{uid}-F*.pdb'))
            if len(pdb_files) == 0:
                continue
            elif len(pdb_files) == 1:
                aa_pos, seq_length = row['Protein_position'].split('/')
                aa_pos = int(aa_pos)
                seq_length = int(seq_length)
                aa = row['Amino_acids'].split('/')[0]
                if check_sequence(pdb_files[0], seq_length, aa, aa_pos):
                    pdbs[i] = pdb_files[0]
                    break
            else:
                # There is more than 1 pdb file for this UniProt ID
                # select the correct one based on the amino acid position

                # First sort the PDB files by the model number
                pdb_files = sorted(pdb_files,
                            key=lambda m: int(re.search(r'-F(\d+)-', m.stem).group(1)))

                aa_pos, seq_length = row['Protein_position'].split('/')
                aa_pos = int(aa_pos)
                seq_length = int(seq_length)
                aa = row['Amino_acids'].split('/')[0]

                # Get the index of the model to use
                mindex = get_model_number(aa_pos, pdb_files)

                # Calculate the expected fragment length
                start = mindex*200
                aa_pos = aa_pos - start
                if (seq_length - start) > 1400:
                    fragment_length = 1400
                else:
                    fragment_length = seq_length - start

                # Check if the sequence matches
                if check_sequence(pdb_files[mindex], fragment_length, aa, aa_pos):
                    pdbs[i] = pdb_files[mindex]
                    break
    
    return pdbs


if __name__ == '__main__':
    
    args = parsing()

    all_variants = pd.read_pickle(args.input)

    start = args.variants_ranges[0]
    end = args.variants_ranges[1]
    variants = all_variants.iloc[start:end].copy()
    del all_variants
    
    ## Search in the MANE database
    af_database = args.alphafold_mane_database
    pdbs = find_pdbs(variants, af_database)

    # Add the PDBs as a column to the variants df
    variants['PDB'] = pd.Series(pdbs)

    # Obtain the missing variants
    missing_variants = variants[variants.PDB.isna()]

    if len(missing_variants) > 0:
        ## Search in the normal AlphaFold database
        af_database = args.alphafold_database
        missing_pdbs = find_pdbs(missing_variants, af_database)

        # Add column with the missing PDBs, and merge with the first one
        variants['PDB_af'] = pd.Series(missing_pdbs)
        variants['PDB_path'] =  variants['PDB'].combine_first(variants['PDB_af'])
        variants = variants.drop(columns=['PDB', 'PDB_af'])

    else:
        variants['PDB_path'] = variants['PDB']
        variants = variants.drop(columns=['PDB'])

    # Save to pickle
    variants.to_pickle(args.output)

    print('Done.')

