import argparse
import pandas as pd
from pathlib import Path
import re
import numpy as np
from typing import List, Dict
from multiprocessing import Pool
from functools import partial
import mmap

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
    

    parser = argparse.ArgumentParser(description=('Obtain the PDB files for the '
                                                  'variants, if existing.'))

    parser.add_argument('input', help=('Input pickle file.'), type=validate_file)
    
    parser.add_argument('alphafold_database', help='Path to the AlphaFold database',
                        type=validate_file)
    
    parser.add_argument('alphafold_mane_database', help='Path to the AlphaFold MANE database',
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    parser.add_argument('--cpus', help='Number of CPUs to use', type=int, default=4)
    
    return parser.parse_args(args)


def check_sequence(seq: str, seq_length: int, aa: str, aa_pos: int) -> bool:
    """
    Check sequence using memory-mapped file for efficiency
    """
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


def process_variant(row: pd.Series, uniprot_pdb_map: Dict[str, List[Path]],
                    af_sequences:Dict[str, str]) -> Path:
    """
    Process a single variant
    """
    for uidp in row['UniProt_IDs']:
        uid = uidp.split('-')[0]
        if uid not in uniprot_pdb_map:
            continue
        
        pdb_files = uniprot_pdb_map[uid]
        aa_pos, seq_length = map(int, row['Protein_position'].split('/'))
        aa = row['Amino_acids'].split('/')[0]
        
        if len(pdb_files) == 1:
            seq = af_sequences[pdb_files[0].name]
            if check_sequence(seq, seq_length, aa, aa_pos):
                return pdb_files[0]
        else:
            # There is more than 1 pdb file for this UniProt ID
            # select the correct one based on the amino acid position

            # First sort the PDB files by the model number
            pdb_files = sorted(pdb_files,
                               key=lambda m: int(re.search(r'-F(\d+)-', m.stem).group(1)))
            
            # Get the index of the model to use
            mindex = get_model_number(aa_pos, pdb_files)
            
            # Calculate the expected fragment length
            start = mindex * 200
            aa_pos = aa_pos - start
            fragment_length = min(1400, seq_length - start)
            
            # Check if the sequence matches
            seq = af_sequences[pdb_files[mindex].name]
            if check_sequence(seq, fragment_length, aa, aa_pos):
                return pdb_files[mindex]
    
    return np.nan


def find_pdbs_parallel(variants: pd.DataFrame, uniprot_pdb_map: Dict[str, List[Path]],
                       af_sequences: Dict[str, str], num_processes: int) -> Dict[int, Path]:
    """
    Find PDBs using parallel processing
    """
    with Pool(num_processes) as pool:
        results = pool.map(partial(process_variant, uniprot_pdb_map=uniprot_pdb_map,
                                   af_sequences=af_sequences), 
                           [row for _, row in variants.iterrows()])
    return {i: path for i, path in enumerate(results)}

if __name__ == '__main__':
    args = parsing()
    variants = pd.read_pickle(args.input)
    
    af_sequences = pd.read_pickle(args.alphafold_database/'alphafold_human_v4_sequences.pkl')
    af_mane_sequences = pd.read_pickle(args.alphafold_mane_database/'alphafold_mane_sequences.pkl')
    af_uidmap = pd.read_pickle(args.alphafold_database/'alphafold_human_v4_uidmap.pkl')
    af_mane_uidmap = pd.read_pickle(args.alphafold_mane_database/'alphafold_mane_uidmap.pkl')
    
    # Search in the MANE database
    pdbs = find_pdbs_parallel(variants, af_mane_uidmap,
                              af_mane_sequences, args.cpus)
    variants['PDB'] = pd.Series(pdbs)
    
    # Process missing variants in the normal AlphaFold database
    missing_variants = variants[variants.PDB.isna()]
    if len(missing_variants) > 0:
        missing_pdbs = find_pdbs_parallel(missing_variants, af_uidmap,
                                          af_sequences, args.cpus)
        variants['PDB_af'] = pd.Series(missing_pdbs)
        variants['PDB_path'] = variants['PDB'].combine_first(variants['PDB_af'])
        variants = variants.drop(columns=['PDB', 'PDB_af'])
    else:
        variants['PDB_path'] = variants['PDB']
        variants = variants.drop(columns=['PDB'])
    
    variants.to_pickle(args.output)
    print('Done.')