"""
Script to create a dictionary mapping UniProt IDs to their corresponding PDB files
"""

import argparse
from pathlib import Path
from multiprocessing import Pool
from typing import Dict, List
import pickle

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
    

    parser = argparse.ArgumentParser(description=('Given a directory with PDB files, '
                    'create a mapping between UniProt ID and file name.'))

    parser.add_argument('input_dir', help=('Directory with the PDB files.'),
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


def create_uniprot_pdb_map(af_database: Path) -> Dict[str, List[Path]]:
    """
    Create a dictionary mapping UniProt IDs to their corresponding PDB files
    """
    pdb_map = {}
    for pdb_file in af_database.glob('AF-*-F*.pdb'):
        uid = pdb_file.stem.split('-')[1]
        if uid not in pdb_map:
            pdb_map[uid] = []
        pdb_map[uid].append(pdb_file)
    return pdb_map


if __name__=='__main__':
    args = parsing()
    
    pdb_map = create_uniprot_pdb_map(args.input_dir)
    
    with open(args.output, 'wb') as f:
        pickle.dump(pdb_map, f)
    
    print(f'UniProt IDs extracted from {len(pdb_map)} PDB files')