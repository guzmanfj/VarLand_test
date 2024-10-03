"""
Script to extract sequences from PDB files
"""

import argparse
from pathlib import Path
from multiprocessing import Pool
import biskit as b
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
                    'extract the sequences from them and save them to a pickle file.'))

    parser.add_argument('input_dir', help=('Directory with the PDB files.'),
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    parser.add_argument('--cpus', help='Number of CPUs to use', type=int, default=4)
    
    return parser.parse_args(args)


def extract_sequence(pdb_file: Path) -> str:
    """
    Extract the sequence from a PDB file
    """
    m = b.PDBModel(str(pdb_file))
    sequence = m.sequence()
    
    return pdb_file, sequence


if __name__=='__main__':
    args = parsing()
    
    pdb_files = list(args.input_dir.glob('*.pdb'))
    
    with Pool(args.cpus) as p:
        results = p.map(extract_sequence, pdb_files)
    
    sequences = {str(p.name): s for p, s in results}
    
    with open(args.output, 'wb') as f:
        pickle.dump(sequences, f)
    
    print(f'Sequences extracted from {len(sequences)} PDB files')