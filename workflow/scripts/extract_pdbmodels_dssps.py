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

    parser.add_argument('input_dir', help=('Directory with the databases.'),
                        type=validate_file)
    
    parser.add_argument('--cpus', help='Number of CPUs to use', type=int, default=4)
    
    return parser.parse_args(args)


def read_PDBModel(pdb_file: Path) -> str:
    """
    Read the PDBModel object from a PDB file
    """
    m = b.PDBModel(str(pdb_file))
    
    return pdb_file, m


def read_DSSP(dssp_file: Path) -> str:
    """
    Read the DSSP secondary structure and accessibility data
    """
    # Get the secondary structure of the residue from DSSP file
    with open (dssp_file, 'r') as f:
        lines = f.readlines()[28:]
    
    ss = [line[16] for line in lines]
    acc = [int(line[35:38].strip()) for line in lines]
    
    return dssp_file, {'ss': ss, 'acc': acc}
    

if __name__=='__main__':
    args = parsing()
    af_database = args.input_dir / 'alphafold_human_v4'
    af_mane_database = args.input_dir / 'alphafold_mane'
    pdb_files = list(af_database.glob('*.pdb')) + list(af_mane_database.glob('*.pdb'))
    
    with Pool(args.cpus) as p:
        results = p.map(read_PDBModel, pdb_files)
    
    pdbmodels = {str(p): m for p, m in results}
    
    with open(args.input_dir / 'pdbmodels.pkl', 'wb') as f:
        pickle.dump(pdbmodels, f)
    
    print(f'PDBModels extracted from {len(pdbmodels)} PDB files')
    
    af_dssp = args.input_dir / 'alphafold_human_v4_dssp'
    af_mane_dssp = args.input_dir / 'alphafold_mane_dssp'
    dssp_files = list(af_dssp.glob('*.dssp')) + list(af_mane_dssp.glob('*.dssp'))
    
    with Pool(args.cpus) as p:
        results = p.map(read_DSSP, dssp_files)
    
    dssps = {str(p): dssp for p, dssp in results}
    
    with open(args.input_dir / 'dssps.pkl', 'wb') as f:
        pickle.dump(dssps, f)
    
    print(f'DSSP data extracted from {len(dssps)} DSSP files')