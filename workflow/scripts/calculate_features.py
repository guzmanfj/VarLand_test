'''
Script to calcualte the features for the variants in the ClinVar dataset:
- Intra-molecular contacts
- Catalytic residues
- pLDDT
- Secondary structure
- Accessibility
'''

import pandas as pd
import numpy as np
import re
from pathlib import Path
import argparse
from typing import Tuple, Dict, List
from biskit import PDBModel
from multiprocessing import Pool
from functools import partial
import copy

import biskit.mathUtils as mu


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
    

    parser = argparse.ArgumentParser(description=('Calculate features for the '
                                                  'variants pickle file.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('catalytic_sites_pickle', help=('Path to the catalytic '
                        'sites pickle file'), type=validate_file)
    
    parser.add_argument('pdbmodels_pickle', help=('Path to the PDBModel objects '
                        'pickle file'), type=validate_file)
    
    parser.add_argument('dssps_pickle', help=('Path to the DSSP data pickle file'),
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


def get_contacts(chain_m, pdb_resNumber, distances=[6]):
    """
    Gets residue contacts 

    Input
    -----
    args : chain_m (biskit.pdbModel.PDBModel)
                'biskit.pdbModel.PDBModel' type variable with information of 
                the chain model
        
    Returns :
              numpy.ndarray
    """
    
    # get index of target amino acid position
    res_numbers = chain_m.atom2resProfile('residue_number')
    res = (res_numbers == pdb_resNumber).nonzero()[0][0]
    
    # get target amino acid indices
    atoms_indices = chain_m.res2atomIndices([res])
    # get target amino acid coordinates
    atoms_coor = np.take(chain_m.xyz, atoms_indices, axis = 0)
    
    # Calculate distance between atoms of interest and all other atoms in chain
    dist = mu.pairwiseDistances(atoms_coor, chain_m.xyz)
    
    contacts = []
    for d in distances:
        # Take atoms with at most `distance` angstroms of distance to the target amino acid atoms
        contacts.append(np.unique((dist < d).nonzero()[1]))
    
    return contacts


def process_variant(row_tup: Tuple[int, pd.Series],
                    pdbmodels: Dict[str, PDBModel],
                    dssps: Dict[str, Dict[str, str]],
                    catalytic_sites: pd.DataFrame,
                    ) -> Tuple[int, float, List[int], bool, str, int]:
    """
    Process a single variant
    """
    ind = row_tup[0]
    var = row_tup[1]

    # Correct the residue number according to the number of the model in the PDB file
    model_no = int(re.search(r'-F(\d+)-', var.PDB_path.stem).group(1))
    resnumber = var.Residue_position - 200 * (model_no - 1)
    
    m = copy.deepcopy(pdbmodels[str(var.PDB_path)])
    
    # Calculate the average of the pLDDT values for a 5-residue window around the residue
    plddt = m.atom2resProfile('temperature_factor')
    if resnumber < 3:
        plddt = plddt[0:resnumber+2].mean()
    elif resnumber > (len(plddt) - 3):
        plddt = plddt[resnumber-3:].mean()
    else:
        plddt = plddt[resnumber-3:resnumber+2].mean()
    
    # Calculate intra-molecular contacts 6 Angstroms around the residue
    atom_contacts_6, atom_contacts_8 = get_contacts(m, resnumber, [6,8])
    # The residues are 1-indexed
    residue_contacts = np.unique(m.atoms['residue_number'][atom_contacts_6])
    # Correct the residue numbers back to the original position
    residue_contacts = residue_contacts + 200 * (model_no - 1)
    
    # Calculate contacts 8 Angstroms around the residue to check for catalytic residues
    catalytic_contacts = np.unique(m.atoms['residue_number'][atom_contacts_8])
    catalytic_contacts = catalytic_contacts + 200 * (model_no - 1)
    
    # See if the residue or any of the contacts are in the catalytic sites
    var_uids = [uid.split('-')[0] for uid in var.UniProt_IDs]
    protein_catalytic_sites = catalytic_sites[catalytic_sites.uniprot_ID.isin(var_uids)]
    if not protein_catalytic_sites.empty:
        if var.Residue_position in protein_catalytic_sites.res_pos.values:
            is_catalytic = True
        elif any([res in protein_catalytic_sites.res_pos.values for res in catalytic_contacts]):
            is_catalytic = True
        else:
            is_catalytic = False
    else:
        is_catalytic = False
    
    # Get the secondary structure of the residue from DSSP file
    dssp = dssps[str(var.DSSP_path)]
    
    ss = dssp['ss']
    ss_res = ss[resnumber-1]
    acc = dssp['acc']
    acc_res = acc[resnumber-1]
    
    return ind, plddt, residue_contacts, is_catalytic, ss_res, acc_res
    

def calculate_features_parallel(variants: pd.DataFrame,
                                pdbmodels: Dict[str, PDBModel],
                                dssps: Dict[str, Dict[str, str]],
                                catalytic_sites: pd.DataFrame
                                ):
    """
    Calculate the features using parallel processing
    """
    results = [process_variant(row_tup, pdbmodels, dssps, catalytic_sites) \
                for row_tup in variants.iterrows()]
    
    features = pd.DataFrame(results, columns=['index', 'pLDDT', 'intra_contacts',
                                                'is_catalytic', 'secondary_structure',
                                                'accessibility'])
    features = features.set_index('index')
    
    return features


if __name__ == '__main__':
    
    args = parsing()

    variants = pd.read_pickle(args.input)
    catalytic_sites = pd.read_pickle(args.catalytic_sites_pickle)
    pdbmodels = pd.read_pickle(args.pdbmodels_pickle)
    dssps = pd.read_pickle(args.dssps_pickle)

    # Calculate the features
    features = calculate_features_parallel(variants, pdbmodels, dssps,
                                           catalytic_sites)

    # Add the dictionaries to the dataframe
    variants['intra_contacts'] = features['intra_contacts']
    variants['is_catalytic'] = features['is_catalytic']
    variants['pLDDT'] = features['pLDDT']
    variants['secondary_structure'] = features['secondary_structure']
    variants['accessibility'] = features['accessibility']

    # Save to pickle
    variants.to_pickle(args.output)

    print('Done.')
