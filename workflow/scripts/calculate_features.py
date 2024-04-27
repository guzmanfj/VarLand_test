'''
Script to calcualte the features for the variants in the ClinVar dataset:
- Intra-molecular contacts
- Catalytic residues
- pLDDT
- Secondary structure
- Accessibility
'''

import biskit as b
import pandas as pd
import numpy as np
import re
from pathlib import Path
import argparse

from structure_features.contactsPdb import get_contacts


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
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


if __name__ == '__main__':
    
    args = parsing()

    # Get the number of the job from the arguments.
    # job = int(sys.argv[1])
    # commands_per_job = int(sys.argv[2])

    # Check if the results for this job have already been calculated
    # if Path(f'./variants_features/15_variants_{job}.pkl').exists():
    #     print(f'Job {job} already done.')
    #     sys.exit()

    # variants = pd.read_pickle(f'./14_variants/14_variants_{job}.pkl')
    variants = pd.read_pickle(args.input)
    # catalytic_sites = pd.read_pickle('./13_catalytic_residues.pkl')
    catalytic_sites = pd.read_pickle(args.catalytic_sites_pickle)

    # Get the variants for this job
    # indices = np.arange(0, len(variants), commands_per_job)
    # start = indices[job]
    # end = start + commands_per_job
    # variants = variants.iloc[start:end]

    # Iterate over all the variants and:
    # 1. Read the PDB file in the PDB_path column
    # 2. Get the intra-molecular contacts for the given residue
    intra_contacts = {}
    is_catalytic = {}
    plddts = {}
    secondary_structure = {}
    accessibility = {}
    for i, var in variants.iterrows():
        
        # Correct the residue number according to the number of the model in the PDB file
        model_no = int(re.search(r'-F(\d+)-', var.PDB_path.stem).group(1))
        resnumber = var.Residue_position - 200 * (model_no - 1)
        
        m = b.PDBModel(str(var.PDB_path))
        
        # Calculate the average of the pLDDT values for a 5-residue window around the residue
        plddt = m.atom2resProfile('temperature_factor')
        if resnumber < 3:
            plddt = plddt[0:resnumber+2].mean()
        elif resnumber > (len(plddt) - 3):
            plddt = plddt[resnumber-3:].mean()
        else:
            plddt = plddt[resnumber-3:resnumber+2].mean()
        
        plddts[i] = plddt
        
        # Calculate intra-molecular contacts 6 Angstroms around the residue
        atom_contacts = get_contacts(m, resnumber, 6)
        # The residues are 1-indexed
        residue_contacts = np.unique(m.atoms['residue_number'][atom_contacts])
        # Correct the residue numbers back to the original position
        residue_contacts = residue_contacts + 200 * (model_no - 1)
        
        intra_contacts[i] = residue_contacts
        
        # Calculate contacts 8 Angstroms around the residue to check for catalytic residues
        atom_contacts = get_contacts(m, resnumber, 8)
        catalytic_contacts = np.unique(m.atoms['residue_number'][atom_contacts])
        catalytic_contacts = catalytic_contacts + 200 * (model_no - 1)
        
        # See if the residue or any of the contacts are in the catalytic sites
        var_uids = [uid.split('-')[0] for uid in var.UniProt_IDs]
        protein_catalytic_sites = catalytic_sites[catalytic_sites.uniprot_ID.isin(var_uids)]
        if not protein_catalytic_sites.empty:
            if var.Residue_position in protein_catalytic_sites.res_pos.values:
                is_catalytic[i] = True
            elif any([res in protein_catalytic_sites.res_pos.values for res in catalytic_contacts]):
                is_catalytic[i] = True
            else:
                is_catalytic[i] = False
        else:
            is_catalytic[i] = False
        
        # Get the secondary structure of the residue from DSSP file
        with open (var.DSSP_path, 'r') as f:
            lines = f.readlines()
        
        ss = [line[16] for line in lines[28:]]
        ss_res = ss[resnumber-1]
        secondary_structure[i] = ss_res
        
        acc = [int(line[35:38].strip()) for line in lines[28:]]
        acc_res = acc[resnumber-1]
        accessibility[i] = acc_res

    # Add the dictionaries to the dataframe
    variants['intra_contacts'] = pd.Series(intra_contacts)
    variants['is_catalytic'] = pd.Series(is_catalytic)
    variants['pLDDT'] = pd.Series(plddts)
    variants['secondary_structure'] = pd.Series(secondary_structure)
    variants['accessibility'] = pd.Series(accessibility)

    # Save to pickle
    # variants.to_pickle(f'./variants_features/15_variants_{job}.pkl')
    variants.to_pickle(args.output)

    print('Done.')

