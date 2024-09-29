"""
Script to annotate physicochemical properties of protein variants
"""

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
    

    parser = argparse.ArgumentParser(description=('Annotate physicochemical properties'
                                                  ' of protein variants.'))

    parser.add_argument('input', help='Input variant pickle file',
                        type=validate_file)
    
    parser.add_argument('panther_db', help='Pickle file with PANTHER database',
                        type=validate_file)
    
    parser.add_argument('output', help='Output file', type=Path)
    
    return parser.parse_args(args)


aa_dict = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E',
    'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
    'Tyr': 'Y', 'Val': 'V'
    }


def getting_ref_aa_properties(variants:pd.DataFrame, AA_GIVEN:bool=False):
    """
    Function to assign the properties threshold to the variants input file

    parameters
    ------------------
    path_variants_input: pd.DataFrame
        DataFrame with the variants
    AA_GIVEN: bool
        if the aminoacid ref is in column REFAA and aa alt in column ALTAA is
        given in the input file as one letter    
    """

    if AA_GIVEN == False:
        # Extract the amino acids and position from the HGVS_p column
        variants[['REFAA', 'pos', 'ALTAA']] = (
                                    variants['HGVS_p'].str
                                    .extract(
                                    r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})'))

        # Map the three letter aminoacid code to one letter
        variants['ALTAA'] = variants['ALTAA'].map(aa_dict)
        variants['REFAA'] = variants['REFAA'].map(aa_dict)

    ## ENERGY
    variants['Stabilizing Energy Change'] = (
        (variants['total energy'])
        .apply(lambda x: 'Stabilizing Energy Change' 
               if x<0 else 
               'not_Stabilizing Energy Change'))
    
    variants['Destabilizing Energy Change (<1 kcal/mol)'] = (
        (variants['total energy'])
        .apply(lambda x: 'Destabilizing Energy Change (<1 kcal/mol)' 
               if x > 0 and x < 1 else 
               'not_Destabilizing Energy Change (<1 kcal/mol)'))
    
    variants['Destabilizing Energy Change (1-3 kcal/mol)'] = (
        (variants['total energy'])
        .apply(lambda x: 'Destabilizing Energy Change (1-3 kcal/mol)' 
               if x > 1 and x < 3 else 
               'not_Destabilizing Energy Change (1-3 kcal/mol)'))
    
    variants['Destabilizing Energy Change (>3 kcal/mol)'] = (
        (variants['total energy'])
        .apply(lambda x: 'Destabilizing Energy Change (>3 kcal/mol)' 
               if x>3 else 
               'not_Destabilizing Energy Change (>3 kcal/mol)'))

    variants['TotalEnergy'] = (
        (variants['total energy'])
        .apply(lambda x: 'TotalEnergy' if x>3 else 'not_TotalEnergy'))
    
    variants['LessTotalEnergy'] = (
        (variants['total energy'])
        .apply(lambda x: 'LessTotalEnergy' if x<0 else 'not_LessTotalEnergy'))
    
    variants['VanDerWaalsClashes'] = (
        (variants['Van der Waals clashes'])
        .apply(lambda x: 'VanDerWaalsClashes' if x>0.9 or x < -0.003 
               else 'not_VanDerWaalsClashes'))
    
    variants['Disulfide'] = (
        (variants['disulfide'])
        .apply(lambda x: 'Disulfide' if x>0.1 or x < 0 else 'not_Disulfide'))
   
    ## CONTACTS
    variants['Contacts'] = (
        (variants['intra_contacts'])
        .apply(lambda x: 'Contacts' if len(x) != 0 else 'not_Contacts'))
    variants['No Contacts'] = (
        (variants['intra_contacts'])
        .apply(lambda x: 'No Contacts' if len(x) == 0 else 'not_No Contacts'))


    ## CONSERVATION
    variants['Conserved'] = (
        (variants['GERP++_RS'])
        .apply(lambda x: 'Conserved' if x > 2 else 'not_Conserved'))

    ## CATALYTIC/DOMAIN NOT USED
    variants['Catalytic'] = (
        (variants['is_catalytic'])
        .apply(lambda x: 'Catalytic' if x == True else 'not_Catalytic'))
    variants['Domain'] = (
        (variants['DOMAINS'])
        .apply(lambda x: 'Catalytic' if x == True else 'not_Domain'))

    ## ORDER
    variants['DisorderpLDDT'] = (
        (variants['pLDDT'])
        .apply(lambda x: 'DisorderpLDDT' if x < 50 else 'not_DisorderpLDDT'))
    variants['OrderpLDDT'] = (
        (variants['pLDDT'])
        .apply(lambda x: 'OrderpLDDT' if x > 50 else 'not_OrderpLDDT'))

    ## EXPOSURE
    variants['Core (<5%)'] = (
        (variants['ACCESIBILITY_NORMALIZED'])
        .apply(lambda x: 'Core (<5%)' if x < 0.05 else 'not_Core (<5%)'))
    variants['Buried (5-25%)'] = (
        (variants['ACCESIBILITY_NORMALIZED'])
        .apply(lambda x: 'Buried (5-25%)' if x > 0.05 and x < 0.25 
               else 'not_Buried (5-25%)'))
    variants['Medium-buried (25-50%)'] = (
        (variants['ACCESIBILITY_NORMALIZED'])
        .apply(lambda x: 'Medium-buried (25-50%)' if x > 0.25 and x < 0.5 
               else 'not_Medium-buried (25-50%)'))
    variants['Medium-exposed (50-75%)'] = (
        (variants['ACCESIBILITY_NORMALIZED'])
        .apply(lambda x: 'Medium-exposed (50-75%)' if x > 0.5 and x < 0.75 
               else 'not_Medium-exposed (50-75%)'))
    variants['Exposed (>75%)'] = (
        (variants['ACCESIBILITY_NORMALIZED'])
        .apply(lambda x: 'Exposed (>75%)' if x > 0.75 else 'not_Exposed (>75%)'))

    ## SECONDARY STRUCTURE
    helices = ['H','G','I']
    βetas = ['B','E']
    coils = ['L','S','T']  

    variants['helices'] = (
        variants.apply(lambda row: 'helices' 
                       if (row['secondary_structure'] in helices) 
                       else 'not_helices', axis=1))
    variants['β-sheet/strand'] = (
        variants.apply(lambda row: 'β-sheet/strand' 
                       if (row['secondary_structure'] in βetas) 
                       else 'not_β-sheet/strand', axis=1))
    variants['coils'] = (
        variants.apply(lambda row: 'coils' 
                       if (row['secondary_structure'] in coils) 
                       else 'not_coils', axis=1))
    variants['α-helix'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'α-helix' if x == "H" else 'not_α-helix'))
    variants['βbridge'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'βbridge' if x == "B" else 'not_βbridge'))
    variants['strandβladder'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'strandβladder' if x == "E" else 'not_strandβladder'))
    variants['310helix'] = (
        (variants['secondary_structure'])
        .apply(lambda x: '310helix' if x == "G" else 'not_310helix'))
    variants['π-helix'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'π-helix' if x == "I" else 'not_π-helix'))
    variants['HBondTurn'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'HBondTurn' if x == "T" else 'not_HBondTurn'))
    variants['Bend'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'Bend' if x == "S" else 'not_Bend'))
    variants['Loop'] = (
        (variants['secondary_structure'])
        .apply(lambda x: 'Loop' if x == "L" else 'not_Loop'))


    ## PHYSICOCHEMICAL CHANGE
    small = ['G','A','S']
    big = ['F','W','Y','R','K','L','I','M']
    nonpolar = ['G','A','V','C','P','L','I','M','W','F']
    positive = ['H','K','R']
    negative = ['D','E']
    charged = positive + negative
    polar_uncharged = ['S','T','Y','N','Q']
    polar = positive + negative + polar_uncharged
    aromatic = ['F','W','Y']
    uncharged = nonpolar + polar_uncharged
    hydrophobic = ['A', 'C', 'I', 'L', 'M', 'F', 'W', 'V']
    hydrophilic = ['R', 'N', 'D', 'Q', 'E', 'K']
    neutral = ['G', 'H', 'P', 'S', 'T', 'Y']

    # create pysicochemical classifications REF aminoacid and ALT aminoacid
    variants['Small to Big'] = (
        variants.apply(lambda row: 'Small to Big' 
                       if row['REFAA'] in small and row['ALTAA'] in big 
                       else 'not_Small to Big', axis=1))
    variants['Big to Small'] = (
        variants.apply(lambda row: 'Big to Small' 
                       if row['REFAA'] in big and row['ALTAA'] in small 
                       else 'not_Big to Small', axis=1))
    variants['Polar to NonPolar'] = (
        variants.apply(lambda row: 'Polar to NonPolar' 
                       if row['REFAA'] in polar and row['ALTAA'] in nonpolar 
                       else 'not_Polar to NonPolar', axis=1))
    variants['NonPolar to Polar'] = (
        variants.apply(lambda row: 'NonPolar to Polar' 
                       if row['REFAA'] in nonpolar and row['ALTAA'] in polar 
                       else 'not_NonPolar to Polar', axis=1))
    variants['Hydrophilic introduced'] = (
        variants.apply(lambda row: 'Hydrophilic introduced' 
                       if (row['REFAA'] in hydrophobic or row['REFAA'] in neutral) 
                       and (row['ALTAA'] in hydrophilic) 
                       else 'not_Hydrophilic introduced', axis=1))
    variants['Hydrophobic introduced'] = (
        variants.apply(lambda row: 'Hydrophobic introduced' 
                       if (row['REFAA'] in hydrophilic or row['REFAA'] in neutral) 
                       and (row['ALTAA'] in hydrophobic or row['ALTAA'] in neutral) 
                       else 'not_Hydrophobic introduced', axis=1))
    variants['Charge switch'] = (
        variants.apply(lambda row: 'Charge switch' 
                       if (row['REFAA'] in positive and row['ALTAA'] in negative) 
                       or (row['REFAA'] in negative and row['ALTAA'] in positive) 
                       else 'not_Charge switch', axis=1))
    variants['Charge lost'] = (
        variants.apply(lambda row: 'Charge lost' 
                       if row['REFAA'] in charged and row['ALTAA'] in uncharged 
                       else 'not_Charge lost', axis=1))
    variants['Charge gain'] = (
        variants.apply(lambda row: 'Charge gain' 
                       if row['REFAA'] in uncharged and row['ALTAA'] in charged 
                       else 'not_Charge gain', axis=1))
    variants['Aromatic to NonAromatic'] = (
        variants.apply(lambda row: 'Aromatic to NonAromatic' 
                       if row['REFAA'] in aromatic and row['ALTAA'] not in aromatic 
                       else 'not_Aromatic to NonAromatic', axis=1))
    variants['Aromatic to polar'] = (
        variants.apply(lambda row: 'Aromatic to polar' 
                       if row['REFAA'] in aromatic and row['ALTAA'] not in polar 
                       else 'not_Aromatic to polar', axis=1))
    
    return variants


def get_protein_classes(uids:list, symbol:str, panther:pd.DataFrame):
    """
    Get the UniProt ID that matches the ProteinID in the PANTHER database.

    Args:
        uids (list): List of UniProt IDs.
        symbol (str): Gene symbol.
        panther (pd.DataFrame): DataFrame with the ProteinID and corresponding
            GeneralProteinClass.
    """
    if isinstance(symbol,str):
        select_matching = (panther.ProteinID.isin(uids)) | (panther.Gene==symbol)
    else:
        select_matching = (panther.ProteinID.isin(uids))
        
    matching_panthers = panther[select_matching]
    protein_classes = list(matching_panthers.GeneralProteinClass.unique())
    if len(protein_classes) == 0:
        return np.nan
    elif len(protein_classes) == 1:
        return protein_classes[0]
    elif len(protein_classes) > 1:
        print(f"Multiple protein classes found for {uids}: {protein_classes}."
              "Returning the first one.")
        return protein_classes[0]



if __name__ == '__main__':
    
    args = parsing()
    
    variants = pd.read_pickle(args.input)
    
    ############### from 01_data_preparation.py ################
    
    # Split the amino acids into two columns
    variants[['REFAA', 'ALTAA']] = (variants['Amino_acids']
                                           .str.split('/', expand=True))

    # Replace variants.secondary_structure == 'nan' with 'L'
    variants.secondary_structure = (variants.secondary_structure
                                           .replace(np.nan, 'L'))
    
    # Normalize residue accessibility
    # Dictionary mapping amino acids to their respective x values
    x_values = {
        'Ala': 129.0,
        'Arg': 274.0,
        'Asn': 195.0,
        'Asp': 193.0,
        'Cys': 167.0,
        'Glu': 223.0,
        'Gln': 225.0,
        'Gly': 104.0,
        'His': 224.0,
        'Ile': 197.0,
        'Leu': 201.0,
        'Lys': 236.0,
        'Met': 224.0,
        'Phe': 240.0,
        'Pro': 159.0,
        'Ser': 155.0,
        'Thr': 172.0,
        'Trp': 285.0,
        'Tyr': 263.0,
        'Val': 174.0
    }

    # Create the 'ACCESIBILITY_NORMALIZED' column
    variants['ACCESIBILITY_NORMALIZED'] = variants.apply(
            lambda row: row['accessibility'] / x_values[row['Residue']], axis=1)
    
    ############### from 02_assignment_wrapper_physicochem.py ################
    
    # Assign binary features to the variants
    AA_GIVEN = True
    variants = getting_ref_aa_properties(variants, AA_GIVEN)

    # Add the protein class to the dataframe
    panther = pd.read_pickle(args.panther_db)
    
    generalproteinclass = []
    for idx, row in variants.iterrows():
        generalproteinclass.append(get_protein_classes(row['UniProt_IDs'],
                                                       row['SYMBOL'],
                                                       panther))
    
    variants['GeneralProteinClass'] = generalproteinclass
    
    # Save the dataframe with the features
    variants.to_pickle(args.output)
    
    print(f"Physicochemical features added to {args.output}")
    print(f"Shape of the dataframe: {variants.shape}")
