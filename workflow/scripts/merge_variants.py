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
    
    parser.add_argument('dbnsfp_dir', help='Input directory with the dbNSFP '
                        'files', type=validate_file)
    
    parser.add_argument('output', help='Output file', type=Path)
    
    return parser.parse_args(args)


def separate_values(row:pd.Series, columns:list):
    """
    Separate the values in the columns of a row by the delimiter, and return
    a dataframe with the separated values into their own rows.

    Args:
        row (pd.Series): Row to separate.
        columns (list): Columns to consider for separation.
    """
    # Separate values
    split_values = []
    for col in columns:
        split_values.append(row[col].split(';'))
    
    # Create dictionary with index values
    index_columns = ['#chr','pos(1-based)','ref','alt','aaref','aaalt']
    df_values = {indcol:indval for indcol, indval in zip(index_columns, row[index_columns])}
    # Add the separated values to the dictionary
    df_values.update({col:val for col, val in zip(columns, split_values)})
    
    return pd.DataFrame(df_values)


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


    ############################### Add dbNSFP features
    dbnsfp_results = []
    for file in args.dbnsfp_dir.glob('*.out'):
        idbnsfp = pd.read_csv(file, sep='\t')
        dbnsfp_results.append(idbnsfp)

    dbnsfp = pd.concat(dbnsfp_results).reset_index(drop=True)
    
    # The dbnsfp dataframe sometimes contains multiple annotations for a single variant,
    # separated by a semicolon. We need to separate these values into their own rows.
    # Separate the values in the columns of a row by the delimiter, and return
    # a dataframe with the separated values into their own rows.
    dbnsfp_exploded = pd.concat(dbnsfp.iloc[:,:13].apply(lambda x: separate_values(
                                        x, list(dbnsfp.columns[6:13])), axis=1).values).reset_index(drop=True)

    # Merge the exploded dataframe with the score columns of the original dbnsfp dataframe
    dbnsfp_merged = pd.merge(dbnsfp_exploded, dbnsfp.drop(columns=dbnsfp.columns[6:13]),
                            on=['#chr','pos(1-based)','ref','alt','aaref','aaalt'], how='left')

    # Rename columns to merge with the variants dataframe
    dbnsfp_merged.rename(columns={'#chr':'#CHROM',
                                'pos(1-based)':'POS',
                                'ref':'REF',
                                'alt':'ALT',
                                'Ensembl_transcriptid':'Feature',
                                'aapos':'Residue_position'}, inplace=True)
    dbnsfp_merged['Residue_position'] = dbnsfp_merged['Residue_position'].astype(int)
    dbnsfp_merged['#CHROM'] = 'chr' + dbnsfp_merged['#CHROM'].astype(str)

    variants_merged = pd.merge(variants, dbnsfp_merged,
                               on=['#CHROM','POS','REF','ALT','Feature','Residue_position'],
                               how='left')

    cols_to_remove = ['Consequence','Feature_type','BIOTYPE','EXON','INTRON',
                      'HGVSc','HGVSp','cDNA_position','DISTANCE','Ensembl_geneid',
                      'Ensembl_proteinid','Uniprot_acc','Uniprot_entry']

    variants_merged.drop(columns=cols_to_remove, inplace=True)

    cols_to_format = ['GERP++_NR','GERP++_RS']
    for col in cols_to_format:
        variants_merged[col] = variants_merged[col].replace('.', np.nan).astype(float)
        variants_merged[col] = variants_merged[col].fillna(np.nan)

    # Save to pickle
    # variants.to_pickle('./AM_pathogenic.pkl')
    variants_merged.to_pickle(args.output)

    print(variants_merged.shape)
    print('Done.')
