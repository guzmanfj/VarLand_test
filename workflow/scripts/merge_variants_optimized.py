import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import os
import multiprocessing as mp

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
    
    parser.add_argument('--cpus', help='Number of CPUs to use', type=int, default=4)
    
    return parser.parse_args(args)

def read_pickle(file):
    return pd.read_pickle(file)

def process_contacts(var):
    res_pos = var['Residue_position']
    contacts = var['intra_contacts']
    select = np.abs(res_pos - contacts) >= 6
    return contacts[select]

def separate_values(df):
    columns = df.columns[6:13]
    return df.apply(lambda row: pd.DataFrame({
        col: row[col].split(';') for col in columns
    }).assign(**{col: row[col] for col in df.columns[:6]}), axis=1).to_list()

def process_dbnsfp_file(file):
    if os.stat(file).st_size == 0:
        return None
    df = pd.read_csv(file, sep='\t')
    df['aapos'] = df['aapos'].astype(str)

    try:
        # Apply separate_values to this individual file
        df_exploded = pd.concat(separate_values(df.iloc[:,:13])).reset_index(drop=True)
        
        # Merge the exploded dataframe with the score columns of the original df
        df_merged = pd.merge(df_exploded, df.drop(columns=df.columns[6:13]),
                            on=['#chr','pos(1-based)','ref','alt','aaref','aaalt'], how='left')
        
        # Rename columns to merge with the variants dataframe
        df_merged.rename(columns={'#chr':'#CHROM',
                                'pos(1-based)':'POS',
                                'ref':'REF',
                                'alt':'ALT',
                                'Ensembl_transcriptid':'Feature',
                                'aapos':'Residue_position'}, inplace=True)
        df_merged['Residue_position'] = df_merged['Residue_position'].astype(int)
        df_merged['#CHROM'] = 'chr' + df_merged['#CHROM'].astype(str)
    except:
        print(file)
        raise
    
    return df_merged

if __name__ == '__main__':
    args = parsing()
    
    # Use multiprocessing to read pickle files
    with mp.Pool(args.cpus) as pool:
        variants = pd.concat(pool.map(read_pickle, args.input.glob('15_variants_*.pkl')))
    
    variants['secondary_structure'] = variants['secondary_structure'].replace(' ', np.nan)
    
    # Vectorize contacts processing
    variants['intra_contacts'] = variants.apply(process_contacts, axis=1)
    
    # Read FoldX energies in parallel
    with mp.Pool(args.cpus) as pool:
        foldx_energies = {}
        for result in pool.map(read_pickle, args.foldx_energies.glob('17_foldx_energies_*.pkl')):
            foldx_energies.update(result)
    
    foldx = pd.DataFrame(foldx_energies).T
    variants = variants.merge(foldx, left_index=True, right_index=True)
    
    select_columns = [
        '#CHROM', 'POS', 'REF', 'ALT', 'RefSeq', 'PDB_path', 'Amino_acids',
        'Residue_position', 'pLDDT'
    ]
    variants = variants.drop_duplicates(subset=select_columns).reset_index(drop=True)
    
    # Process dbNSFP files in parallel
    with mp.Pool(args.cpus) as pool:
        dbnsfp_results = pool.map(process_dbnsfp_file, args.dbnsfp_dir.glob('*.out'))
    
    dbnsfp_merged = pd.concat([df for df in dbnsfp_results if df is not None]).reset_index(drop=True)

    # Concatenate all processed dbNSFP results
    dbnsfp_merged = pd.concat(dbnsfp_results).reset_index(drop=True)
    
    # Merge with variants
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
    variants_merged.to_pickle(args.output)
    
    print(variants_merged.shape)
    print('Done.')
