import argparse
import pandas as pd
import numpy as np
from pathlib import Path


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
            raise ValueError("The specified file doesn't exist.")
            
        return d
    

    parser = argparse.ArgumentParser(description=('Takes the pickle with the '
                'variant annotations from VEP, and separates the INFO field into '
                'individual columns. Also filters missense variants. '
                'Saves the resulting DataFrame into a pickle.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


if __name__=='__main__':

    args = parsing()


    # variants = pd.read_pickle('./02_variants.pkl')
    variants = pd.read_pickle(args.input)

    # Process VEP fields
    variants['CSQ'] = variants['CSQ'].str.split(',')
    variants = variants.explode('CSQ')
    variants = variants.reset_index(drop=True)

    # Read columns for VEP fields
    columns_vep = ('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|'
                'BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|'
                'Protein_position|Amino_acids|Codons|Existing_variation|'
                'DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|'
                'CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|ENSP|SWISSPROT|'
                'TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|SIFT|PolyPhen|DOMAINS|'
                'CLIN_SIG|SOMATIC|PHENO').split('|')

    vep = variants.CSQ.str.split('|', expand=True)
    vep.columns = columns_vep

    # Replace empty strings with NaN
    vep = vep.replace('', np.nan)
    variants = pd.merge(variants.drop(columns=['CSQ']),vep, left_index=True,
                        right_index=True)

    # Select variants that contain the 'missense_variant' consequence in the
    # 'Consequence' column
    select = variants['Consequence'].str.contains('missense_variant')
    variants = variants[select]
    variants = variants.reset_index(drop=True)

    # variants.to_pickle('./03_missense.pkl')
    variants.to_pickle(args.output)

    print('Done.')

