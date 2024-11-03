import argparse
import pandas as pd
from pathlib import Path
from typing import List
import re
import numpy as np


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
    
    def validate_dir(d:str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            raise ValueError("The specified directory doesn't exist.")
            
        return d
    
    def validate_file(d:str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            raise ValueError("The specified file doesn't exist.")
            
        return d

    parser = argparse.ArgumentParser(description=('Takes a directory containing '
                'subdirectories with the .vep files for each chromosome. Reads '
                'all .vep files and concatenates them into a single dataframe.'))

    parser.add_argument('input', help='Input directory', type=validate_dir)
    
    parser.add_argument('info_fields', help='File with the INFO fields definitions',
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


def read_vcfs(chr_dirs: List[Path]) -> pd.DataFrame:
    """
    Read all .vep files in the chromosome directories and concatenate them into
    a single dataframe

    Args:
        chr_dirs (List[Path]): List of directories containing the .vep files

    Returns:
        pd.DataFrame
    """
    # Read one of the files and find the first line that starts with '#CHROM'
    with open(chr_dirs[0] / 'chunk_aa.vcf.vep','r') as f:
        for i, line in enumerate(f):
            if line.startswith('#CHROM'):
                header_line = i
    
    # Read all files ending with .vep for each chromosome directory
    chrom_dfs = []
    for d in chr_dirs:
        print(f'Reading {d}')
        # Data starts at line 38 (1-based)
        df = pd.concat([pd.read_csv(f, sep='\t', header=header_line) for f in d.glob('*.vep')],
                    ignore_index=True)
        chrom_dfs.append(df)

    # Concatenate all dataframes into one
    df = pd.concat(chrom_dfs, ignore_index=True)

    return df


def process_infos(variants: pd.DataFrame,
                  info_fields: Path) -> pd.DataFrame:
    """
    Process the INFO field from the VEP output

    Args:
        variants (pd.DataFrame): DataFrame with the VEP output
        info_fields (Path): Path to the file with the INFO fields definitions

    Returns:
        pd.DataFrame
    """
    # Modify dtypes
    variants['#CHROM'] = variants['#CHROM'].astype(str)
    variants['REF'] = variants['REF'].astype(str)
    variants['ALT'] = variants['ALT'].astype(str)
    variants['QUAL'] = variants['FILTER'].astype(str)
    variants['FILTER'] = variants['FILTER'].astype(str)
    variants['INFO'] = variants['INFO'].astype(str)

    # Erase useless columns
    variants = variants.drop(columns=['QUAL', 'FILTER'])

    # Process INFO fields
    # Read the field definitions from the info_fields.txt file
    # Get the field name and field type from each line
    # Save the names and types in a dictionary
    fields = {}
    with open(info_fields, 'r') as f:
        for line in f:
            match = re.match(r'##INFO=<ID=(\w+),Number=.,Type=(\w+),Description',
                            line)
            try:
                fields[match.group(1)] = match.group(2)
            except AttributeError:
                print(line)
                raise
            
    # Split INFO fields into lists of fields
    infos = variants.INFO.str.split(';')

    # Create a dictionary that will have all the field's contents for all of the
    # variants. The fields with a type of Flag will be stored into a list (no fields)
    info_dict = {f : [] for f in fields if fields[f] != 'Flag'}

    # No fields have a type of `Flag`, but just copy pasted from the gnomad notebook
    info_dict['flags'] = []

    # Process fields for all the variants
    # Go through the info fields for each variant
    for info_list in infos:
        variant_infos = {'flags': []}
        # split each info field into key and value
        for field in info_list:
            key = field.split('=')[0]

            # Check if the field has a type of 'Flag'
            is_flag = (fields[key] == 'Flag') if key in fields else False
            if is_flag:
                variant_infos['flags'].append(key)
            else:
                variant_infos[key] = field.split('=')[1]

        # Save the values only for the info fields that were defined in info_fields.txt
        for key in info_dict.keys():
            if key in variant_infos.keys():
                info_dict[key].append(variant_infos[key])
            else:
                info_dict[key].append(np.nan)

    # Merge DataFrames
    variants = pd.merge(variants.drop(columns=['INFO']),
                        pd.DataFrame(info_dict), left_index=True, right_index=True)

    return variants

def process_vep_fields(variants: pd.DataFrame) -> pd.DataFrame:
    """
    Process the VEP fields

    Args:
        variants (pd.DataFrame): DataFrame with the VEP output

    Returns:
        pd.DataFrame
    """
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
    
    return variants

def filter_missense(variants: pd.DataFrame) -> pd.DataFrame:
    """
    Filter the variants that contain the 'missense_variant' consequence in the
    'Consequence' column

    Args:
        variants (pd.DataFrame): DataFrame with the VEP output

    Returns:
        pd.DataFrame
    """
    # Select variants that contain the 'missense_variant' consequence in the
    # 'Consequence' column
    select = variants['Consequence'].str.contains('missense_variant')
    variants = variants[select]

    return variants

def filter_mane(variants: pd.DataFrame) -> pd.DataFrame:
    """
    Filter MANE variants

    Args:
        variants (pd.DataFrame): DataFrame with the VEP output

    Returns:
        pd.DataFrame
    """
    # Filter MANE variants
    select = np.logical_or(variants.MANE_SELECT.notna(),variants.MANE_PLUS_CLINICAL.notna())
    variants = variants[select]

    return variants

def filter_refseq(variants: pd.DataFrame) -> pd.DataFrame:
    """
    Filter RefSeq variants

    Args:
        variants (pd.DataFrame): DataFrame with the VEP output

    Returns:
        pd.DataFrame
    """
    # Split the `RefSeq` column and explode to create a row for each transcript
    variants['RefSeq_list'] = variants.RefSeq.str.split('&')
    variants = variants.explode('RefSeq_list')
    variants = variants.reset_index(drop=True)

    # Select the variants where the RefSeq ID is in either the MANE_SELECT or
    # MANE_PLUS_CLINICAL columns
    select_mane = variants['RefSeq_list'] == variants['MANE_SELECT']
    select_mane_plus = variants['RefSeq_list'] == variants['MANE_PLUS_CLINICAL']
    select = np.logical_or(select_mane, select_mane_plus)

    variants = variants[select]

    # Replace `RefSeq` column with column with unique RefSeq ID
    variants['RefSeq'] = variants['RefSeq_list']
    variants = variants.drop(columns=['RefSeq_list'])
    
    return variants


if __name__=='__main__':
    
    args = parsing()

    # Define the names of the chromosome directories
    chr_dirs = list(args.input.glob('*_split'))

    variants = read_vcfs(chr_dirs)
    variants = process_infos(variants, args.info_fields)
    variants = process_vep_fields(variants)
    variants = filter_missense(variants)
    variants = filter_mane(variants)
    variants = filter_refseq(variants)
    
    # Filter out variants with NA in the 'Amino_acids' column
    variants = variants.dropna(subset=['Amino_acids'])
    
    variants = variants.reset_index(drop=True)
    print(f'Number of variants: {len(variants)}')
    print(variants.shape)
    
    variants.to_csv(args.output, index=True)
    
    print(f'Variants saved to {args.output}')
