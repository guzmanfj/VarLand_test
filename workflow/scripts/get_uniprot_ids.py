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
    
    def validate_dir(d:str) -> Path:
        """
        Validate that the directory exists, creating it if necessary
        """
        d = Path(d)
        if not d.exists():
            d.mkdir(parents=True)
        
        return d

    parser = argparse.ArgumentParser(description=('Obtain the UniProt IDs.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('idmapping', help='UniProt ID mapping file',
                        type=validate_file )
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    parser.add_argument('output_ranges', help=('Output directory with the ranges '
                        'of variants per job. It will contain files named with the '
                        'start and end of the range of variants to process, separated '
                        'by a dash.'), type=validate_dir)
    
    return parser.parse_args(args)


def calculate_variants_per_job(nvariants:int, out_dir:Path, min_variants:int=50,
                               max_jobs:int=100) -> None:
    """
    I have to split the variants in smaller chunks to run the next script in
    parallel. This function calculates the number of variants to run in each
    job, and writes to a file a list with pairs of integers indicating the
    start and end of the range of variants to run in each job. E.g.:
    0,100
    101,200
    201,300
    
    Args:
        nvariants (int): total number of variants
        out_dir (Path): output directory to write the ranges
        min_variants (int): minimum number of variants to run in each job
        max_jobs (int): maximum number of jobs to run in parallel
    """
    
    variants_per_job = int(np.ceil(nvariants/max_jobs))
    variants_per_job = max(variants_per_job, min_variants)
    
    njobs = int(np.ceil(nvariants/variants_per_job))
    
    # Create empty files with the ranges as the name
    for i in range(njobs):
        start = i*variants_per_job
        end = min((i+1)*variants_per_job, nvariants)
        fname = out_dir/f'range_{start}-{end}'
        fname.touch()


if __name__=='__main__':

    args = parsing()

    # variants = pd.read_pickle('./05_variants.pkl')
    variants = pd.read_pickle(args.input)

    print(f'{variants.shape} initial variants')

    # Remove empty FLAGS column
    variants = variants.drop(columns=['FLAGS'])

    # See how many missing IDs in the UniProt database
    # variants.UNIPROT_ISOFORM.isna().sum()
    # Assume that there are missing ids in the UniProt database

    # Search for IDs in the UniProt database
    # Read human ID mapping from UniProt
    # human_idmapping = pd.read_csv(
    #     '/ibex/scratch/projects/c2102/databases/uniprot/2023_02/current_release/'
    #     'knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz',
    #     sep='\t', names=['UniProtKB-AC', 'ID_type', 'ID'])
    human_idmapping = pd.read_csv(args.idmapping, sep='\t',
                                  names=['UniProtKB-AC', 'ID_type', 'ID'])

    # Extract RefSeq and Ensembl IDs
    transcripts_refseq = human_idmapping[human_idmapping.ID_type == 'RefSeq_NT'].copy()
    transcripts_refseq = transcripts_refseq.drop(columns=['ID_type'])

    transcripts_ensemble = human_idmapping[human_idmapping.ID_type == 'Ensembl_TRS'].copy()
    transcripts_ensemble = transcripts_ensemble.drop(columns=['ID_type'])

    # Remove transcript version numbers from the 'RefSeq' column
    variants['RefSeq_noversion'] = variants['RefSeq'].str.split('.').str[0]

    # Remove transcript version numbers from the 'ID' columns
    transcripts_refseq['ID_noversion'] = transcripts_refseq['ID'].str.split('.').str[0]
    transcripts_ensemble['ID_noversion'] = transcripts_ensemble['ID'].str.split('.').str[0]


    # Make a dictionary to store all the UniProt IDs for each variant
    uids = {}

    # For each variant, get the UniProt IDs from:
    # 1. The 'UNIPROT_ISOFORM' column
    # 2. The transcripts_refseq dataframe, by looking up with the 'RefSeq_noversion' column
    # 3. The transcripts_ensemble dataframe, by looking up with the 'Feature' column
    for i, row in variants.iterrows():
        # Use a list instead of a set, because we want to preserve the order
        l = []
        if not pd.isna(row['UNIPROT_ISOFORM']):
            l.append(row['UNIPROT_ISOFORM'])

        l.extend(transcripts_refseq[
            transcripts_refseq.ID_noversion == row['RefSeq_noversion']]['UniProtKB-AC'].values)

        l.extend(transcripts_ensemble[
            transcripts_ensemble.ID_noversion == row['Feature']]['UniProtKB-AC'].values)
        
        uids[i] = list(dict.fromkeys(l))

    variants['UniProt_IDs'] = pd.Series(uids)

    # Remove variants with no UniProt IDs
    select = variants['UniProt_IDs'].apply(lambda x: len(x)>0)
    variants = variants[select]
    variants = variants.reset_index(drop=True)


    # Remove variants that have a range of amino acid positions changed (e.g. 1-3)
    select = (variants.Protein_position.str.split('/').str[0].str.split('-')
            .apply(lambda x:len(x)>1))

    variants = variants[~select]
    variants = variants.reset_index(drop=True)


    print(f'{variants.shape} variants with UniProt IDs')

    # variants.to_pickle('./07_variants.pkl')
    variants.to_pickle(args.output)
    
    # Calcualte the number of variants per job
    calculate_variants_per_job(variants.shape[0], args.output_ranges)

    print('Done!')

