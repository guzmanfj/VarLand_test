import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing

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
    
    parser.add_argument('output_ranges', help=('Output directory with the ranges '
                        'of variants per job. It will contain files named with the '
                        'start and end of the range of variants to process, separated '
                        'by a dash.'), type=validate_dir)
    
    parser.add_argument('--cpus', help='Number of CPUs to use', type=int, default=1)
    
    parser.add_argument('--chunk_size', help='Chunk size to read the input file',
                        type=int, default=100000)
    
    parser.add_argument('--max_jobs', help=('Maximum number of files to divide '
                                'the variants into.'), type=int, default=100)
    
    return parser.parse_args(args)


def calculate_variants_per_job(nvariants:int, min_variants:int=50,
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
        min_variants (int): minimum number of variants to run in each job
        max_jobs (int): maximum number of jobs to run in parallel
    """
    
    variants_per_job = int(np.ceil(nvariants/max_jobs))
    variants_per_job = max(variants_per_job, min_variants)
    
    njobs = int(np.ceil(nvariants/variants_per_job))
    
    # Create empty files with the ranges as the name
    ranges = []
    for i in range(njobs):
        start = i*variants_per_job
        end = min((i+1)*variants_per_job, nvariants)
        ranges.append((start, end))
    
    return ranges


def process_chunk(chunk, transcripts_refseq, transcripts_ensemble):
    chunk['RefSeq_noversion'] = chunk['RefSeq'].str.split('.').str[0]
    
    def get_uniprot_ids(row):
        uids = []
        if not pd.isna(row['UNIPROT_ISOFORM']):
            uids.append(row['UNIPROT_ISOFORM'])
        
        uids.extend(transcripts_refseq[
            transcripts_refseq.ID_noversion == row['RefSeq_noversion']]['UniProtKB-AC'].values)
        
        uids.extend(transcripts_ensemble[
            transcripts_ensemble.ID_noversion == row['Feature']]['UniProtKB-AC'].values)
        
        return list(dict.fromkeys(uids))
    
    chunk['UniProt_IDs'] = chunk.apply(get_uniprot_ids, axis=1)
    
    # Return only the rows with UniProt IDs
    return chunk[chunk['UniProt_IDs'].apply(len) > 0]


if __name__=='__main__':

    args = parsing()

    # Read the input file in chunks
    chunksize = args.chunk_size  # Adjust based on available memory
    chunks = pd.read_csv(args.input, chunksize=chunksize, index_col=0)

    # Preprocess ID mapping data
    human_idmapping = pd.read_csv(args.idmapping, sep='\t',
                                  names=['UniProtKB-AC', 'ID_type', 'ID'])

    # Extract RefSeq and Ensembl IDs
    transcripts_refseq = human_idmapping[human_idmapping.ID_type == 'RefSeq_NT'].copy()
    transcripts_refseq = transcripts_refseq.drop(columns=['ID_type'])
    transcripts_refseq['ID_noversion'] = transcripts_refseq['ID'].str.split('.').str[0]

    transcripts_ensemble = human_idmapping[human_idmapping.ID_type == 'Ensembl_TRS'].copy()
    transcripts_ensemble = transcripts_ensemble.drop(columns=['ID_type'])
    transcripts_ensemble['ID_noversion'] = transcripts_ensemble['ID'].str.split('.').str[0]

    # Process chunks in parallel
    with multiprocessing.Pool(args.cpus) as pool:
        results = pool.starmap(process_chunk, [(chunk, transcripts_refseq, transcripts_ensemble)
                                                  for chunk in chunks])

    # Combine results
    variants = pd.concat(results, ignore_index=True)

    # Remove variants with a range of amino acid positions changed (e.g. 1-3)
    variants = variants[~variants.Protein_position.str.split('/').str[0].str.contains('-')]
    variants = variants.reset_index(drop=True)

    print(f'{variants.shape[0]} variants with UniProt IDs')
    
    # Calcualte the number of variants per job
    ranges = calculate_variants_per_job(variants.shape[0], max_jobs=args.max_jobs)
    
    # Save split variants according to the ranges
    for i, (start, end) in enumerate(ranges):
        variants.iloc[start:end].to_pickle(args.output_ranges/f'07_variants_{start}-{end}.pkl')

    print('Done!')

