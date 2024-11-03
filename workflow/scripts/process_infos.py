'''
Process the INFO field of the variants
'''
import argparse
import re
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
                'individual columns. Saves the resulting DataFrame into a pickle.'))

    parser.add_argument('input', help='Input pickle file', type=validate_file)
    
    parser.add_argument('info_fields', help='File with the INFO fields definitions',
                        type=validate_file)
    
    parser.add_argument('output', help='Output pickle file', type=Path)
    
    return parser.parse_args(args)


if __name__=='__main__':

    args = parsing()

    variants = pd.read_pickle(args.input)

    print(f'Number of variants: {len(variants)}')
    print(variants.shape)

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
    with open(args.info_fields, 'r') as f:
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

        # Check if there are fields that are not present in the variant
        for key in info_dict.keys():
            if key in variant_infos.keys():
                info_dict[key].append(variant_infos[key])
            else:
                info_dict[key].append(np.nan)

    # Merge DataFrames
    variants = pd.merge(variants.drop(columns=['INFO']),
                        pd.DataFrame(info_dict), left_index=True, right_index=True)


    # variants.to_pickle('./02_variants.pkl')
    variants.to_pickle(args.output)

    print('Done.')

