'''
Script to compress vcf file
'''

import pysam
import argparse
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
    
    def validate_vcf(d:str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            raise ValueError("The specified file doesn't exist.")
        
        # Check .vcf extension
        if d.suffix != '.vcf':
            raise ValueError("The specified file is not a .vcf file.")
            
        return d
    

    parser = argparse.ArgumentParser(description=('Takes a .vcf file and '
            'compresses it into a .vcf.bgz file.'))

    parser.add_argument('input', help='Input .vcf file', type=validate_vcf)
    
    parser.add_argument('output', help='Output .vcf.bgz file', type=Path)
    
    return parser.parse_args(args)


if __name__=='__main__':
    
    args = parsing()
    fin = str(args.input.absolute())
    fout = str(args.output.absolute())
    
    # Compress
    pysam.tabix_compress(fin, fout)

    print(f'File {fin} compressed into {fout}')