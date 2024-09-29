"""
Script to calculate FoldX energies for the variants
"""

import re
import sys
import pandas as pd
from FoldX.FoldX import FoldX
import pickle
import argparse
from pathlib import Path


def parsing(args: list = None) -> argparse.Namespace:
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

    def validate_file(d: str) -> Path:
        """
        Validate that the directory with the features exists
        """
        d = Path(d)
        if not d.exists():
            raise ValueError("The specified file or directory doesn't exist.")

        return d

    parser = argparse.ArgumentParser(
        description=("Calculate FoldX energies for" " the variants pickle file.")
    )

    parser.add_argument("input", help="Input pickle file", type=validate_file)

    parser.add_argument("output", help="Output pickle file", type=Path)
    
    parser.add_argument('--bin_dir', help='FoldX binary directory', default=None,
                        type=validate_file)

    return parser.parse_args(args)


if __name__ == "__main__":

    args = parsing()
    
    if args.bin_dir is None:
        bin_dir = Path.cwd() / 'resources/foldx'
    else:
        bin_dir = args.bin_dir

    variants = pd.read_pickle(args.input)

    # Iterate over all the variants and:
    # 1. Get the PDB file in the PDB_path column
    # 2. Get the FoldX energy for the variant
    foldx_energies = {}
    for i, var in variants.iterrows():

        # Correct the residue number according to the number of the model in the PDB file
        model_no = int(re.search(r"-F(\d+)-", var.PDB_path.stem).group(1))
        resnumber = var.Residue_position - 200 * (model_no - 1)

        wt, mut = var["Amino_acids"].split("/")
        mut_string = f"{wt}{resnumber}{mut}"

        exe = FoldX(var.PDB_path, [[mut_string]], ["A"], verbose=False, bin_dir=bin_dir)

        foldx_energies[i] = exe.run()

    # Save the results to pickle
    # with open(f'./foldx_energies/17_foldx_energies_{job}.pkl', 'wb') as f:
    #     pickle.dump(foldx_energies, f)
    with open(args.output, "wb") as f:
        pickle.dump(foldx_energies, f)

    print(f"Done.")
