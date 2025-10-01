# VarLand

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/guzmanfj/variant_annotator_snakemake.git
   ```

2. Navigate to the project directory, download, and extract the resources and data archives:
   ```bash
   cd variant_annotator_snakemake
   wget resources.tar.gz
   tar -xvzf resources.tar.gz
   rm resources.tar.gz
   wget data.tar.gz
   tar -xvzf data.tar.gz
   rm data.tar.gz
   ```

3. Download the dbNSFP database. Go to the legacy dbNSFP website (https://sites.google.com/site/jpopgen/dbNSFP) and download dbNSFP4.9a (dbNSFP4.9a.zip) into the `resources/dbnsfp` directory. Unzip the file.

4. Download FoldX. Register for an academic license at https://foldxsuite.crg.eu/academic-license-info and download `foldx5Linux64.zip` into the `resources/foldx` directory. Unzip the file.

5. Install Miniconda if you haven't already. See installation instructions at https://www.anaconda.com/docs/getting-started/miniconda/main

6. Create and activate the Conda environment:
   ```bash
   conda env create --file ./environment.yml --prefix ./env
   conda activate ./env
   ```

7. Install Singularity or Apptainer. This pipeline was tested with [Singularity version 3.9.7](https://docs.sylabs.io/guides/3.9/user-guide/quick_start.html).

## Usage

### Configure Snakemake profile

#### SLURM

The profile located in `profile/slurm/config.v8+.yaml` is set up to use SLURM as the job scheduler and Apptainer for containerization. If you want to use this profile, you need to modify the following fields:

- `slurm_account` : Set this to your SLURM account name. You can find this information by running the command `sacctmgr show user $USER` on your SLURM login node, under the column `Def Acct` (or something similar).

- `apptainer-args` : Point this to the directory where you have the VEP cache. This is the full path to the `resources/vep_data` directory you extracted earlier.

- The rest of the fields you can leave as they are, or modify according to the resources that you have available.

#### Local execution

The profile located in `profile/local/config.v8+.yaml` is set up to use local execution without a job scheduler. Modify the fields to adjust to your local resources.

### Annotate sample VCF file

#### SLURM

You can run snakemake with the SLURM profile as follows: