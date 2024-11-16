#!/usr/bin/bash -l


#SBATCH --time=24:00:00
#SBATCH --mem=10g                                                                                  
#SBATCH --cpus-per-task=1

source /home/hsiehph/shared/bin/initialize_conda.sh

conda activate snakemake_plus



snakemake -s Snakefile --use-conda --jobs 100 --profile profile -w 100  -p -k

