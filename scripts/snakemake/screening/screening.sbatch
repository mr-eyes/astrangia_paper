#! /bin/bash -login
#SBATCH -p med2
#SBATCH -J snakeJob
#SBATCH -t 07-03:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1gb
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err

# activate conda in general
. "/home/mhussien/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate astrangia

# go to a particular directory
cd /home/mhussien/astrangia/astrangia_paper/scripts/snakemake/screening


# make things fail on errors
set -o nounset
set -o errexit
set -x


snakemake -s screening.smk -j 999 --cluster-config cluster_config.json --cluster "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} -N {cluster.N} -p {cluster.partition} -c {cluster.threads} -J {cluster.name} --output {cluster.output} --error {cluster.error}"

# Print out values of the current jobs SLURM environment variables
env | grep SLURM