#! /bin/bash -login
#SBATCH -p med2
#SBATCH -J fastqc
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=25gb
#SBATCH --output=slurm_%x.%j.out
#SBATCH --error=slurm_%x.%j.err

# activate conda in general
. "/home/mhussien/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate astrangia

# go to a particular directory
cd /home/mhussien/astrangia/astrangia_paper/runs/fastp_new/

# make things fail on errors
set -o nounset
set -o errexit
set -x

TRIMMED_DIR=/home/mhussien/astrangia/astrangia_paper/runs/fastp_new
OUTPUT_DIR=${TRIMMED_DIR}
SAMPLES=/home/mhussien/astrangia/astrangia_paper/sample_data

R1=${SAMPLES}/Ast28B_R1_001.fastq.gz
R2=${SAMPLES}/Ast28B_R2_001.fastq.gz
SAMPLE=Ast28B

### run your commands here!
/usr/bin/time -v fastqc -f fastq -o ${OUTPUT_DIR} -t 32 --noextract ${R1} ${R2} ${TRIMMED_DIR}/*gz ${TRIMMED_DIR}/*fastq

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
