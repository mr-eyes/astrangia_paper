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
cd /home/mhussien/loretta_decontamination/astrangia_paper/runs/fastqc

# make things fail on errors
set -o nounset
set -o errexit
set -x

TRIMMED_DIR=/home/mhussien/loretta_decontamination/astrangia_paper/runs/trimming_test
OUTPUT_DIR=/home/mhussien/loretta_decontamination/astrangia_paper/runs/fastqc

SAMPLES=/home/mhussien/loretta_decontamination/astrangia_paper/sample_data

R1=${SAMPLES}/Ast28B_R1_001.fastq.gz
R2=${SAMPLES}/Ast28B_R2_001.fastq.gz
SAMPLE=Ast28B

fastp_OP_R1_SE=${TRIMMED_DIR}/fastp_${SAMPLE}_R1_SE.fastq
fastp_OP_R2_SE=${TRIMMED_DIR}/fastp_${SAMPLE}_R2_SE.fastq
fastp_OP_R1_PE=${TRIMMED_DIR}/fastp_${SAMPLE}_R1_PE.fastq
fastp_OP_R2_PE=${TRIMMED_DIR}/fastp_${SAMPLE}_R2_PE.fastq

trim_OP_R1_SE=${TRIMMED_DIR}/trimmomatic_${SAMPLE}_R1_SE.fastq
trim_OP_R2_SE=${TRIMMED_DIR}/trimmomatic_${SAMPLE}_R2_SE.fastq
trim_OP_R1_PE=${TRIMMED_DIR}/trimmomatic_${SAMPLE}_R1_PE.fastq
trim_OP_R2_PE=${TRIMMED_DIR}/trimmomatic_${SAMPLE}_R2_PE.fastq

### run your commands here!
/usr/bin/time -v fastqc -f fastq -o ${OUTPUT_DIR} -t 32 --noextract ${R1} ${R2} ${fastp_OP_R1_SE} ${fastp_OP_R2_SE} ${fastp_OP_R1_PE} ${fastp_OP_R2_PE} ${trim_OP_R1_SE} ${trim_OP_R2_SE} ${trim_OP_R1_PE} ${trim_OP_R2_PE}

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
