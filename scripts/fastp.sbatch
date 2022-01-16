#! /bin/bash -login
#SBATCH -p med2
#SBATCH -J fastp
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
cd /home/mhussien/loretta_decontamination/astrangia_paper/runs/trimming_test

# make things fail on errors
set -o nounset
set -o errexit
set -x

OUTPUT=/home/mhussien/loretta_decontamination/astrangia_paper/runs/trimming_test
SAMPLES=/home/mhussien/loretta_decontamination/astrangia_paper/sample_data

R1=${SAMPLES}/Ast28B_R1_001.fastq.gz
R2=${SAMPLES}/Ast28B_R2_001.fastq.gz

SAMPLE=Ast28B

OP_R1_SE=fastp_${SAMPLE}_R1_SE.fastq
OP_R2_SE=fastp_${SAMPLE}_R2_SE.fastq
OP_R1_PE=fastp_${SAMPLE}_R1_PE.fastq
OP_R2_PE=fastp_${SAMPLE}_R2_PE.fastq

### run your commands here!
/usr/bin/time -v fastp --in1 ${R1} --in2 ${R2} \
             --detect_adapter_for_pe  --qualified_quality_phred 4 \
             --length_required 25 --thread 32 \
             --json ${SAMPLE}.json --html ${SAMPLE}.html \
             --out1 ${OP_R1_PE} --out2 ${OP_R2_PE} \
             --unpaired1 ${OP_R1_SE} --unpaired2 ${OP_R2_SE}

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
