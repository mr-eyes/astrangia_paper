#! /bin/bash -login
#SBATCH -p med2
#SBATCH -J trmmtk
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

OP_R1_SE=trimmomatic_${SAMPLE}_R1_SE.fastq
OP_R2_SE=trimmomatic_${SAMPLE}_R2_SE.fastq
OP_R1_PE=trimmomatic_${SAMPLE}_R1_PE.fastq
OP_R2_PE=trimmomatic_${SAMPLE}_R2_PE.fastq
adap="$CONDA_PREFIX/share/trimmomatic/adapters"

CMD="trimmomatic PE -threads 32 -phred33 ${R1} ${R2} ${OP_R1_PE} ${OP_R1_SE} ${OP_R2_PE} ${OP_R2_SE} ILLUMINACLIP:${adap}/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:25"

### run your commands here!
/usr/bin/time -v ${CMD}

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
