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
cd /home/mhussien/astrangia/astrangia_paper/runs/fastp_new

# make things fail on errors
set -o nounset
set -o errexit
set -x

OUTPUT=/home/mhussien/astrangia/astrangia_paper/runs/fastp_new
SAMPLES=/home/mhussien/astrangia/astrangia_paper/sample_data

R1=${SAMPLES}/Ast28B_R1_001.fastq.gz
R2=${SAMPLES}/Ast28B_R2_001.fastq.gz

sample=Ast28B

fastp --in1 ${R1} --in2 ${R2} --cut_front --cut_right --cut_window_size 4 --cut_mean_quality 5 --trim_poly_x --length_required 25 --low_complexity_filter --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --merge --merged_out ${sample}.merged.fq.gz --out1 ${sample}.R1.pe.fq.gz --out2 ${sample}.R2.pe.fq.gz --unpaired1 ${sample}.R1.se.fq.gz --unpaired2 ${sample}.R2.se.fq.gz --failed_out ${sample}.failed.fq.gz

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
