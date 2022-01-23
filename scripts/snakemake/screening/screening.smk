import os
import snakemake.io
import glob


SAMPLES_DIR="/home/mhussien/astrangia/astrangia_paper/data/samples"
TRIMMED_SAMPLES="/home/mhussien/astrangia/astrangia_paper/data/trimmed"
SIGS_OUTDIR="/home/mhussien/astrangia/astrangia_paper/data/sigs"

SAMPLES, = glob_wildcards(SAMPLES_DIR + "/{sample}_1.fastq.gz")

rule all:
    input:
        expand("{OUTDIR}" + "/trimmed_{sample}.sig", OUTDIR = SIGS_OUTDIR, sample=SAMPLES),
        SIGS_OUTDIR + "/all_samples_k31.sig",


rule merge_signatures:
    input: expand("{OUTDIR}" + "/trimmed_{sample}.sig", OUTDIR = SIGS_OUTDIR, sample=SAMPLES),
    resources:
        time = 120
    output:
        merged_sig = SIGS_OUTDIR + "/all_samples_k31.sig"
    shell: """
        sourmash signature merge {input} -o {output.merged_sig}
    """

rule compute_signatures:
    threads: 32

    input:
        MERGED_SAMPLE = TRIMMED_SAMPLES + "/trimmed_{sample}.fastq",

    output:
        SIG = SIGS_OUTDIR + "/trimmed_{sample}.sig",

    shell: """
        sourmash sketch dna -p k=31,scaled=1000,abund \
           {input.MERGED_SAMPLE} -o {output.SIG} \
           --name {wildcards.sample}
    """

rule merge_fastp_output:
    threads: 32

    input:
        OP_R1_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_SE.fastq.gz",
        OP_R2_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_SE.fastq.gz",
        OP_R1_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq.gz",
        OP_R2_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq.gz",
        OP_FAILED=TRIMMED_SAMPLES + "/trimmed_{sample}_failed.fastq.gz",
        OP_MERGED=TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz",

    output:
        MERGED_SAMPLE = TRIMMED_SAMPLES + "/trimmed_{sample}.fastq",

    shell: """
        zcat {input.OP_R1_PE} {input.OP_R2_PE} {input.OP_MERGED} > {output.MERGED_SAMPLE}
    """


rule fastp:
    threads: 32
    resources: 
        mem_mb=30000, 
        time_min=60, 
    input:
        r1_pe = SAMPLES_DIR + "/{sample}_1.fastq.gz",
        r2_pe = SAMPLES_DIR + "/{sample}_2.fastq.gz",

    output:
        OP_R1_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_SE.fastq.gz",
        OP_R2_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_SE.fastq.gz",
        OP_R1_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq.gz",
        OP_R2_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq.gz",
        OP_FAILED=TRIMMED_SAMPLES + "/trimmed_{sample}_failed.fastq.gz",
        OP_MERGED=TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz",

    shell: """
        fastp --in1 {input.r1_pe} --in2 {input.r2_pe} \
        --cut_front --cut_right --cut_window_size 4 \
        --cut_mean_quality 5 --trim_poly_x --length_required 25 \
        --low_complexity_filter --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --merge \
        --merged_out {output.OP_MERGED} --out1 {output.OP_R1_PE} \
        --out2 {output.OP_R2_PE} --unpaired1 {output.OP_R1_SE} \
        --unpaired2 {output.OP_R2_SE} --failed_out {output.OP_FAILED}
    """