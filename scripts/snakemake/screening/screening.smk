import os
import snakemake.io
import glob


SAMPLES_DIR="/home/mhussien/astrangia/astrangia_paper/data/samples"
TRIMMED_SAMPLES="/home/mhussien/astrangia/astrangia_paper/data/trimmed"
SIGS_OUTDIR="/home/mhussien/astrangia/astrangia_paper/data/sigs"
cDBG_OUTDIR="/home/mhussien/astrangia/astrangia_paper/data/cDBGk75_samples"
cDBG_ASSEMBLY_DIR="/home/mhussien/astrangia/astrangia_paper/data/assembled_cDBGk75_samples"

SAMPLES, = glob_wildcards(SAMPLES_DIR + "/{sample}_1.fastq.gz")

rule all:
    input:
        expand("{OUTDIR}" + "/trimmed_{sample}.sig", OUTDIR = SIGS_OUTDIR, sample=SAMPLES),
        expand("{OUTDIR}" + "/trimmed_{sample}.fastq", OUTDIR = TRIMMED_SAMPLES, sample=SAMPLES),
        expand("{OUTDIR}" + "/cDBG_k75_all_samples.{EXT}", OUTDIR = cDBG_OUTDIR, EXT = ["histo", "unitigs.fa"]),
        SIGS_OUTDIR + "/all_samples_k31.sig",
        cDBG_ASSEMBLY_DIR + "/transcripts.fasta",

rule cDBG_assembly:
    threads: 32

    input:
        cDBG_OUTDIR + "/cDBG_k75_all_samples.unitigs.fa",
    
    output:
        transcripts = cDBG_ASSEMBLY_DIR + "/transcripts.fasta",
        rnaspades_tmp_dir = temp(directory("/scratch/mhussien/rnaSpades_cDBG_assembly/"))
    
    params:
        spades_output_dir = cDBG_ASSEMBLY_DIR,
        cores = 32
    
    resources:
        mem_mb = 300000,
        nodes = 1,
        time = 2000,
        partition = "bmm"
    
    shell: """
        /usr/bin/time -v \
        rnaspades.py -s {input} -t {params.cores} \
        --tmp-dir {output.rnaspades_tmp_dir} \
        -o {params.spades_output_dir}
    """

# might need some enhacements for handling tmp files on scratch
rule bcalm:
    input: expand("{OUTDIR}" + "/trimmed_{sample}.fastq", OUTDIR = TRIMMED_SAMPLES, sample=SAMPLES),
    threads: 16
    resources:
        time = 1000,
        mem_mb = 200000,
        nodes=2,
        partition="bmh",
    params:
        bcalm_app = "/home/mhussien/astrangia/astrangia_paper/src/bcalm/build/bcalm",
        bcalm_tmp_dir = "/scratch/mhussien/bcalm_k75/",
        k = 75,
        max_ram = 180000,
        cores = 32,
        min_abund = 3,
        out_dir = cDBG_OUTDIR,
        out_prefix = "cDBG_k75_all_samples"
    output:
        tmp_dir = temp(directory("/scratch/mhussien/bcalm_k75/")),
        histo = cDBG_OUTDIR + "/cDBG_k75_all_samples.histo",
        unitigs = cDBG_OUTDIR + "/cDBG_k75_all_samples.unitigs.fa",
        bcalm_list = temp("bcalm.list")
    shell: """
        ls -d -1 {input} > {output.bcalm_list} && \
        mkdir -p {output.tmp_dir} && \
        {params.bcalm_app} -kmer-size {params.k} -nb-cores {params.cores} \
        -max-memory {params.max_ram} -abundance-min {params.min_abund} \
        -out-tmp {output.tmp_dir} -out-dir {params.out_dir} \
        -in {output.bcalm_list} \
        -out {params.out_prefix} -verbose 1  -histo 1 && \
        mv {params.out_prefix}.unitigs {params.out_prefix}.histo {params.out_dir} 
    """


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