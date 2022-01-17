import os
import snakemake.io
import glob


SAMPLES_DIR="/home/mhussien/astrangia/astrangia_paper/data/samples"
TRIMMED_SAMPLES="/home/mhussien/astrangia/astrangia_paper/workflow/grist_data/trimmed"
ABUNDTRIM_SIGS="/home/mhussien/astrangia/astrangia_paper/workflow/grist_data/sigs"
ABUNDTRIM_PROTEIN_SIGS="/home/mhussien/astrangia/astrangia_paper/workflow/grist_data/aa_sigs"

SAMPLES, = glob_wildcards(SAMPLES_DIR + "/{sample}_1.fastq.gz")

rule all:
    input: expand("{OUTDIR}" + "/{sample}.abundtrim.sig", sample=SAMPLES, OUTDIR = [ABUNDTRIM_SIGS, ABUNDTRIM_PROTEIN_SIGS])

rule compute_dna_sigs:
    input:
        trimmed_sample = TRIMMED_SAMPLES + "/{sample}.trim.fq.gz"
    output:
        sample_sig = ABUNDTRIM_SIGS + "/{sample}.abundtrim.sig",
    params:
        dna_outdir = ABUNDTRIM_SIGS + "/{sample}.abundtrim.sig",
    
    shell: """
        sourmash sketch dna -p k=21,scaled=1000,abund \
           {input} -o {params.dna_outdir} \
           --name {wildcards.sample}
    """

rule compute_prot_sigs:
    input:
        trimmed_sample = TRIMMED_SAMPLES + "/{sample}.trim.fq.gz"
    output:
        sample_aa_sig = ABUNDTRIM_PROTEIN_SIGS + "/{sample}.abundtrim.sig"
    params:
        outdir = ABUNDTRIM_PROTEIN_SIGS + "/{sample}.abundtrim.sig",
    
    shell: """
        sourmash sketch translate -p k=10,scaled=100,abund \
           {input} -o {params.outdir} \
           --name {wildcards.sample}
    """

rule Trimmomatic:
    threads: 32

    input:
        r1_pe = SAMPLES_DIR + "/{sample}_1.fastq.gz",
        r2_pe = SAMPLES_DIR + "/{sample}_2.fastq.gz",
        
    params:
        OP_R1_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_SE.fastq",
        OP_R2_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_SE.fastq",
        OP_R1_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq",
        OP_R2_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq",
        adap= os.environ['CONDA_PREFIX'] + "/share/trimmomatic/adapters",
        OUTDIR = TRIMMED_SAMPLES,

    output:
        trimmed_sample = TRIMMED_SAMPLES + "/{sample}.trim.fq.gz"

    shell: """
        trimmomatic PE -threads 32 -phred33 {input.r1_pe} {input.r2_pe} {params.OP_R1_PE} {params.OP_R1_SE} {params.OP_R2_PE} {params.OP_R2_SE} ILLUMINACLIP:{params.adap}/TruSeq3-PE-2.fa:2:30:0:1 SLIDINGWINDOW:4:2 MINLEN:25
        cat {params.OP_R1_PE} {params.OP_R1_SE} {params.OP_R2_PE} | gzip - > {params.OUTDIR}/{wildcards.sample}.trim.fq.gz
    """