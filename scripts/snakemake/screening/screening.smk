import os
import snakemake.io
import glob
import sys

def prepare_rnaSpades(param, reads):
    # print(f"original: {reads}", file=sys.stderr)
    # print(f"str: {str(reads)}", file=sys.stderr)
    reads = list(reads)
    # print(f"list: {reads}", file=sys.stderr)
    all_params = list()
    for read in reads:
        all_params.append(f"{param} {read}")
    print(' '.join(all_params), file=sys.stderr)
    return ' '.join(all_params)


ROOT_DIR = "/home/mhussien/astrangia/astrangia_paper/data/"
SAMPLES_DIR = ROOT_DIR + "samples"
TRIMMED_SAMPLES = ROOT_DIR + "trimmed"
SIGS_OUTDIR = ROOT_DIR + "sigs"
cDBG_OUTDIR = ROOT_DIR + "cDBGk75_samples"
cDBG_ASSEMBLY_DIR = ROOT_DIR + "assembled_cDBGk75_samples"
TRIMMED_ASSEMBLY_DIR = ROOT_DIR + "assembled_trimmed_samples"
ITS_DB_DIR = ROOT_DIR + "its_db"
SPLITTED_FASTA_DIR = ROOT_DIR + "splitted_fasta"
BLASTN_RESULTS_DIR = ROOT_DIR + "blastn_results"
SEQCLEAN_TRIMMED_ASSEMBLY_DIR = ROOT_DIR + "seqclean_assembled_trimmed_samples"
BUSCO_DATASET_DIR = ROOT_DIR + "BUSCO_DATASET"
BUSCO_REPORTS = ROOT_DIR + "BUSCO_REPORTS"


SAMPLES, = glob_wildcards(SAMPLES_DIR + "/{sample}_1.fastq.gz")

rule all:
    input:
        # Generate signatures for trimmed samples {R1_PE, R2_PE, Merged}
        expand("{OUTDIR}" + "/trimmed_{sample}.sig", OUTDIR = SIGS_OUTDIR, sample=SAMPLES),
        # fastp trimming, this can be turned off in favor of previous or next rules
        expand("{OUTDIR}" + "/trimmed_{sample}.fastq", OUTDIR = TRIMMED_SAMPLES, sample=SAMPLES),
        # Generate cDBG for all trimmed samples {R1_PE, R2_PE, Merged}
        expand("{OUTDIR}" + "/cDBG_k75_all_samples.{EXT}", OUTDIR = cDBG_OUTDIR, EXT = ["histo", "unitigs.fa"]),
        # Merge all signatures
        SIGS_OUTDIR + "/all_samples_k31.sig",
        # Assemble the cDBG of all trimmed samples
        cDBG_ASSEMBLY_DIR + "/transcripts.fasta",
        TRIMMED_ASSEMBLY_DIR + "/transcripts.fasta",
        # Download and create blastdb of the ITS database
        # expand(ITS_DB_DIR + "/its_8.3.{EXT}", EXT = ['nhr', 'nin', 'nsq']),
        ITS_DB_DIR + "/its_8.3.fa.ndb",
        # cDBG blastn query on ITS
        BLASTN_RESULTS_DIR + "/its_cDBG_all_samples.blastn",
        # BLASTN_RESULTS_DIR + "/its_assembled_trimmed_samples.blastn",
        # seqclean of assembled trimmed samples output
        TRIMMED_ASSEMBLY_DIR + "/cleaned_transcripts.fasta",

        BUSCO_REPORTS + "/CLEANED_TRANSCRIPT_TRIMMED/short_summary.specific.eukaryota_odb10.CLEANED_TRANSCRIPT_TRIMMED.txt",

        # BLASTn assembled trimmed samples
        # expand(BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/its_cleaned_assembled_trimmed_samples{SPLIT_PART}.blastn", SPLIT_PART = ["%03d" % (num,) for num in range(1,33,1)])
        BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/merged_its_cleaned_assembled_trimmed_samples.blastn",
        TRIMMED_ASSEMBLY_DIR + "/busco_config.ini"



# rule blast_query_its_assembled_trimmed_samples:
#     input:
#         query_fasta = TRIMMED_ASSEMBLY_DIR + "/cleaned_transcripts.fasta",
#         # one file is enough to detect that blast has created the db
#         blast_db = ITS_DB_DIR + "/its_8.3.fa.ndb"
#     output:
#         blast_results = BLASTN_RESULTS_DIR + "/its_cleaned_assembled_trimmed_samples.blastn"
#     params:
#         blastn_out_dir = BLASTN_RESULTS_DIR,
#         its_db_prefix = ITS_DB_DIR + "/its_8.3.fa"
#     threads: 32
#     resources:
#         time = 3000,
#         partition = "bmm",
#         mem_mb = 100000,
#     shell: """
#         mkdir -p {params.blastn_out_dir} && \
#         blastn -query {input.query_fasta} \
#         -db {params.its_db_prefix} \
#         -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
#         -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
#         -max_target_seqs 10  -out {output.blast_results}
#     """

rule merge_blast_query_its_assembled_trimmed_samples:
    input:
        expand(BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/its_cleaned_assembled_trimmed_samples{SPLIT_PART}.blastn", SPLIT_PART = ["%03d" % (num,) for num in range(1,33,1)])
    output:
        BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/merged_its_cleaned_assembled_trimmed_samples.blastn"
    resources:
        mem_mb = 100,
        time = 60,
        partition = "low2"
    threads: 10
    shell: """
        cat {input} > {output}
    """

rule blast_query_its_assembled_trimmed_samples:
    input:
        query_fasta = SPLITTED_FASTA_DIR + "/cleaned_transcripts/cleaned_transcripts.part_{SPLIT_PART}.fasta",
        # one file is enough to detect that blast has created the db
        blast_db = ITS_DB_DIR + "/its_8.3.fa.ndb"
    output:
        blast_results = BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/its_cleaned_assembled_trimmed_samples{SPLIT_PART}.blastn"
    params:
        blastn_out_dir = BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed",
        its_db_prefix = ITS_DB_DIR + "/its_8.3.fa"
    threads: 2
    resources:
        time = 5 * 24 * 60,
        partition = "bmm",
        mem_mb = 40 * 1024,
    shell: """
        mkdir -p {params.blastn_out_dir} && \
        blastn -query {input.query_fasta} \
        -db {params.its_db_prefix} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
        -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
        -max_target_seqs 10  -out {output.blast_results}
    """

rule split_assembled_for_blast:
    threads: 8
    input:
        TRIMMED_ASSEMBLY_DIR + "/cleaned_transcripts.fasta",
    output:
        expand(SPLITTED_FASTA_DIR + "/cleaned_transcripts/cleaned_transcripts.part_{SPLIT_PART}.fasta", SPLIT_PART = ["%03d" % (num,) for num in range(1,33,1)])
    resources:
        mem_mb = 50 * 1024,
        time = 60,
        partition = "low2"
    params:
        parts = 32,
        splitted_dir = SPLITTED_FASTA_DIR + "/cleaned_transcripts",
        splitted_fastas_dir = SPLITTED_FASTA_DIR
    shell: """
        mkdir -p {params.splitted_fastas_dir} && \
        seqkit split -p {params.parts} {input} -O {params.splitted_dir}
    """


# This needs an offline busco dataset (I preferred it that way)
# https://busco-data.ezlab.org/v5/data/
rule busco_assembly_trimmed_reads:
    threads: 32
    input:
        busco_config = TRIMMED_ASSEMBLY_DIR + "/busco_config.ini"
    output:
        BUSCO_REPORTS + "/CLEANED_TRANSCRIPT_TRIMMED/short_summary.specific.eukaryota_odb10.CLEANED_TRANSCRIPT_TRIMMED.txt"
    resources:
        partition = "high2",
        nodes = 2,
        mem_mb = 100*1024,
        time = 5 * 24 * 60
    shell: """
        /usr/bin/time -v busco --config {input.busco_config}
    """

rule prepare_config_busco_assembly_trimmed_reads:
    input: 
        TRIMMED_ASSEMBLY_DIR + "/cleaned_transcripts.fasta"
    output:
        TRIMMED_ASSEMBLY_DIR + "/busco_config.ini"
    params:
        busco_output_dir = BUSCO_REPORTS,
        busco_output = "CLEANED_TRANSCRIPT_TRIMMED",
        cores = 64,
        busco_dataset = BUSCO_DATASET_DIR
    threads: 1
    resources:
        partition = "bmh",
        time = 10,
        mem_mb = 1 * 1024
    run:
        import inspect
        config = f"""
            [busco_run]
            in = {input}
            mode = transcriptome
            out = {params.busco_output}
            out_path = {params.busco_output_dir}
            auto-lineage-euk = True
            cpu = {params.cores}
            force = False
            restart = True
            download_path = {params.busco_dataset}
            offline=True
        """
        config = inspect.cleandoc(config)
        with open(f"{output}", 'w') as OUT:
            OUT.write(config)


rule seqclean_assembly_trimmed_reads:
    input: 
        TRIMMED_ASSEMBLY_DIR + "/transcripts.fasta",
    threads: 16
    output:
        cleaned_output = protected(TRIMMED_ASSEMBLY_DIR + "/cleaned_transcripts.fasta"),
    params:
        seqclean_app = "/home/mhussien/astrangia/astrangia_paper/src/seqcleaner/seqclean",
        tmp_dir = temp(directory("/scratch/mhussien/seqclean_assembly_trimmed_reads")),
        cores = 16,
    resources:
        partition = "med2",
        time = 60*10,
        mem_mb = 50 * 1024,
        tmpdir = "/scratch/mhussien/tmp_seqclean_assembly_trimmed_reads",
    shell: """
        mkdir -p {params.tmp_dir} && cd {params.tmp_dir} && \
        {params.seqclean_app} {input} \
        -c {params.cores} -o {output.cleaned_output} && \
        rm -rf {params.tmp_dir}
    """



rule assembly_trimmed_reads:
    input:
        # Just to make it wait
        R1_PE=TRIMMED_SAMPLES + "/merged_trimmed_R1_PE.fastq",
        R2_PE=TRIMMED_SAMPLES + "/merged_trimmed_R2_PE.fastq",
        MERGED=TRIMMED_SAMPLES + "/merged_trimmed_merged.fastq",
    threads: 32
    output:
        transcripts = TRIMMED_ASSEMBLY_DIR + "/transcripts.fasta"
    params:
        cores = 32,
        spades_output_dir = TRIMMED_ASSEMBLY_DIR,
        rnaspades_tmp_dir = "/scratch/mhussien/rnaSpades_trimmed_reads_assembly/",
    resources:
        mem_mb = 200000,
        nodes = 1,
        time = 3000,
        partition = "med2"
    shell: """
        /usr/bin/time -v \
        rnaspades.py -1 {input.R1_PE} -2 {input.R2_PE} --merged {input.MERGED} -t {params.cores} \
        --tmp-dir {params.rnaspades_tmp_dir} \
        -o {params.spades_output_dir}
    """

rule prepare_assembly_trimmed_reads:
    input:
        # Just to make it wait until they all done
        merged_reads = expand("{OUTDIR}" + "/trimmed_{sample}.fastq", OUTDIR = TRIMMED_SAMPLES, sample=SAMPLES),
        R1_PE = expand(TRIMMED_SAMPLES  + "/trimmed_{sample}_R1_PE.fastq.gz", sample=SAMPLES),
        R2_PE = expand(TRIMMED_SAMPLES  + "/trimmed_{sample}_R2_PE.fastq.gz", sample=SAMPLES),
        MERGED = expand(TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz", sample=SAMPLES),
    output:
        R1_PE=TRIMMED_SAMPLES + "/merged_trimmed_R1_PE.fastq",
        R2_PE=TRIMMED_SAMPLES + "/merged_trimmed_R2_PE.fastq",
        MERGED=TRIMMED_SAMPLES + "/merged_trimmed_merged.fastq",
    shell: """
        zcat {input.R1_PE} > {output.R1_PE} && \
        zcat {input.R2_PE} > {output.R2_PE} && \
        zcat {input.MERGED} > {output.MERGED}
    """


rule blast_query_its_cDBG_allSamples:
    input:
        query_fasta = cDBG_ASSEMBLY_DIR + "/transcripts.fasta",
        # one file is enough to detect that blast has created the db
        blast_db = ITS_DB_DIR + "/its_8.3.fa.ndb"
    output:
        blast_results = BLASTN_RESULTS_DIR + "/its_cDBG_all_samples.blastn"
    params:
        blastn_out_dir = BLASTN_RESULTS_DIR,
        its_db_prefix = ITS_DB_DIR + "/its_8.3.fa"
    threads: 10
    resources:
        time = 500,
        partition = "bmm",
        mem_mb = 100000
    shell: """
        mkdir -p {params.blastn_out_dir} && \
        blastn -query {input.query_fasta} \
        -db {params.its_db_prefix} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
        -dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
        -max_target_seqs 10  -out {output.blast_results}
    """


rule create_blast_db_its:
    input: 
        its_fasta = ITS_DB_DIR + "/its_8.3.fa"
    threads: 12
    resources:
        time = 200,
        mem_mb = 10000,
        partition = "med2"
    params:
        db_name = ITS_DB_DIR + "/its_8.3.fa"
    output:
        expand(ITS_DB_DIR + "/its_8.3.fa.{EXT}", EXT = ['nhr', 'nin', 'nsq', 'ndb'])
    shell: """
         makeblastdb -in {input.its_fasta} -input_type fasta -dbtype nucl  -out {params.db_name}
    """
    

rule download_its:
    threads: 1
    resources:
        mem_mb = 500,
        partition = "med2"
    output:
        fasta_file = ITS_DB_DIR + "/its_8.3.fa",
        output_dir = protected(directory(ITS_DB_DIR))
    params:
        # https://doi.org/10.15156/BIO/1281567
        URL = "https://files.plutof.ut.ee/public/orig/28/88/28881015F8784D2A68C3F8C6CB851EAAE3803A8BFDA735C6390CE05DB4E34851.gz"

    shell: """
        mkdir -p {params.output_dir} && \
        curl {params.URL} --output its_8.3.fa.gz && \
        gunzip its_8.3.fa.gz && \
        mv its_8.3.fa {output.fasta_file}
    """

rule cDBG_assembly:
    threads: 32

    input:
        cDBG_OUTDIR + "/cDBG_k75_all_samples.unitigs.fa",
    
    output:
        transcripts = cDBG_ASSEMBLY_DIR + "/transcripts.fasta",
    
    params:
        spades_output_dir = cDBG_ASSEMBLY_DIR,
        cores = 64,
        rnaspades_tmp_dir = "/scratch/mhussien/rnaSpades_cDBG_assembly/"
    
    resources:
        mem_mb = 150000,
        nodes = 2,
        time = 2000,
        partition = "med2"
    
    shell: """
        /usr/bin/time -v \
        rnaspades.py -s {input} -t {params.cores} \
        --tmp-dir {params.rnaspades_tmp_dir} \
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