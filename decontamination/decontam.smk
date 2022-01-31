import os
import snakemake.io
import glob
import sys

ROOT_DIR = "/home/mhussien/astrangia/astrangia_paper/decontamination"
CONTAMINANTS_DIR = ROOT_DIR + "/contaminants"
CONTAMINANTS_SKETCHES = ROOT_DIR + "/genomes_sketches"
CONTAMINANTS_IDX = ROOT_DIR + "/index"
DECONTAMINATED_PARTITIONS = ROOT_DIR + "/genomes_partitions"
REF_FASTA_DIR = ROOT_DIR + "/ref"
BUSCO_DATASET_DIR = "/home/mhussien/astrangia/astrangia_paper/data/BUSCO_DATASET"
BUSCO_REPORTS = ROOT_DIR + "/BUSCO_REPORTS"

GENOMES, = glob_wildcards(CONTAMINANTS_DIR + "/{genome}")

rule all:
    input:
        expand("{OUTDIR}" + "/{genome}.{ext}", OUTDIR = CONTAMINANTS_SKETCHES, genome=GENOMES, ext = ["phmap", "extra"]),
        idx = CONTAMINANTS_IDX + "/genomes_sketches.phmap",
        namesMap = CONTAMINANTS_IDX + "/genomes_sketches.namesMap",
        unmapped_partition = DECONTAMINATED_PARTITIONS + "/cleaned_assembled_trimmed_samples_genomes_sketches" + "/unmapped_partition.fa",
        decontaminated_astrangia_busco_config = BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA_busco_config.ini",
        decontaminated_astrangia_busco_result = BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA/short_summary.specific.eukaryota_odb10.DECONTAMINATED_ASTRANGIA.txt"

        decontaminated_astrangia_with_unmapped_busco_config = BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA_with_unmapped_busco_config.ini",
        decontaminated_astrangia_with_unmapped_busco_result = BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA_WITH_UNMAPPED/short_summary.specific.eukaryota_odb10.DECONTAMINATED_ASTRANGIA_WITH_UNMAPPED.txt"


rule prepare_config_busco_decontaminated_astrangia_with_unmapped:
    input:
        astrangia_decontaminated = DECONTAMINATED_PARTITIONS + "/cleaned_assembled_trimmed_samples_genomes_sketches" + "/unmapped_with_genome_apoculata.assembly.scaffolds_chromosome_level.fa"
    output:
        BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA_with_unmapped_busco_config.ini"
    params:
        busco_output_dir = BUSCO_REPORTS,
        busco_output = "DECONTAMINATED_ASTRANGIA_WITH_UNMAPPED",
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


# This needs an offline busco dataset (I preferred it that way)
# https://busco-data.ezlab.org/v5/data/
rule busco_decontaminated_astrangia:
    threads: 32
    input:
        busco_config = BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA_busco_config.ini"
    output:
        BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA/short_summary.specific.eukaryota_odb10.DECONTAMINATED_ASTRANGIA.txt"
    resources:
        partition = "high2",
        nodes = 2,
        mem_mb = 100*1024,
        time = 5 * 24 * 60
    shell: """
        /usr/bin/time -v busco --config {input.busco_config}
    """

rule prepare_config_busco_decontaminated_astrangia:
    input:
        astrangia_decontaminated = DECONTAMINATED_PARTITIONS + "/cleaned_assembled_trimmed_samples_genomes_sketches" + "/genome_apoculata.assembly.scaffolds_chromosome_level.fa_partition.fa"
    output:
        BUSCO_REPORTS + "/DECONTAMINATED_ASTRANGIA_busco_config.ini"
    params:
        busco_output_dir = BUSCO_REPORTS,
        busco_output = "DECONTAMINATED_ASTRANGIA",
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


rule decontaminate_astragnia:
    input:
        idx = CONTAMINANTS_IDX + "/genomes_sketches.phmap",
        ref_fasta = REF_FASTA_DIR + "/cleaned_assembled_trimmed_samples.fasta",
    params:
        idx_prefix_path = CONTAMINANTS_IDX + "/genomes_sketches",
        ref_name = "cleaned_assembled_trimmed_samples",
        idx_prefix = "genomes_sketches",
        out_dir = DECONTAMINATED_PARTITIONS,
        ref_decontm = "/home/mhussien/astrangia/astrangia_paper/src/refDecontam/build/decontaminate",
    output:
        unmapped_partition = DECONTAMINATED_PARTITIONS + "/cleaned_assembled_trimmed_samples_genomes_sketches" + "/unmapped_partition.fa"
    resources:
        mem_mb = 500 * 1024,
        time = 3 * 24 * 60,
        partition = "bmm",
    shell: """
        mkdir -p {params.out_dir}/{params.ref_name}_{params.idx_prefix} && \
        cd {params.out_dir}/{params.ref_name}_{params.idx_prefix} && \
        /usr/bin/time -v {params.ref_decontm} {input.ref_fasta} {params.idx_prefix_path}
    """
    

rule index_cotaminants_sketches:
    threads: 1
    input:
        expand("{OUTDIR}" + "/{genome}.{ext}", OUTDIR = CONTAMINANTS_SKETCHES, genome=GENOMES, ext = ["phmap", "extra"]),
    params:
        sketches_dir = CONTAMINANTS_SKETCHES,
        idx_dir = CONTAMINANTS_IDX,
        idx_prefix = "genomes_sketches"
    resources:
        mem_mb = 800 * 1024,
        time = 3 * 24 * 60,
        partition = "bmm",
    output:
        idx = CONTAMINANTS_IDX + "/genomes_sketches.phmap",
        namesMap = CONTAMINANTS_IDX + "/genomes_sketches.namesMap"
    shell: """
        mkdir -p {params.idx_dir} && \
        cd {params.idx_dir} && \
        /usr/bin/time -v kSpider index_datasets --dir {params.sketches_dir}
    """


rule sketch_contaminants:
    threads: 1
    resources: 
        mem_mb=50000, 
        time=2 * 24 * 60,
        partition="med2"
    input:
        genome_fasta = CONTAMINANTS_DIR + "/{genome}"
    output:
        kframe = CONTAMINANTS_SKETCHES + "/{genome}.phmap",
        extra = CONTAMINANTS_SKETCHES + "/{genome}.extra"
    params:
        working_dir = CONTAMINANTS_SKETCHES,
        kSize = 21,

    shell: """
        mkdir -p {params.working_dir} && \
        cd {params.working_dir} && \
        kSpider sketch --fastx {input.genome_fasta} -k {params.kSize}
    """
