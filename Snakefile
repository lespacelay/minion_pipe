# Configuration file
configfile: "config.yml"

# Rule for running all
rule all:
    input:
        "virus_search/viral_contigs.fa"

# Read quality control. Trim and size select reads.
rule read_qc:
    input:
        config["reads_in"]
    output:
        "read_prep/reads_trim.fastq"
    log:
        "logs/read_prep/reads_trim.logs"
    params:
        "-q 12 -l 50"
    threads: 8
    shell:
        "porechop -t {threads} -i {input} -o read_prep/tmp.fastq && "
        "NanoFilt {params} read_prep/tmp.fastq > {output}; "
        "rm read_prep/tmp.fastq"

# DIAMOND BLASTx of reads
rule diamond_blast:
    input:
        "read_prep/reads_trim.fastq"
    output:
        "virus_search/diamond_blastx.txt"
    log:
        "logs/virus_search/diamond_blastx.log"
    shell:
        "#diamond blastx --query {input} --outfmt 6 --remote -db nr > {output}"

# Remove host reads by reference mapping
rule host_removal:
    input:
        config["host_ref"]["monkey"],
        "read_prep/reads_trim.fastq"
    output:
        "read_prep/reads_trim_hostless.fastq"
    log:
        "logs/read_prep/reads_trim_hostless.log"
    threads: 8
    shell:
        "minimap2 -ax map-ont -t {threads} {input} | "
        #"samtools view -f 4 -bS - > reads_prep/unmapped.bam && "
        "samtools fastq -f 4 -O {output} reads_prep/unmapped.bam && "
        "rm reads_prep/unmapped.bam"

# De novo assembly of reads
rule flye_assemble:
    input:
        "read_prep/reads_trim.fastq" 
#        config["reads_in"]
    output:
        "assembly/assembly.fasta",
        "assembly/assembly_graph.gfa",
        "assembly/assembly_info.txt"
    log:
        "logs/assembly/assembly.log"
    params:
        "--nano-hq"
    threads: 10
    shell:
        "flye {params} {input} --out-dir assembly --threads {threads}"

# BLASTx similarity search
rule blastx_search:
    input:
        "assembly/assembly.fa"
    output:
        "virus_search/blastx.txt"
    log:
        "logs/virus_search/blastx.log"
    shell:
        ""

# Interproscan similarity search
rule interproscan_search:
    input:
        "assembly/assembly.fa"
    output:
        "virus_search/interproscan.out"
    log:
        "logs/virus_search/interproscan.log"
    shell:
        ""

# Select virus hitting contigs
rule get_viral_contigs:
    input:
        "virus_search/interproscan.out",
        "virus_search/blastx.txt"
    output:
        "virus_search/viral_contigs.fa"
    log:
        "logs/virus_search/viral_contigs.log"
    script:
        ""
