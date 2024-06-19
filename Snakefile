# Configuration file
configfile: "config.yml"

# Rule for running all
rule all:
    input:
        "virus_search/{}_assem_blastx.txt".format(config["name"]),
        "virus_search/{}_dmd_blastx.txt".format(config["name"])

# Read quality control. Trim and size select reads.
rule read_qc:
    input:
        "data/reads/{iden}.fastq.gz"
    output:
        "read_prep/{iden}_trim.fastq"
    log:
        "logs/read_prep/{iden}_trim.logs"
    params:
        "-q 12 -l 50"
    threads: config["threads"]
    shell:
        "porechop -t {threads} -i {input} -o read_prep/tmp.fastq && "
        "NanoFilt {params} read_prep/tmp.fastq > {output}; "
        "rm read_prep/tmp.fastq"

#
rule read_assess:
    input:
        "read_prep/{iden}_trim.fastq.gz"
    output:
        "read_assess/{iden}"
    log:
        "logs/read_assess/{iden}_assess.log"
    params:
        ""
    threads: config["threads"]
    shell:
        "python longQC.py sampleqc -x ont-rapid -o read_assess/{wildcards.iden} {input}" #ont-rapid or ont-ligation

# DIAMOND BLASTx of reads
rule diamond_blast:
    input:
        "read_prep/{iden}_trim.fastq"
    output:
        "virus_search/{iden}_dmd_blastx.txt"
    log:
        "logs/virus_search/{iden}_dmd_blastx.log"
    params:
        "-d {} -e 1E-5 --max-target-seqs 1 -p 24 -b 5 -c 1 --outfmt 6 qseqid qlen sseqid stitle pident length evalue bitscore qstart qend".format(config["dmd_db"])
    shell:
        "diamond blastx -q {input} -o {output} {params}"

# DIAMOND BLASTn of reads
rule diamond_blastn:
    input:
        "read_prep/{iden}_trim.fastq"
    output:
        "virus_search/{iden}_dmd_blastn.txt"
    log:
        "logs/virus_search/{iden}_dmd_blastn.log"
    params:
        "-db {} -evalue 0.005 -max_target_seqs 5 -num_threads 24 -outfmt '6 qseqid sseqid stitle pident length evalue sscinames scomnames sskingdoms staxids'".format(config["nt_db"])
    shell:
        "blastn -query {input} -out {output} {params}"


# Remove host reads by reference mapping
rule host_removal:
    input:
        config["host_ref"]["monkey"],
        "read_prep/{iden}_trim.fastq"
    output:
        "read_prep/{iden}_trim_hostless.fastq"
    log:
        "logs/read_prep/{iden}_trim_hostless.log"
    threads: config["threads"]
    shell:
        "minimap2 -ax map-ont -t {threads} {input} | "
        "samtools view -f 4 -bS - > reads_prep/{iden}_unmapped.bam && "
        "samtools fastq -f 4 -O {output} reads_prep/{iden}_unmapped.bam && "
        "rm reads_prep/{iden}_unmapped.bam"

# De novo assembly of reads
rule flye_assemble:
    input:
        "read_prep/{iden}_trim.fastq" 
    output:
        "assembly/{iden}/assembly.fasta",
        "assembly/{iden}/assembly_graph.gfa",
        "assembly/{iden}/assembly_info.txt"
    log:
        "logs/assembly/{iden}_assembly.log"
    params:
        "--meta -i 10"
    threads: config["threads"]
    shell:
        "flye --nano-hq {input} --out-dir assembly/{wildcards.iden} --threads {threads} {params}"

# BLASTx similarity search
rule blastx_search:
    input:
        "assembly/{iden}/assembly.fasta"
    output:
        "virus_search/{iden}_assem_blastx.txt"
    log:
        "logs/virus_search/{iden}_assem_blastx.log"
    params:
        "-db {} -evalue 0.005 -max_target_seqs 5 -num_threads 24 -outfmt '6 qseqid sseqid stitle pident length evalue sscinames scomnames sskingdoms staxids'".format(config["nr_db"])
    shell:
        "blastx -out {output} -query {input} {params}"

# BLASTn similarity search
rule blastn_search:
    input:
        "assembly/{iden}/assembly.fasta"
    output:
        "virus_search/{iden}_assem_blastn.txt"
    log:
        "logs/virus_search/{iden}_assem_blastn.log"
    params:
        "-db {} -evalue 0.005 -max_target_seqs 5 -num_threads 24 -outfmt '6 qseqid sseqid stitle pident length evalue sscinames scomnames sskingdoms staxids'".format(config["nt_db"])
    shell:
        "blastn -out {output} -query {input} {params}"

# Interproscan similarity search
rule interproscan_search:
    input:
        "assembly/{iden}/assembly.fasta"
    output:
        "virus_search/{iden}_interproscan.out"
    log:
        "logs/virus_search/{iden}_interproscan.log"
    shell:
        ""

# Select virus hitting contigs
rule get_viral_contigs:
    input:
        "virus_search/{iden}_interproscan.out",
        "virus_search/{iden}_assem_blastx.txt",
        "virus_search/{iden}_assem_blastn.txt",
        "assembly/{iden}/assembly.fasta"
    output:
        "virus_search/{iden}_viral_contigs.fa"
    log:
        "logs/virus_search/{iden}_viral_contigs.log"
    script:
        ""
