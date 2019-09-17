#!/usr/bin/env python

# trim adapter sequences from 3' end of read
# if RNA-seq, also trim poly-A sequences
# (the end the poly-A is on may be different depending direction of sequencing)
# also do quality trimming, assuming Nextseq 2-color chemistry
# (sequence loss of poly-G sequences from non-Nextseq platforms should be negligible)
rule clean_reads:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    output:
        fastq = "fastq/cleaned/{sample}_rnaseq-clean.fastq.gz",
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"],
    conda:
        "../envs/cutadapt.yaml"
    threads:
        config["threads"]
    shell: """
        (cutadapt \
                --cores {threads} \
                --adapter {params.adapter} \
                --adapter "A{{100}}" \
                --times 2 \
                --cut=-1 \
                --trim-n \
                --quality-cutoff {params.trim_qual} \
                --length-tag 'length=' \
                --minimum-length 12 \
                --output {output.fastq} \
                {input}) &> {output.log}
        """

