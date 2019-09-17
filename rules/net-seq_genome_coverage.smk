#!/usr/bin/env python

# the strand of the coverage is the opposite of the aligned strand
# NOTE: '5end' coverage files generated are always the 5' end of the read.
rule genome_coverage:
    input:
        "alignment/{sample}_rnaseq-uniquemappers.bam"
    output:
        "coverage/counts/{sample}_rnaseq-{readtype}-counts-{strand}.bedgraph",
    params:
        strand = lambda wc: {"plus":  "-",
                             "minus": "+"}.get(wc.strand),
        split = lambda wc: {"5end": "",
                            "wholeread": "-split"}.get(wc.readtype),
        end = lambda wc: {"5end": "-5",
                          "wholeread": ""}.get(wc.readtype)
    wildcard_constraints:
        strand="plus|minus"
    log:
        "logs/genome_coverage/genome_coverage_{sample}-{readtype}-{strand}.log"
    shell: """
        (bedtools genomecov \
                -bga {params.end} \
                -strand {params.strand} \
                {params.split} \
                -ibam {input} | \
         LC_COLLATE=C sort \
            -k1,1 \
            -k2,2n > \
         {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_rnaseq-{readtype}-counts-{strand}.bedgraph",
        bam = "alignment/{sample}_rnaseq-uniquemappers.bam"
    output:
        normalized = "coverage/libsizenorm/{sample}_rnaseq-{readtype}-libsizenorm-{strand}.bedgraph",
    wildcard_constraints:
        strand="plus|minus"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{readtype}-{strand}.log"
    shell: """
        (awk \
            -v norm_factor=$(samtools view -c {input.bam} | \
                paste -d "" - <(echo "/1000000") | \
                bc -l) \
            'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = lambda wc: "coverage/{norm}/{sample}_rnaseq-{readtype}-{norm}-".format(**wc) + ("plus" if wc.strand=="SENSE" else "minus") +".bedgraph",
        minus = lambda wc: "coverage/{norm}/{sample}_rnaseq-{readtype}-{norm}-".format(**wc) + ("minus" if wc.strand=="SENSE" else "plus") +".bedgraph",
    output:
        "coverage/{norm}/{sample}_rnaseq-{readtype}-{norm}-{strand}.bedgraph",
    wildcard_constraints:
        strand="SENSE|ANTISENSE"
    log:
        "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{readtype}-{norm}-{strand}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}_rnaseq-{readtype}-{norm}-{strand}.bedgraph",
        fasta = lambda wc: config["genome"]["fasta"]
    output:
        "coverage/{norm}/{sample}_rnaseq-{readtype}-{norm}-{strand}.bw",
    params:
        stranded = lambda wc: [] if wc.strand in ["plus", "minus"] else """| awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2; print $1"-minus", $2}}' | LC_COLLATE=C sort -k1,1"""
    log:
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{readtype}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes {params.stranded}) {output}) &> {log}
        """

