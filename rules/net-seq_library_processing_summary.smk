#!/usr/bin/env python

localrules:
    aggregate_read_numbers,

rule aggregate_read_numbers:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("alignment/{sample}/align_summary.txt", sample=SAMPLES),
        # nodups = expand("alignment/{sample}_rnaseq-uniquemappers.bam", sample=SAMPLES)
    output:
        "qual_ctrl/read_processing/rnaseq_read-processing-summary.tsv"
    log:
        "logs/aggregate_read_numbers.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map" > {output}) &> {log}""")
        for sample, adapter, align in zip(SAMPLES.keys(), input.adapter, input.align):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(awk 'NR==3 {{ORS="\t"; print $3}} NR==4 {{ORS="\\n"; print $3}}' {align} >> {output}) &> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$5=$4-$5; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing/rnaseq_read-processing-summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing/rnaseq_read-processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing/rnaseq_read-processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing/rnaseq_read-processing-loss.svg",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/processing_summary.R"

