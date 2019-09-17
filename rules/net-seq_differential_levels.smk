#!/usr/bin/env python

localrules:
    map_counts_to_annotation,
    combine_annotation_counts

rule map_counts_to_annotation:
    input:
        bed = lambda wc: "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_merged-transcripts-annotated.bed" if wc.annotation=="transcripts" else config["differential_expression"]["annotations"][wc.annotation],
        bg = "coverage/counts/{sample}_rnaseq-5end-counts-SENSE.bedgraph"
    output:
        temp("diff_exp/{annotation}/{condition}-v-{control}/{sample}_rnaseq-counts-{annotation}.tsv")
    log:
        "logs/map_counts_to_annotation/map_counts_to_annotation_{condition}-v-{control}_{sample}-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input.bed} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule combine_annotation_counts:
    input:
        lambda wc: ["diff_exp/{annotation}/{condition}-v-{control}/" + x + "_rnaseq-counts-{annotation}.tsv" for x in get_samples("passing", [wc.control, wc.condition])]
    output:
        "diff_exp/{annotation}/{condition}-v-{control}/{condition}-v-{control}_rnaseq-counts-{annotation}.tsv"
    params:
        n = lambda wc: 7*len(get_samples("passing", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", [wc.control, wc.condition]))
    log:
        "logs/combine_annotation_counts/combine_annotation_counts-{condition}-v-{control}-{annotation}.log"
    shell: """
        (paste {input} | \
         cut -f$(paste -d, <(echo "1-6") <(seq -s, 7 7 {params.n})) | \
         cat <(echo -e "chrom\tstart\tend\tname\tscore\tstrand\t{params.names}") - > {output}) &> {log}
        """

rule differential_expression:
    input:
        exp_counts = "diff_exp/{annotation}/{condition}-v-{control}/{condition}-v-{control}_rnaseq-counts-{annotation}.tsv",
    output:
        results_all = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-all.tsv",
        bed_all = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-all.bed",
        results_up = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-up.tsv",
        bed_up = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-up.bed",
        results_down = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-down.tsv",
        bed_down = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-down.bed",
        results_unchanged = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-unchanged.tsv",
        bed_unchanged = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-diffexp-results-unchanged.bed",
        counts_norm = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-counts-sizefactornorm.tsv",
        counts_rlog = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-counts-rlogtransformed.tsv",
        qc_plots = "diff_exp/{annotation}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{annotation}-qc-plots.svg"
    params:
        samples = lambda wc: get_samples("passing", [wc.control, wc.condition]),
        groups = lambda wc: [PASSING[x]["group"] for x in get_samples("passing", [wc.control, wc.condition])],
        alpha = config["differential_expression"]["fdr"],
        lfc = log2(config["differential_expression"]["fold-change-threshold"])
    conda:
        "../envs/diff_exp.yaml"
    script:
        "../scripts/differential_expression_netseq.R"

# rule summarise_diffexp_results:
#     input:
#         total = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-all.tsv",
#         genic = "diff_exp/transcripts/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-genic-all.tsv",
#         antisense = "diff_exp/transcripts/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-antisense-all.tsv",
#         convergent = "diff_exp/transcripts/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-convergent-all.tsv",
#         divergent = "diff_exp/transcripts/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-divergent-all.tsv",
#         intergenic = "diff_exp/transcripts/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-intergenic-all.tsv",
#     output:
#         summary_table = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-summary.tsv",
#         mosaic = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-mosaic.svg",
#         maplot = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-maplot.svg",
#         volcano = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-volcano.svg",
#         volcano_free = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-volcano-freescale.svg",
#     params:
#         lfc = config["differential_expression"]["fold-change-threshold"],
#         alpha = config["differential_expression"]["fdr"]
#     conda:
#         "../envs/tidyverse.yaml"
#     script:
#         "../scripts/plot_diffexp_summary.R"

