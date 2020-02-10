#!/usr/bin/env python

import os
from math import log2
import itertools

configfile: "config.yaml"

SAMPLES = config["samples"]
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"]]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"]]))

COUNTTYPES = ["counts"]
NORMS = ["libsizenorm"]

FIGURES = config["figures"]

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys())),
    group = "|".join(set(re.escape(v["group"]) for k,v in SAMPLES.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + ["all"])),
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys())),
    annotation = "|".join(re.escape(x) for x in set(list(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES])) + list(config["differential_expression"]["annotations"].keys() if config["differential_expression"]["annotations"] else []) + ["transcripts"])),
    status = "all|passing",
    counttype= "counts",
    norm = "counts|libsizenorm",
    readtype = "5end|wholeread",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",

status_norm_sample_dict = {
    "all": SAMPLES,
    "passing": PASSING}

def get_samples(status, groups):
    if "all" in groups:
        return(list(status_norm_sample_dict[status].keys()))
    else:
        return([k for k,v in status_norm_sample_dict[status].items() if v["group"] in groups])

def cluster_samples(status, cluster_groups, cluster_strands):
    ll = []
    for group, strand in zip(cluster_groups, cluster_strands):
        sublist = [k for k,v in status_norm_sample_dict[status].items() if v["group"] in cluster_groups]
        if strand in ["sense", "both"]:
            ll.append([f"{sample}-sense" for sample in sublist])
        if strand in ["antisense", "both"]:
            ll.append([f"{sample}-antisense" for sample in sublist])
    return(list(itertools.chain(*ll)))

include: "rules/net-seq_clean_reads.smk"
include: "rules/net-seq_alignment.smk"
include: "rules/net-seq_genome_coverage.smk"
include: "rules/net-seq_fastqc.smk"
include: "rules/net-seq_library_processing_summary.smk"
include: "rules/net-seq_sample_similarity.smk"
include: "rules/net-seq_datavis.smk"
include: "rules/net-seq_differential_levels.smk"
include: "rules/net-seq_transcript_annotation.smk"
include: "rules/tss-seq_gene_ontology.smk"
# include: "rules/net-seq_transcript_classification.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #FastQC
        'qual_ctrl/fastqc/rnaseq-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_rnaseq-uniquemappers.bam",
                sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_rnaseq-{readtype}-{norm}-{strand}.bw",
                norm=["counts","libsizenorm"],
                sample=SAMPLES,
                readtype=["5end", "wholeread"],
                strand=["SENSE", "ANTISENSE", "plus", "minus"]),
        #quality control
        "qual_ctrl/read_processing/rnaseq_read-processing-loss.svg",
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_rnaseq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg",
            zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)),
            status=statuscheck(SAMPLES, PASSING),
            windowsize=config["scatterplot_binsizes"]),
        #datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/rnaseq-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bygroup-sense.svg",
            zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)),
            figure=FIGURES,
            status=statuscheck(SAMPLES, PASSING),
            readtype=["5end", "wholeread"]) if config["plot_figures"] else [],
        #differential expression
        expand(expand("diff_exp/transcripts/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-transcripts-diffexp-results-{{direction}}.tsv",
            zip, condition=conditiongroups, control=controlgroups),
            direction=["all", "up", "down", "unchanged"]),
        expand(expand("diff_exp/{{annotation}}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-{{annotation}}-diffexp-results-{{direction}}.tsv",
            zip, condition=conditiongroups, control=controlgroups),
            direction=["all", "up", "down", "unchanged"],
            annotation=list(config["differential_expression"]["annotations"].keys()) if config["differential_expression"]["annotations"] else []),
        #gene ontology
        expand(expand("gene_ontology/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-{{direction}}-gene-ontology-enriched-all.svg",
            zip, condition=conditiongroups, control=controlgroups),
            direction=["up", "down", "unchanged"]),


