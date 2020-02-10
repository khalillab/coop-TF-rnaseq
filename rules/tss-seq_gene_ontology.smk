#!/usr/bin/env python

rule gene_ontology:
    input:
        universe = config["genome"]["transcript_annotation"],
        diffexp_path = "diff_exp/transcripts/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_rnaseq-libsizenorm-transcripts-diffexp-results-{direction}.tsv",
        go_anno_path = config["genome"]["gene_ontology_mapping_file"]
    output:
        results = "gene_ontology/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-{direction}-gene-ontology-results.tsv",
        enriched_combined = "gene_ontology/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-{direction}-gene-ontology-enriched-all.svg",
        depleted_combined = "gene_ontology/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-{direction}-gene-ontology-depleted-all.svg",
        enriched_facet = "gene_ontology/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-{direction}-gene-ontology-enriched-facetted.svg",
        depleted_facet = "gene_ontology/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-{direction}-gene-ontology-depleted-facetted.svg",
    params:
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated" if wc.direction=="down" else wc.direction
    conda:
        "../envs/gene_ontology.yaml"
    script:
        "../scripts/gene_ontology.R"

