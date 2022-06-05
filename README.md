
# RNA-seq analysis pipeline, [Bragdon *et. al.* 2022](https://doi.org/10.1101/2022.05.22.492993)

[Snakemake](https://snakemake.github.io/) workflow used to analyze RNA-seq data for the 2022 publication [*Cooperative assembly confers regulatory specificity and long-term genetic circuit stability*](https://doi.org/10.1101/2022.05.22.492993). For the pipeline with raw data, see the Zenodo archive (coming).

## workflow summary
The workflow has the following major steps:

- 3' adapter trimming, 3' quality trimming, and poly-A tail trimming with [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)
- alignment with [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml)
- selection of unique mappers
- summaries of quality statistics from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- summaries of library processing statistics
- generation of coverage tracks
- library size and spike-in normalization of coverage
- genome-wide scatterplots and correlations
- *ab initio* transcript annotation with [StringTie](https://ccb.jhu.edu/software/stringtie/)
- differential expression analysis over reference transcripts and discovered transcripts with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- data visualization (heatmaps and metagenes)

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

### required files

- FASTQ files of RNA-seq libraries.
- FASTA file of the genome
- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files:
    - transcript annotation
    - optional: other annotations for data visualization (i.e. heatmaps and metagenes)

