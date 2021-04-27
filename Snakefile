configfile: "2020-09_UW_M016.yaml"
datasets = config["datasets"]

rule all:
    input:
        "reports/template-report.html",
        expand("postproc/{dataset}/{dataset}_qc_bins.csv", dataset = datasets)

rule filter_and_trim_primers:
    input:
        "fastq/{dataset}.fastq"
    output:
        "trimmed/{dataset}.fastq",
        temporary("tmp/{dataset}_filt.fastq")
    params:
        rev_primer = lambda wc: config["datasets"][wc.dataset]["Reverse_Primer_2ndRd_Sequence"],
        fwd_primer = lambda wc: config["datasets"][wc.dataset]["Forward_Primer_2ndRd_Sequence"],
        min_length = 1440,
        max_length = 3360,
        error_rate = 0.01
    script:
        "src/filter-trim-primers.jl"

rule porpid:
    input:
        "trimmed/{dataset}.fastq"
    output:
        "consensus/{dataset}.fasta",
        "tagged/{dataset}.fastq/UMI1_family_tags.csv",
        directory("tagged/{dataset}.fastq/UMI1/"),
        directory("tagged/{dataset}.fastq/UMI1_keeping/"),
        directory("tagged/{dataset}.fastq/UMI2/"),
        directory("tagged/{dataset}.fastq/dUMI/"),
        "tagged/{dataset}.fastq/dUMI_ranked.csv"
    params:
        dir = "tagged/",
        sUMI_primer = lambda wc: config["datasets"][wc.dataset]["sUMI_Primer_Sequence"],
        dUMI_primer = lambda wc: config["datasets"][wc.dataset]["dUMI_Primer_Sequence"],
    script:
        "src/processing.jl"

rule qc_bins:
    input:
        "tagged/{dataset}.fastq/UMI1_family_tags.csv"
    output:
        report("postproc/{dataset}/{dataset}_qc_bins.png", category = "PORPID QC"),
        report("postproc/{dataset}/{dataset}_qc_bins.csv", category = "PORPID QC")
    script:
        "src/qc-bins-snakemake.jl"

rule R_aggregate_tags:
    input:
        expand("tagged/{dataset}.fastq/UMI1_family_tags.csv", zip, dataset = datasets)
    output:
        "reports/template-report.html"
    script:
        "src/aggregate.Rmd"
