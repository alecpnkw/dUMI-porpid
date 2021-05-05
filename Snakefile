configfile: "config.yaml"
DATASETS = [d for d in config for s in config[d]]
TEMPLATES = [s for d in config for s in config[d]]
print(DATASETS)
print(TEMPLATES)

rule all:
    input: 
        "reports/template-report.html",
        expand("postproc/{dataset}/{template}_qc_bins.csv", zip, dataset = DATASETS, template = TEMPLATES)

rule filter_fastq:
    input:
        "fastq/{dataset}.fastq"
    output:
        "fastq/{dataset}_filt.fastq",
        #temporary("tmp/{dataset}_filt.fastq")
    params:
        #rev_primer = lambda wc: config[wc.dataset]["Reverse_Primer_2ndRd_Sequence"],
        #fwd_primer = lambda wc: config[wc.dataset]["Forward_Primer_2ndRd_Sequence"],
        min_length = 1440,
        max_length = 3360,
        error_rate = 0.01
    script:
        "src/filter.jl"

rule demux:
    input:
        "fastq/{dataset}_filt.fastq"
    output:
        directory("demux/{dataset}")
    params:
        config = lambda wc: config[wc.dataset]
    script:
        "src/demux.jl"

rule porpid:
    input:
        "demux/{dataset}"
    output:
        "consensus/{dataset}/{template}.fasta",
        "tagged/{dataset}/{template}.fastq/UMI1_family_tags.csv",
        directory("tagged/{dataset}/{template}.fastq/UMI1/"),
        directory("tagged/{dataset}/{template}.fastq/UMI1_keeping/"),
        directory("tagged/{dataset}/{template}.fastq/UMI2/"),
        directory("tagged/{dataset}/{template}.fastq/dUMI/"),
        "tagged/{dataset}/{template}.fastq/dUMI_ranked.csv"
    params:
        dir = "tagged/{dataset}",
        sUMI_primer = lambda wc: config[wc.dataset][wc.template]["sUMI_Primer_Sequence"],
        dUMI_primer = lambda wc: config[wc.dataset][wc.template]["dUMI_Primer_Sequence"],
    script:
        "src/processing.jl"

rule qc_bins:
    input:
        "tagged/{dataset}/{template}.fastq/UMI1_family_tags.csv"
    output:
        report("postproc/{dataset}/{template}_qc_bins.png", category = "PORPID QC"),
        report("postproc/{dataset}/{template}_qc_bins.csv", category = "PORPID QC")
    script:
        "src/qc-bins-snakemake.jl"

rule R_aggregate_tags:
    input:
        expand("tagged/{dataset}/{template}.fastq/UMI1_family_tags.csv", zip, dataset = DATASETS, template = TEMPLATES)
    output:
        "reports/template-report.html"
    script:
        "src/aggregate.Rmd"
