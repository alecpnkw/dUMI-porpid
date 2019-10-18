from snakemake.remote.SFTP import RemoteProvider
SFTP = RemoteProvider(username="appankow", private_key="/home/appankow/.ssh/id_rsa_themis")
configfile: "config.yaml"
paired = [(ds, config["datasets"][ds]["Run"]) for ds in config["datasets"]]
datasets = [ds for (ds,run) in paired]
runIDs = [run for (ds,run) in paired]

rule all:
    input:
        "reports/template-report.html"

rule trim_primers:
    input:
        SFTP.remote("hercules/opt/shared/PacBio_PipelineData/{runID}/Analysis/Demultiplexing/IRF3/{dataset}.fastq")
    output:
        "trimmed/{runID}/{dataset}.fastq", #add temporary() flag
        temporary("tmp/{runID}{dataset}_filt.fastq")
    params:
        fwd_primer = lambda wc: config["datasets"][wc.dataset]["Reverse_Primer_2ndRd_Sequence"],
        rev_primer = lambda wc: config["datasets"][wc.dataset]["Forward_Primer_2ndRd_Sequence"],
        target_size = 2400
    script:
        "src/trim-primers.jl"

rule porpid:
    input:
        "trimmed/{runID}/{dataset}.fastq"
    output:
        "consensus/{runID}/{dataset}.fasta",
        "tagged/{runID}/{dataset}.fastq/family_tags.csv",
        directory("tagged/{runID}/{dataset}.fastq/{dataset}"),
        directory("tagged/{runID}/{dataset}.fastq/{dataset}_keeping")
    params:
        sUMI_primer = lambda wc: config["datasets"][wc.dataset]["sUMI_Primer_Sequence"]
    script:
        "src/processing.jl"

rule qc_bins:
    input:
        "tagged/{runID}/{dataset}.fastq/family_tags.csv"
    output:
        report("postproc/{runID}/{dataset}/{dataset}_qc_bins.png", category = "PORPID QC"),
        report("postproc/{runID}/{dataset}/{dataset}_qc_bins.csv", category = "PORPID QC")
    script:
        "src/qc-bins-snakemake.jl"

rule aggregate_tags:
    input:
        expand("tagged/{runID}/{dataset}.fastq/family_tags.csv", zip, dataset = datasets, runID = runIDs)
    output:
        "reports/template-report.html"
    script:
        "src/aggregate.Rmd"
