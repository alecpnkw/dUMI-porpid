from snakemake.remote.SFTP import RemoteProvider
SFTP = RemoteProvider(username="appankow", private_key="/home/appankow/.ssh/id_rsa_themis")
configfile: "config.yaml"
runIDs = [config["datasets"][ds]["Run"] for ds in config["datasets"]]

rule all:
    input:
        expand("consensus/{runID}/{dataset}.fasta", dataset = config["datasets"], runID = runIDs)

rule trim_primers:
    input:
        SFTP.remote("hercules/opt/shared/PacBio_PipelineData/2019-03_NIH_B670-IRF3_M4363/Analysis/Demultiplexing/IRF3/{dataset}.fastq")
    output:
        "trimmed/{runID}/{dataset}.fastq" #add temporary() flag
    params:
        fwd_primer = lambda wc: config["datasets"][wc.dataset]["Reverse_Primer_2ndRd_Sequence"],
        rev_primer = lambda wc: config["datasets"][wc.dataset]["Forward_Primer_2ndRd_Sequence"]
    script:
        "src/trim-primers.jl"

rule porpid:
    input:
        "trimmed/{runID}/{dataset}.fastq"
    output:
        "consensus/{runID}/{dataset}.fasta",
        directory("tagged/{runID}/{dataset}.fastq")
    params:
        sUMI_primer = lambda wc: config["datasets"][wc.dataset]["sUMI_Primer_Sequence"]
    script:
        "src/processing.jl"
