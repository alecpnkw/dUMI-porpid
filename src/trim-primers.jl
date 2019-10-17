using NextGenSeqUtils

fwd_primer = snakemake.params["fwd_primer"]
rev_primer = snakemake.params["rev_primer"]

seqs, phreds, seq_names = read_fastq(snakemake.input[1])
trimmed = [double_primer_trim(seqs[i],phreds[i],fwd_primer,rev_primer) for i in 1:length(seqs)]
write_fastq(snakemake.output[1],
    [s for (s,p) in trimmed],
    [p for (s,p) in trimmed];
    names = seq_names
    )
