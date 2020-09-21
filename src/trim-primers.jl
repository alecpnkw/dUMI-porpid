using NextGenSeqUtils

fwd_primer = snakemake.params["fwd_primer"]
rev_primer = snakemake.params["rev_primer"]

#filter .fastq
filtered_path = snakemake.output[2] #julia temp paths don't work
println("Filtering .fastq file...")
@time fastq_filter(snakemake.input[1],
                   filtered_path, #path here
                   error_rate = 0.01,
                   min_length = snakemake.params["target_size"]*0.6,
                   max_length = snakemake.params["target_size"]*1.4)

seqs, phreds, seq_names = read_fastq(filtered_path)
need to reverse complement collections..
get_phred_revc = x -> [p for p in x[end:-1:1]]
seqs_revc = reverse_complement.(seqs)
phreds_revc = get_phred_revc.(phreds)
trimmed = [double_primer_trim(seqs[i],phreds[i],fwd_primer,rev_primer) for i in 1:length(seqs)]
write_fastq(snakemake.output[1],
    [s for (s,p) in trimmed],
    [p for (s,p) in trimmed];
    names = seq_names
    )
