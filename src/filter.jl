using NextGenSeqUtils

#fwd_primer = snakemake.params["fwd_primer"]
#rev_primer = snakemake.params["rev_primer"]

#filter .fastq
filtered_path = snakemake.output[1] #julia temp paths don't work
println("Filtering .fastq file...")
@time fastq_filter(snakemake.input[1],
                   filtered_path, #path here
                   error_rate = snakemake.params["error_rate"],
                   min_length = snakemake.params["min_length"],
                   max_length = snakemake.params["max_length"])
 
#=
seqs, phreds, seq_names = read_fastq(snakemake.input[1])
#need to reverse complement collections...
get_phred_revc = x -> [p for p in x[end:-1:1]]
seqs_revc = reverse_complement.(seqs)
phreds_revc = get_phred_revc.(phreds)
trimmed = [double_primer_trim(seqs_revc[i],phreds_revc[i],fwd_primer,rev_primer) for i in 1:length(seqs)]
write_fastq(snakemake.output[1],
    [s for (s,p) in trimmed],
    [p for (s,p) in trimmed];
    names = seq_names
    )
=#