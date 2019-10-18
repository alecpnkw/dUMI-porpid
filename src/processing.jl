using Distributed
include("../../src/functions.jl")

sUMI_primer = snakemake.params["sUMI_primer"] #fill in here.
filtered_data_file = snakemake.input[1]
template_name = split(basename(filtered_data_file),'.')[1]
outdir = dirname(dirname(snakemake.output[2]))

println("Extracting UMIs...")
t1 = time()
templates = Dict()
templates[template_name] = replace(sUMI_primer, "N" => "n")*"*"

cfg = Configuration()
cfg.files = [filtered_data_file]
cfg.filetype = fastq
cfg.start_inclusive = 0
cfg.end_inclusive = length(first(values(templates))) + 3
cfg.try_reverse_complement = false #The sequences are already oriented
for (name, template) in templates
    println("Using template $(template)")
    push!(cfg.templates, Template(name, template))
end

println("$(outdir)/$(basename(filtered_data_file))/")
if isdir("$(outdir)/$(basename(filtered_data_file))/")
    @error "The directory already exists! Please delete"
end
dir_dict = Dict()
my_output_func(source_file_name, template, tag, output_sequence, score) = my_write_to_file_count_to_dict(dir_dict, source_file_name, template, tag, output_sequence, score, outdir)
say_print_func = function(count)
    println("Processed $(count) sequences")
end
# This is the slow bit
extract_tags_from_file(cfg.files[1], cfg, my_output_func, print_every=1000, print_callback=say_print_func)
directories = collect(keys(dir_dict))

tag_dfs = []
for (ix,template) in enumerate(templates)
    template_name = basename(directories[ix])
    tag_dict = dir_dict[directories[ix]]
    delete!(tag_dict, "REJECTS")
    tag_counts = tag_dict
    tags = collect(keys(tag_dict))

    #Builts matrix of conditional probabilities.
    tag_to_index, index_to_tag = tag_index_mapping(tags)
    pacbio_error_rate = 0.005
    recurse = 1
    probabilities_array = prob_observed_tags_given_reals(tag_to_index, PORPID.PacBioErrorModel(pacbio_error_rate), recurse)
    indexed_counts = index_counts(tag_counts, tag_to_index);
    #Runs the LDA inference
    most_likely_real_for_each_obs = LDA(probabilities_array, indexed_counts);

    path = "$(outdir)/$(basename(filtered_data_file))/"*template_name*"/"

    likely_real = []
    UMI = []
    bin_sizes = []
    tags = []
    probs = []

    #Filter and copy to "_keeping"
    tag_df = filterCCSFamilies(most_likely_real_for_each_obs, path, index_to_tag, tag_counts, template_name)
    push!(tag_dfs, tag_df)
end
CSV.write("$(outdir)/$(basename(filtered_data_file))/family_tags.csv", sort!(vcat(tag_dfs...), [:Sample, :tags, :fs], rev = (false, false, true))); #io step
t2 = time()
println("UMI identification took $(t2-t1) seconds.")

## Calculate consensus sequences for each family.
println("Generating consensus...")
t1 = time()
for (template_name, template) in templates
    println("Processing $(template_name)")
    direc = template_name*".fastq"
    base_dir = "$(outdir)/"*direc*"/"*template_name*"_keeping"
    @time seq_collection, seqname_collection = generateConsensusFromDir(base_dir, template_name)
    trimmed_collection = [primer_trim(s,sUMI_primer) for s in seq_collection];
    write_fasta(snakemake.output[1],
        reverse_complement.(trimmed_collection),
        names = seqname_collection)
end
t2 = time()
println("Consensus generation took $(t2-t1) seconds.")
