ENV["MPLBACKEND"] = "Agg"
using DataFrames, DataFramesMeta, Seaborn, CSV

tag_df = CSV.read(snakemake.input[1])
sample = snakemake.wildcards["dataset"]

"""
    family_size_umi_len_stripplot
Draws a stripplot of family sizes vs. UMI length from an input
DataFrame. Returns the figure object.
"""
function family_size_umi_len_stripplot(data)
    tight_layout()
    fig = figure(figsize = (6,2))
    ax = PyPlot.axes()

    stripplot(y = [length(ix) for ix in data[!,:UMI]],
        x = data[!,:fs],
        hue = data[!,:tags],
        alpha = 0.2, dodge = false, jitter = 0.3, orient = "h")
        labels = xlabel("fs"), ylabel("UMI len")

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # Summary for plot title
    t = @transform(data, is_likely_real = :tags .== "likely_real")
    g = groupby(t,:is_likely_real)
    cts = @based_on(g, CCS = sum(:fs))
    cts = sort!(cts, [:is_likely_real], rev = true)
    title("likely_real: $(cts[1, :CCS]) CCS, rejected: $(cts[2, :CCS]) CCS")

    return fig
end

selected = @linq tag_df |> where(:Sample .== sample)
fig = family_size_umi_len_stripplot(selected)
fig.savefig(snakemake.output[1];
    transparent = true,
    dpi = 200,
    bbox_inches = "tight")
summary = @linq selected |>
    groupby(:tags) |>
    based_on(n=length(:fs),CCS=sum(:fs))
CSV.write(snakemake.output[2],summary)
