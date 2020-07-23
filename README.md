# IRF3-porpid

### Setup

For this project, libraries are expected to contain both a UMI primer and a pair of sequencing index primers. The config file is expected to follow this formatting scheme:

```
datasets:
  <dataset_name>:
    Forward_Primer_2ndRd_Sequence: <fwd_primer_seq>
    N7_Index: N7XX
    Reverse_Primer_2ndRd_Sequence: <rev_primer_seq>
    Run: <pacbio_run_id>
    S5_Index: S5XX
    sUMI_Primer_Sequence: <umi_primer_seq>
```

The file `config.yaml` contains an example.

The pipeline is also designed to download .fastq files from our storage server hercules, which requires an ssh key for authentication. Generate your public/private key pair and edit the following lines  at the top of the `Snakefile`.

```
from snakemake.remote.SFTP import RemoteProvider
#Need to configure remote before running
SFTP = RemoteProvider(username="", private_key="")
```

`username` should be your hercules user and `private_key` should point to the location of the key on your system.

THIS PIPELINE DOES NOT PERFORM DEMULTIPLEXING, and instead expects separately demux'ed .fastq files to be in place on hercules at `/opt/shared/PacBio_PipelineData/{pacbio_run_id}/Analysis/Demultiplexing/IRF3/{dataset_name}.fastq` where "pacbio_run_id" and "dataset_name" match the config.

### Usage

Ensure all dependencies are installed, navigate to the project directory and run

```{bash}
snakemake --cores n
```

where "n" is the number of concurrent jobs you wish to run. For more on snakemake, look [here](https://snakemake.readthedocs.io/en/stable/).

### Dependencies

Julia >v1.0 must be installed on your system, which you can download from [here](https://julialang.org/downloads/). If they are not already present on your system, fire up the Julia interpreter and run the following cell to install the required Julia packages:

```{julia}
using Pkg

#MurrellGroup packages
Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="1.0", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))
Pkg.add(PackageSpec(name="DPMeansClustering", rev="1.0", url = "https://github.com/MurrellGroup/DPMeansClustering.jl.git"))
Pkg.add(PackageSpec(name="RobustAmpliconDenoising", rev="1.0", url = "https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git"))
Pkg.add(PackageSpec(name="PORPID", rev="master", url = "https://github.com/MurrellGroup/PORPID.jl.git"))

#Others
Pkg.add("BioSequences")
Pkg.add("Compose")
Pkg.add("PyCall")
Pkg.add("PyPlot")
Pkg.add("Seaborn")
Pkg.add("MultivariateStats")
Pkg.add("StatsBase")
Pkg.add("HypothesisTests")
Pkg.add("IterTools")
Pkg.add("DataFrames")
Pkg.add("DataFramesMeta")
Pkg.add("CSV")
Pkg.add("Distributions")
Pkg.add("Distances")
```

This pipeline also provides a conda environment to set up the other dependencies, including snakemake. Install conda and run

```
conda env create --file environment.yaml
conda env activate IRF3-porpid
```

Finally, your Julia executable needs to be available inside the environment. Once your environment is activated, type `which julia` to see if it is available. If it isn't, you may need to symlink it
to your environment path. On a typical system this looks something like

```
ln -s /usr/local/bin/julia ~/anaconda3/envs/IRF3-porpid/bin/
```
