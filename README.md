# lcWGS pipes

This repo contains Snakemake workflows for performing population genomics analyses from low-coverage WGS data. The pipelines assume that you already have processed your raw sequencing reads into sorted `.bam` files for each individual sample in your analysis.  

## Software installation

These pipelines require installing the following software:

1. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), for all pipelines.
2. [angsd](http://www.popgen.dk/angsd/index.php/ANGSD), for all pipelines. 
3. [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd), for the genomic PCA pipeline.

See the websites for each of these software packages for installation. Depending on your preference, you may wish to install some or all of these software packages via [conda](https://docs.conda.io/en/latest/).

Snakemake pipelines can be run locally (on laptops or desktops), but these pipelines will perform better on a high-performance computing cluster. Snakemake makes deploying pipelines on a cluster [relatively simple](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). The `SLURM_profile` folder contains an example profile for running these pipelines on a SLURM cluster, inspired by [this profile([https://github.com/jdblischak/smk-simple-slurm). You'll need to modify it to work with your cluster. 

## Current pipelines

This repo provides the following pipelines:

1. `angsd_GL_genome_wide.smk`: calculate genotype likelihoods (GLs) across the genome. This pipeline parallelizes across the scaffolds/contigs in the genome to speed computation and reduce peak memory usage. This pipeline is a prerequisite for the other two pipelines. 
2. `angsd_window_fst.smk`: calculates Fst in sliding windows across the genome. Also parallelizes across scaffolds/contigs.
3. `pcangsd_thin_and_PCA.smk`: Thins markers by position and then performs a genomic PCA on the thinned markers using PCAngsd. 

## Running pipelines

The snakemake `.smk` files describe the pipeline, and must be paired with a configuration file that specifies various options for running the pipeline (input data, filtering options, etc). There are template configuration files for each pipeline in the `sm_config` folder, which explain the various parameters that need to be specified.

After you edit the configuration files for a given pipeline, run the pipeline by invoking snakemake:

```
snakemake -s path/to/snakefile.smk \
          --configfile path/to/configfile.yml
```

When running on a SLURM cluster, use the `--profile` flag to tell snakemake the folder in which your (edited) `config.yaml` profile is saved:

```
snakemake -s path/to/snakefile.smk \
          --configfile path/to/configfile.yml \
          --profile SLURM_profile_folder/
```

In general, I recommend checking that the pipeline is doing what you're expecting by adding the `--dryrun` flag to your call to smakemake. Then, run the pipleline. You may also wish to run the pipeline with the `-k` option, to allow independent jobs to continue running if a single job fails. 
