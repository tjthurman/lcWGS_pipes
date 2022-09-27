# A generalized snakemake pipeline
# for using PCANGSD
# to do a genomic PCA.

# Takes a set of per-scaffold genotype likelihoods, previously calculated
# Does SNP thinning (position based, not LD-based) if you want.
# Combines the per-scaffold stuff.
# And then runs PCANGSD. 

import os

#################
##   GLOBALS   ##
################# 

# List the reference genome you want to map to:
REF = config["ref_genome"]


# The bamlist containing the list of bam files for input into angsd
BAM_LIST_FILE = config["bamlist"]
BAM_LIST_NAME = os.path.splitext(os.path.basename(BAM_LIST_FILE))[0]

# Batch names
# One for the caclulation of genotype likelihoods thats already been done
# And one for the output files
BATCH_NAME_GL = config["batch_name_GL"]
BATCH_NAME_PCA = config["batch_name_PCA"]

# Thinning interval
INTERVAL = config["thin_interval"]

# Path of pcangsd executable
PCANGSD_path=config["path_to_pcangsd"]

###############
##   SETUP   ##
############### 

# Get scaffold list from reference genome index
with open(REF + ".fai", "r") as f:
    lines =  f.read().splitlines()
    scaff_and_length = list(zip([line.split("\t")[0] for line in lines], [int(line.split("\t")[1]) for line in lines]))


scaffolds = [tup[0] for tup in scaff_and_length]


######################
## HELPER FUNCTIONS ##
######################
# A short function to add the bamlist and batch name to the results filepath

def bd_in(filepath):
    return os.path.normpath(os.path.join("results", "angsd_GLs", BAM_LIST_NAME, BATCH_NAME_GL, filepath))

def bd_out(filepath):
    return os.path.normpath(os.path.join("results", "pcangsd", BAM_LIST_NAME, BATCH_NAME_PCA, filepath))


####################
## PIPELINE START ##
####################

localrules: all


rule all:
    input:
        bd_out("pcangsd_out.cov")

rule thin_GLs:
    input:
        raw_GL = bd_in("scaffolds/{scaffold}_GL.beagle.gz")
    output:
        thinned_GL = bd_in("thinned_scafs/thin{interval}/{scaffold}_GL_thin.beagle.gz")
    shell:
        """
        zcat {input.raw_GL} | awk '{{if (NR == 1) {{print $0;}} else if (NR == 2) {{print $0; split($1, a, "_"); prev = a[length(a)];}} else {{split($1, b, "_"); if (b[length(b)] > prev + {wildcards.interval}) {{print $0; prev = b[length(b)];}}}}}}' | gzip > {output.thinned_GL}
        """


rule gather_per_scaff_GLS:
    input:
        expand(bd_in("thinned_scafs/thin{interval}/{scaffold}_GL_thin.beagle.gz"), scaffold = scaffolds, interval = INTERVAL)
    output:
        glf=bd_in("genome_wide_GL_filtered_thin{interval}.beagle.gz"),
        sites=bd_in("genome_wide_GL_filtered_thin{interval}.sites")
    shell:
        """
        cat {input} > {output.glf}
        zcat {output.glf} | awk -F '\t' '{{print $1}}' > {output.sites}
        """


rule pcangsd:
    input:
        expand(bd_in("genome_wide_GL_filtered_thin{interval}.beagle.gz"), interval = INTERVAL)
    output:
        bd_out("pcangsd_out.cov")
    params:
        basename=bd_out("pcangsd_out")
        exec_path=PCANGSD_path
    resources:
        cpus=config["cpus"]
    shell:
        """
        python {params.exec_path} -beagle {input} -out {params.basename} -threads {resources.cpus} -pcadapt -selection -sites_save -tree -minMaf 0.05
        """



