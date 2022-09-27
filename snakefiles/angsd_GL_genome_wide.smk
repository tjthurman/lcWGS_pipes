# A generalized snakemake pipeline
# which calculates genotype likelihoods genome-wide using ANGSD
# for a given bamlist. 

import os

#################
##   GLOBALS   ##
################# 

# List the reference genome you want to map to:
REF = config["ref_genome"]

# The bamlist containing the list of bam files for input into angsd
BAM_LIST_FILE = config["bamlist"]
BAM_LIST_NAME = os.path.splitext(os.path.basename(BAM_LIST_FILE))[0]

# Batch name
BATCH_NAME = config["batch_name"]

# do we gather into one genome-wide file?
GATHER = config["gather"]

# Are we testing the pipeline or not?
TEST = config["test"]

###############
##   SETUP   ##
############### 

# Get scaffold list from reference genome index
with open(REF + ".fai", "r") as f:
    lines =  f.read().splitlines()
    scaff_and_length = list(zip([line.split("\t")[0] for line in lines], [int(line.split("\t")[1]) for line in lines]))

# For testing on a small number of scafs, if you want to make sure its working
if TEST:
    sorted_by_length = sorted(scaff_and_length, key=lambda tup: tup[1])
    scaffolds = [tup[0] for tup in sorted_by_length[0:3]]
else:
    scaffolds = [tup[0] for tup in scaff_and_length]

######################
## HELPER FUNCTIONS ##
######################
# A short function to add the bamlist and batch name to the results filepath
def bd(filepath):
    return os.path.normpath(os.path.join("results", "angsd_GLs", BAM_LIST_NAME, BATCH_NAME, filepath))

####################
## PIPELINE START ##
####################

onsuccess:
    smk_copy_command = 'cp snakefiles/angsd_GL_genome_wide.smk ' + str(bd("snakefile_used_copy.smk"))
    os.popen(smk_copy_command) 

# After rule "all", pipeline order goes from top to bottom
localrules: all

if TEST or GATHER == False:
    rule all:
        input:
            expand(bd("scaffolds/{scaffold}_GL.beagle.gz"), scaffold = scaffolds)
else: 
    rule all:
        input:
            bd("genome_wide_GL_filtered.beagle.gz"),
            sites=bd("genome_wide_GL_filtered.sites")


# calculate genotype likelihoods separately for each scaffold
# Simultanesouly does filering,
# calculating of genotype
rule calc_beagle_GLs_by_scaff:
    input:
        bamlist=BAM_LIST_FILE,
        ref=REF
    output:
        bd("scaffolds/{scaffold}_GL.beagle.gz")
    params:
        basename=lambda wildcards: bd("scaffolds/" + wildcards.scaffold + "_GL"),
        region=lambda wildcards: wildcards.scaffold + ":",
        SNP_pval = config["SNP_pval"],
        minMapQ = config["minMapQ"],
        minQ = config["minQ"],
        uniqueOnly = config["uniqueOnly"],
        remove_bads = config["remove_bads"],
        only_proper_pairs = config["only_proper_pairs"],
        trim = config["trim"],
        C = config["C"],
        baq = config["baq"],
        minInd = config["minInd"],
        setMinDepth = config["setMinDepth"],
        setMaxDepth = config["setMaxDepth"],
        setMinDepthInd = config["setMinDepthInd"],
        setMaxDepthInd = config["setMaxDepthInd"]
    priority: 10
    log:
        bd("logs/calc_GL_by_scaff/{scaffold}.txt")
    resources:
        cpus=config["cpus"]
    shell:
        """
        angsd -bam {input.bamlist} \
            -ref {input.ref} \
            -anc {input.ref} \
            -out {params.basename} \
            -r {params.region} \
            -gl 1 -doGlf 2 -doCounts 1 -dosaf 1 \
            -doMajorMinor 1 -doMaf 2 \
            -SNP_pval {params.SNP_pval} \
            -minMapQ {params.minMapQ} \
            -minQ {params.minQ} \
            -uniqueOnly {params.uniqueOnly} \
            -remove_bads {params.remove_bads} \
            -only_proper_pairs {params.only_proper_pairs} \
            -trim {params.trim} \
            -C {params.C} \
            -baq {params.baq} \
            -minInd {params.minInd} \
            -setMinDepth {params.setMinDepth} \
            -setMaxDepth {params.setMaxDepth} \
            -setMinDepthInd {params.setMinDepthInd} \
            -setMaxDepthInd {params.setMaxDepthInd} \
            -nThreads  {resources.cpus} > {log} 2>&1
        """

rule gather_per_scaff_GLS:
    input:
        expand(bd("scaffolds/{scaffold}_GL.beagle.gz"), scaffold = scaffolds)
    output:
        glf=bd("genome_wide_GL_filtered.beagle.gz"),
        sites=bd("genome_wide_GL_filtered.sites")
    shell:
        """
        cat {input} > {output.glf}
        zcat {output.glf} | awk -F '\t' '{{print $1}}' > {output.sites}
        """