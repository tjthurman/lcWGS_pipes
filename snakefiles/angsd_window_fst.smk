# A generalized snakemake pipeline
# for doing windowed Fst analysis

# Assumes .saf files for both populations have
# already been created with the angsd_GL_genome_wide.smk pipeline



# for a given bamlist. 

import os

#################
##   GLOBALS   ##
################# 

# List the reference genome you want to map to:
REF = config["ref_genome"]

# The bamlist containing the list of bam files for input into angsd
# Pop 1
BAM_LIST_FILE_POP1 = config["pop1_bamlist"]
BAM_LIST_NAME_POP1 = os.path.splitext(os.path.basename(BAM_LIST_FILE_POP1))[0]
# Pop 2
BAM_LIST_FILE_POP2 = config["pop2_bamlist"]
BAM_LIST_NAME_POP2 = os.path.splitext(os.path.basename(BAM_LIST_FILE_POP2))[0]

# Batch names
BATCH_NAME_POP1 = config["pop1_batch_name"]
BATCH_NAME_POP2 = config["pop2_batch_name"]
BATCH_NAME_FST = config["fst_batch_name"]

# Scaffold used to get prior for 2d SFS
SCAFF_FOR_2d_SFS = config["scaffold_prior"]

# window metrics
SIZE = config["win_size"]
STEP = config["win_step"]

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

def bd_pop1_in(filepath):
    return os.path.normpath(os.path.join("results", "angsd_GLs", BAM_LIST_NAME_POP1, BATCH_NAME_POP1, filepath))

def bd_pop2_in(filepath):
    return os.path.normpath(os.path.join("results", "angsd_GLs", BAM_LIST_NAME_POP2, BATCH_NAME_POP2, filepath))

def bd_out(filepath):
    return os.path.normpath(os.path.join("results", "angsd_Fst", BAM_LIST_NAME_POP1  + "_" + BAM_LIST_NAME_POP2, BATCH_NAME_FST, filepath))



####################
## PIPELINE START ##
####################
localrules: all


rule all:
    input:
        bd_out("sfs_2d/" + SCAFF_FOR_2d_SFS + ".ml"), # prior for 2d SFS
        expand(bd_out("by_scaff_idx/{scaffold}.fst.idx"), scaffold = scaffolds), # fst index for each scaffold
        expand(bd_out("fst/winSize{size}/winStep{step}/fst_{scaffold}.txt"), size = SIZE, step = STEP, scaffold = scaffolds)

# First, calculate a 2dSFS. 
# To keep things fast, just using the scaffold with the most SNPs

rule sfs_2d:
    input:
        pop1_saf= bd_pop1_in("scaffolds/" + SCAFF_FOR_2d_SFS + "_GL.saf.idx"),
        pop2_saf= bd_pop2_in("scaffolds/" + SCAFF_FOR_2d_SFS + "_GL.saf.idx")
    output:
        bd_out("sfs_2d/" + SCAFF_FOR_2d_SFS + ".ml")
    resources:
        cpus=config["cpus"]
    log:
        bd_out("logs/sfs_2d/" + SCAFF_FOR_2d_SFS + ".txt")
    shell:
        """
        realSFS -P {resources.cpus} {input.pop1_saf} {input.pop2_saf}  > {output} 2> {log}
        """

rule fst_index_by_scaf:
    input:
        pop1_saf = bd_pop1_in("scaffolds/{scaffold}_GL.saf.idx"),
        pop2_saf = bd_pop2_in("scaffolds/{scaffold}_GL.saf.idx"),
        prior= bd_out("sfs_2d/" + SCAFF_FOR_2d_SFS + ".ml")
    output:
        bd_out("by_scaff_idx/{scaffold}.fst.idx")
    resources:
        cpus=config["cpus"]
    log:
        bd_out("logs/fst_index_by_scaf/{scaffold}.txt")
    params:
        basename=lambda wildcards: bd_out("by_scaff_idx/" + wildcards.scaffold)
    shell:
        """
        realSFS fst index {input.pop1_saf} {input.pop2_saf} -sfs {input.prior} -P {resources.cpus} -fstout {params.basename} -whichFst 1 > {log} 2>&1
        """

rule fst_sliding_window_by_scaff:
    input:
        bd_out("by_scaff_idx/{scaffold}.fst.idx")
    output:
        bd_out("fst/winSize{size}/winStep{step}/fst_{scaffold}.txt")
    log:
        bd_out("logs/fst_sliding_win_by_scaf/winSize{size}/winStep{step}/{scaffold}.txt")
    resources:
        cpus=config["cpus"]
    params:
        stepsize = lambda wildcards: wildcards.step,
        winsize = lambda wildcards: wildcards.size
    shell:
        """
        realSFS fst stats2 {input} -win {params.winsize} -step {params.stepsize} -P {resources.cpus} > {output}  2> {log}
        """
