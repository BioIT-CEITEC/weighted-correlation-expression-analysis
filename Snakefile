import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalReferences"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok = True)

#### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

config = BR.load_organism()

wildcards_constraints:
    sample = "|".join(sample_tab.sample_name)

rule all:
    input: "results/cluster_modules_eigengenes.pdf"

##### Modules #####

include: "rules/wgcna_analysis.smk"