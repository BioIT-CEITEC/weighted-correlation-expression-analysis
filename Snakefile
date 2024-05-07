import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
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

analysis = []
if config["featureCount"]:
    count_over_list = config['count_over'].split(",")
    if ("exon" in count_over_list):
        config["featureCount_exon"] = True
        analysis.append("featureCount_exon")
    if ("gene" in count_over_list):
        config["featureCount_gene"] = True
        analysis.append("featureCount_gene")
    if ("transcript" in count_over_list):
        config["featureCount_transcript"] = True
        analysis.append("featureCount_transcript")
    if ("three_prime_UTR" in count_over_list):
        config["featureCount_3pUTR"] = True
        analysis.append("featureCount_3pUTR")
    if ("five_prime_UTR" in count_over_list):
        config["featureCount_5pUTR"] = True
        analysis.append("featureCount_5pUTR")

if config["RSEM"]:
    analysis.append("RSEM")
if config["salmon_align"]:
    analysis.append("salmon_align")
if config["salmon_map"]:
    analysis.append("salmon_map")
if config["kallisto"]:
    analysis.append("kallisto")
if len(analysis) == 0:
    raise ValueError("There was no RSEM or featureCount used in previous analysis!")


config["analysis_type"] = "|".join(analysis)


wildcard_constraints:
     sample = "|".join(sample_tab.sample_name), 
     lib_name="[^\.\/]+",
     analysis_type= "featureCount_exon|featureCount_gene|featureCount_transcript|featureCount_3pUTRn|featureCount_5pUTR|RSEM|salmon_map|salmon_align|kallisto"


rule all:
    input: "results/cluster_modules_eigengenes.pdf"

##### Modules #####

include: "rules/wgcna_analysis.smk"