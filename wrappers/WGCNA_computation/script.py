#############################################################
# wrapper for rule: WGCNA_computation
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: WGCNA_computation \n##\n")
f.close()

snakemake.params.sample_tab.to_csv(snakemake.params.experiment_design, index=False,sep="\t")


command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/wgcna_computation.R " +\
            snakemake.params.experiment_design + " " +\
            snakemake.input.counts + " " +\
            snakemake.input.gtf + " " +\
            str(snakemake.params.signed) + " " +\
            str(snakemake.params.min_module_size) + " " +\
            str(snakemake.params.merge_threshold) + " " +\
            str(snakemake.params.power_selection) + " " +\
            str(snakemake.params.trait_names) + \
            " >> " + log_filename + " 2>&1"

f = open(log_filename, "a+")
f.write("##COMMAND: " + command + "\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
shell(command)