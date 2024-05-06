rule WGCNA_computation:
    input: experiment_design = "experiment_design.tsv",
           gtf = config["organism_gtf"],
           counts = "normalized_counts.tsv",
           traits = "traits.tsv"
    output: "rules/cluster_modules_eigengenes.pdf"
    params: signed = config["signed_analysis"],
            min_module_size = config["min_mod_size"],
            merge_threshold = config["merg_thresh"],
            power_selection = config["power_sel"]
    log: "logs/WGCNA/WGCNA_summary.log"
    conda: "../wrappers/WGCNA_computation/env.yaml"
    script: "../wrappers/WGCNA_computation/script.py"