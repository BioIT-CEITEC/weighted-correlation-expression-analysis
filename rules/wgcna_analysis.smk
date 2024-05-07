rule WGCNA_computation:
    input: counts = expand("DE_{analysis_type}/count_dt.RDS", analysis_type = config["analysis_type"])[0],
           gtf = config["organism_gtf"]
    output: "results/cluster_modules_eigengenes.pdf"
    params: sample_tab = sample_tab,
            signed = config["signed_analysis"],
            min_module_size = config["min_mod_size"],
            merge_threshold = config["merg_thresh"],
            power_selection = config["power_sel"],
            experiment_design = "results/WGCNA_experiment_design.tsv",
            trait_names = config["trait_names"]
    log: "logs/WGCNA/WGCNA_summary.log"
    conda: "../wrappers/WGCNA_computation/env.yaml"
    script: "../wrappers/WGCNA_computation/script.py"