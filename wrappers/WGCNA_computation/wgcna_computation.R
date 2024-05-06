library(WGCNA)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)


run_all <- function(args) {
  experiment_design_file <- args[1]
  counts_file <- args[2]
  traits_table <- args[3]
  relative_output_dir <- dirname(experiment_design_file)
  gtf_filename <- args[4]
  signed <- as.logical(toupper(args[5]))
  min_module_size <- as.numeric(args[6])
  merging_threshold <- as.numeric(args[7])
  power_threshold <- as.numeric(args[8])

  allowWGCNAThreads(nThreads = 30)
  options(stringsAsFactors = FALSE)

  output_dir <- paste0(getwd(), "/results")

  if (signed) {
    type <- "signed"
  } else {
    type <- "unsigned"
  }

  normalized_counts <- t(fread(counts_file) %>% column_to_rownames("rn"))
  traits_dt <- fread(traits_table) %>% column_to_rownames("rn")

  final_counts <- first_dendro_and_colors(counts_table = normalized_counts,
                                          traits_table = traits_dt,
                                          output_directory = output_dir,
                                          signed = signed)

  allIDs <- as.data.table(colnames(final_counts))
  feat_type <- "gene"
  annotate_by <- c("gene_id", "gene_name")
  gtf_gene_tab <- as.data.table(rtracklayer::import(gtf_filename, feature.type = feat_type))[,annotate_by,with=F]
  setnames(gtf_gene_tab, c("Geneid", "gene_name"))
  merged_table <- merge.data.table(allIDs, gtf_gene_tab, by.x = "V1",
                                   by.y = "Geneid", all.x = TRUE)

  setnames(merged_table, "V1", "Geneid")

  res_list <- network_construction(output_directory = output_dir,
                                   counts_table = final_counts,
                                   traits_table = traits_dt,
                                   minModuleSize = min_module_size,
                                   merg_thresh = merging_threshold,
                                   signed = type,
                                   merged_table = merged_table,
                                   power_threshold = power_threshold)
  MEs <- res_list[[1]]
  moduleTraitCor <- res_list[[2]]
  moduleTraitPvalue <- res_list[[3]]
  moduleColors <- res_list[[4]]
  TOM <- res_list[[5]]
  dissTOM <- res_list[[6]]

  trait_relationship(output_dir = output_dir,
                     counts_dt = final_counts,
                     traits_table = traits_dt,
                     eigenvalues = MEs,
                     colors =  moduleColors,
                     merged_table = merged_table,
                     gtf = gtf_gene_tab,
                     tom = TOM)


}

script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/first_dendro_and_colors.R"))
source(paste0(script.dir, "/network_construction.R"))
source(paste0(script.dir, "/trait_relationships.R"))

# run as Rscript
args <- commandArgs(trailingOnly = TRUE)
run_all(args)
