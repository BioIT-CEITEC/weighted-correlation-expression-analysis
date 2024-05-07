library(WGCNA)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)


run_all <- function(args) {
  experiment_design_file <- args[1]
  counts_file <- args[2]
  gtf_filename <- args[3]
  signed <- as.logical(toupper(args[4]))
  min_module_size <- as.numeric(args[5])
  merging_threshold <- as.numeric(args[6])
  power_threshold <- as.numeric(args[7])
  trait_names <- args[8]


  allowWGCNAThreads(nThreads = 30)
  options(stringsAsFactors = FALSE)

  output_dir <- paste0(getwd(), "/results")

  if (signed) {
    signed_type <- "signed"
  } else {
    signed_type <- "unsigned"
  }

  res_list <- load_counts_and_traits(counts_file,
                                     experiment_design_file,
                                     trait_names)

  normalized_counts <- t(res_list[[1]])
  traits_dt <- res_list[[2]]

  filtered_counts <- first_dendro_and_colors(normalized_counts,
                                             traits_dt,
                                             output_dir,
                                             signed)

  allIDs <- as.data.table(colnames(filtered_counts))
  feat_type <- "gene"
  annotate_by <- c("gene_id", "gene_name")
  gtf_gene_tab <- as.data.table(rtracklayer::import(gtf_filename, feature.type = feat_type))[,annotate_by,with=F]
  setnames(gtf_gene_tab, c("Geneid", "gene_name"))
  merged_table <- merge.data.table(allIDs, gtf_gene_tab, by.x = "V1",
                                   by.y = "Geneid", all.x = TRUE)

  setnames(merged_table, "V1", "Geneid")

  res_list <- network_construction(output_dir,
                                   filtered_counts,
                                   traits_dt,
                                   min_module_size,
                                   merging_threshold,
                                   signed_type,
                                   merged_table,
                                   power_threshold)
  MEs <- res_list[[1]]
  moduleTraitCor <- res_list[[2]]
  moduleTraitPvalue <- res_list[[3]]
  moduleColors <- res_list[[4]]
  TOM <- res_list[[5]]
  dissTOM <- res_list[[6]]

  trait_relationship(output_dir,
                     filtered_counts,
                     traits_dt,
                     MEs,
                     moduleColors,
                     merged_table,
                     TOM)


}

script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/load_counts_traits.R"))
source(paste0(script.dir, "/first_dendro_and_colors.R"))
source(paste0(script.dir, "/network_construction.R"))
source(paste0(script.dir, "/trait_relationships.R"))

# run as Rscript
args <- commandArgs(trailingOnly = TRUE)
run_all(args)