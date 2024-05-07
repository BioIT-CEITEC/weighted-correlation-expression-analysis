load_counts_and_traits <- function(counts_file, experiment_design_file, trait_names) {

  counts_dt <- readRDS(counts_file) %>%
    dcast(., Ensembl_Id + Feature_name ~ sample_name, value.var = "vstcounts", fill = 0) %>%
    .[, -2] %>%
    tibble::column_to_rownames(counts_dt, "Ensembl_Id")


  traits_dt <- fread(experiment_design_file)
  traits_dt <- tibble::column_to_rownames(traits_dt, "sample_name")
  new_names <- strsplit(trait_names, split = ",")
  colnames(traits_dt) <- new_names[[1]]

  return(list(counts_dt, traits_dt))
}