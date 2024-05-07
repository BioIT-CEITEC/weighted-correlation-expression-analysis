load_counts_and_traits <- function(counts, traits, new_colnames) {

  counts_dt <- readRDS(counts)
  counts_dt <- dcast(counts_dt, Ensembl_Id + Feature_name ~ sample_name,
                     value.var = "vstcounts",
                     fill = 0)
  counts_dt <- counts_dt[, -2]
  counts_dt <- tibble::column_to_rownames(counts_dt, "Ensembl_Id")


  traits_dt <- fread(traits)
  traits_dt <- tibble::column_to_rownames(traits_dt, "sample_name")
  new_names <- strsplit(new_colnames, split = ",")
  colnames(traits_dt) <- new_names[[1]]

  return(list(counts_dt, traits_dt))
}