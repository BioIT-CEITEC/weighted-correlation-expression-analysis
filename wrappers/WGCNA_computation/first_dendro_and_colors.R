first_dendro_and_colors <- function(normalized_counts, traits_dt,
                                    output_dir, signed) {

  orig_dir <- getwd()
  dir.create(path = output_dir, showWarnings = FALSE, recursive = TRUE)
  setwd(output_dir)

  gsg <- goodSamplesGenes(datExpr = normalized_counts, verbose = 3)

  badGenes <- colnames(normalized_counts)[!gsg$goodGenes]
  badSamples <- rownames(normalized_counts)[!gsg$goodSamples]

  fwrite(as.data.frame(badGenes), "removed_genes.tsv", sep = "\t")
  fwrite(as.data.frame(badSamples), "removed_samples.tsv", sep = "\t")

  print(paste("Filtering... ", length(badGenes), "genes from a total of",
              length(gsg$goodGenes), ", and", length(badSamples),
              "samples from a total of", length(gsg$goodSamples),
              "will be removed, flagged as low quality"))

  filtered_counts <- normalized_counts[gsg$goodSamples, gsg$goodGenes]

  sampleTree <- hclust(dist(filtered_counts), method = "average")
  pdf(file = "sample_clustering.pdf",
      width = 12,
      height = 9)

  sizeGrWindow(12, 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  plot(sampleTree, main = "Sample clustering to detect outliers",
       sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()

  traitColors <- numbers2colors(traits_dt, signed = signed)

  pdf(file = "sample_trait_dendrogram.pdf",
      width = 12,
      height = 9)

  sizeGrWindow(12, 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))

  plotDendroAndColors(dendro = sampleTree,
                      colors = traitColors,
                      groupLabels = names(traits_dt),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()

  setwd(orig_dir)

  return(filtered_counts)
}