network_construction <- function(output_directory,counts_table,traits_table,minModuleSize,power_threshold = 0, merg_thresh, signed, merged_table) {

  orig_dir <- getwd()
  setwd(output_directory)

  if (power_threshold != 0) {
    power <- power_threshold
  } else {

  print("Selecting power threshold for the scale free topology model")
  powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
  sft <- pickSoftThreshold(data = counts_table,
                           powerVector = powers,
                           verbose = 3)
  power <- sft$powerEstimate
  pdf(file = "power_selection.pdf",
      width = 9,
      height = 5)
  sizeGrWindow(9, 5)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit,signed R^2",
       type = "n",
       main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 0.9, col = "red")
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
       cex = 0.9, col = "red")
  dev.off()

  }

  print(paste("The selected power threshold for building the scale free topology model is ",power))

  print("Creating adjacency and topology overlap matrix (TOM)...")

  temp_cor <- cor
  cor <- WGCNA::cor

  adjacency <- adjacency(datExpr = counts_table, power = power, type = signed)
  TOM <- TOMsimilarity(adjacency)
  print("Finished TOM creation...")
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  pdf(file = "TOM_dendrogram.pdf", width = 12, height = 9)
  sizeGrWindow(12, 9)
  plot(x = geneTree,
       xlab = "",
       sub = "",
       main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE,
       hang = 0.04)
  dev.off()

  cor <- temp_cor

  dynamicMods <- cutreeDynamic(dendro = geneTree,
                               distM = dissTOM,
                               deepSplit = 2,
                               pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)

  dynamicColors <- labels2colors(dynamicMods)
  pdf(file = "dendrogram_dynamic_tree_cut.pdf",
      width = 8,
      height = 6)
  sizeGrWindow(8, 6)
  plotDendroAndColors(dendro = geneTree,
                      colors = dynamicColors,
                      "Dynamic Tree Cut",
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()

  MEList <- moduleEigengenes(expr = counts_table, colors = dynamicColors)
  MEs <- MEList$eigengenes
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  pdf(file = "cluster_modules_eigengenes.pdf",
      width = 7,
      height = 6)
  sizeGrWindow(7, 6)
  plot(x = METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = merg_thresh, col = "red")
  dev.off()

  print("Started merging modules...")

  merge <- mergeCloseModules(exprData = counts_table,
                             colors = dynamicColors,
                             cutHeight = merg_thresh,
                             verbose = 3)

  mergedColors <- merge$colors
  num_colors <- as.data.table(table(mergedColors))
  fwrite(num_colors, "num_genes_per_module.tsv")

  mergedMEs <- merge$newMEs
  module_genes <- cbind(merged_table,mergedColors)
  fwrite(module_genes, "genes_name_per_module.tsv")

  print("Finished merging modules... Starting plotting...")

  pdf(file = "gene_dendro_merged.pdf",
      width = 12,
      height = 9)
  sizeGrWindow(12, 9)
  plotDendroAndColors(dendro = geneTree,
                      colors = cbind(dynamicColors,mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  moduleColors <- mergedColors
  colorOrder <- c("grey", standardColors(50))
  moduleLabels <- match(moduleColors, colorOrder) - 1
  ME <- mergedMEs
  MEs0 <- moduleEigengenes(counts_table, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, traits_table, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(counts_table))
  pdf("heatmap_traits.pdf",
      width = 10,
      height = 6)
  sizeGrWindow(10, 6)
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(traits_table),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1, 1),
                 main = paste("Module-trait relationship"))
  dev.off()

  print("Finished plotting...")
  #plotTOM <- dissTOM^7
  #diag(plotTOM) <- NA
  #pdf(file = "network_heatmap.pdf", width = 9, height = 9)
  #sizeGrWindow(9, 9)
  #TOMplot(dissim = plotTOM,
  #        dendro = geneTree,
  #        Colors = moduleColors,
  #        main = "Network heatmap plot, all genes")
  #dev.off()

  setwd(orig_dir)

  return(list(MEs, moduleTraitCor, moduleTraitPvalue, moduleColors, TOM, dissTOM))
}
