trait_relationship <- function(output_dir, counts_dt, traits_table, eigenvalues, colors, merged_table, gtf, tom) {

  orig_dir <- getwd()
  working_dir <- paste0(output_dir, "/trait_relationship")
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  setwd(working_dir)

  cor <- WGCNA::cor

  print("Starting trait relationship analysis...")

  nSamples <- nrow(counts_dt)

  gene_symbol <- as.character(merged_table$gene_name)
  moduleColors <- colors
  geneModuleMembership <- as.data.frame(cor(counts_dt, eigenvalues,
                                            use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),
                                             nSamples))
  modNames <- substring(names(eigenvalues), 3)
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  mod_memb <- cbind(geneModuleMembership, MMPvalue)
  mod_memb <- tibble::rownames_to_column(mod_memb, var = "Geneid")
  mod_memb_symb <- merge.data.table(mod_memb, merged_table,
                                    by.x = "Geneid",
                                    by.y = "Geneid",
                                    all.x = TRUE)
  setcolorder(mod_memb_symb, c("Geneid", "gene_name"))
  fwrite(as.data.frame(mod_memb_symb), file = "gene-module-membership.tsv")

  for (trait in colnames(traits_table)) {
    temp_dir <- paste0(working_dir, "/", trait)
    dir.create(path = temp_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(temp_dir)
    trait_data <- as.data.frame(traits_table[[trait]])
    GS1 <- as.numeric(cor(counts_dt, trait_data, use = "p"))
    gene_significance <- abs(GS1)

    pdf(file = paste0(trait, "-gene_significance.pdf"), width = 16, height = 14)
    sizeGrWindow(16, 14)
    par(mfrow = c(1, 1))
    plotModuleSignificance(gene_significance, moduleColors,
                           ylab = paste("Gene Significance for", trait))
    dev.off()

    geneTraitSignificance <- as.data.frame(cor(counts_dt, trait_data,
                                               use = "p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),
                                               nrow(counts_dt)))
    names(geneTraitSignificance) <- paste0("GS.", trait)
    names(GSPvalue) <- paste0("p.GS.", trait, sep = "")
    trait_sig <- cbind(geneTraitSignificance, GSPvalue)
    trait_sig <- tibble::rownames_to_column(trait_sig, var = "Geneid")
    gen_symb_merg <- merge.data.table(trait_sig, gtf,
                                      by.x = "Geneid",
                                      by.y = "Geneid",
                                      all.x = TRUE)
    setcolorder(gen_symb_merg, c("Geneid", "gene_name"))
    filename <- paste0(trait, "-significance.tsv")
    fwrite(as.data.frame(gen_symb_merg), file = filename)

    modOrder <- order(-abs(cor(eigenvalues, trait_data, use = "p")))
    MET <- orderMEs(cbind(eigenvalues, trait_data))
    pdf(file = paste0(trait, "-network.pdf"), width = 10, 15)
    sizeGrWindow(10, 15)
    par(cex = 0.9)
    plotEigengeneNetworks(multiME = MET,
                          setLabels = "",
                          marDendro = c(0, 4, 1, 2),
                          marHeatmap = c(3, 4, 1, 2),
                          cex.lab = 0.8,
                          xLabelsAngle = 90)
    dev.off()

    pdf(file = paste0(trait, "-eigengene_dendrogram.pdf"),
        width = 12,
        height = 12)
    sizeGrWindow(12, 12)
    par(cex = 1.0)
    plotEigengeneNetworks(multiME = MET,
                          "Eigengene dendrogram",
                          marDendro = c(0, 4, 2, 0),
                          plotHeatmaps = FALSE)
    dev.off()

    pdf(file = paste0(trait, "-eigengene_adjacency_heatmap.pdf"),
        width = 12,
        height = 12)
    sizeGrWindow(12, 12)
    par(cex = 1.0)
    plotEigengeneNetworks(multiME = MET,
                          "Eigengene adjacency heatmap",
                          marHeatmap = c(3, 4, 2, 2),
                          plotDendrograms = FALSE,
                          xLabelsAngle = 90)
    dev.off()

    for (mod in unique(moduleColors)) {
      print(paste("The module membership and cytoscape networks for", mod,
                  "module in the condition:", trait, "is being processed",
                  sep = " "))
      column <- match(mod, modNames)
      moduleGenes <- moduleColors == mod
      pdf(file = paste0(trait, "-", mod, "-scatterplot.pdf"),
          width = 14,
          height = 14)
      sizeGrWindow(14, 14)
      par(mfrow = c(1, 1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module membership in", mod, "module"),
                         ylab = paste("Gene significance for", trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
      dev.off()

      modGenes <- (moduleColors == mod)
      fileName <- paste(trait, "-IDs-", mod, ".tsv", sep = "")
      fwrite(as.data.frame(merged_table[modGenes]), file = fileName)

      gene_name <- colnames(counts_dt)
      inModule <- is.finite(match(moduleColors, mod))
      modProbes <- gene_name[inModule]
      modGenes <- gene_symbol[inModule]
      modTOM <- tom[inModule, inModule]
      dimnames(modTOM) <- list(modProbes, modProbes)
      cyt <- exportNetworkToCytoscape(adjMat = modTOM,
                                      edgeFile = paste0(trait,"-cytoscapeInput-edges-", paste(mod,collapse="-"), ".tsv"),
                                      nodeFile = paste0(trait,"-cytoscapeInput-nodes-", paste(mod,collapse="-"),".tsv"),
                                      weighted = TRUE,
                                      threshold = 0.02,
                                      nodeNames = modProbes,
                                      altNodeNames = modGenes,
                                      nodeAttr = moduleColors[inModule])

    }

    setwd(working_dir)
  }
  setwd(orig_dir)
}
