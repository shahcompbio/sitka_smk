library(tidyverse)
library(data.table)
library(cowplot)
library(ComplexHeatmap)
library(ape)
library(signals)

snakemake@source("utils.R")

tree <- ape::read.tree(snakemake@input$tree)

tip.loci <- grep('locus', tree$tip.label, value = T)
while (length(tip.loci) > 0) {
  tree <- ape::drop.tip(tree, tip.loci, trim.internal = FALSE, collapse.singles = FALSE)
  tip.loci <- grep('locus', tree$tip.label, value = T)
}
tree$tip.label <- str_remove(tree$tip.label, "cell_")

tree <- remove_fingers(tree,  frac_cells = snakemake@config$frac_cells)
ape::write.tree(tree, file = snakemake@output$tree)