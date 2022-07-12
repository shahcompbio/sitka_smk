library(tidyverse)
library(data.table)
library(cowplot)
library(ComplexHeatmap)
library(ape)
snakemake@source("tree_utils.R")
devtools::load_all(snakemake@config$schnapps)

cn <- fread(snakemake@input$hscn)
tree <- ape::read.tree(snakemake@input$tree)

tip.loci <- grep('locus', tree$tip.label, value = T)
while (length(tip.loci) > 0) {
  tree <- ape::drop.tip(tree, tip.loci, trim.internal = FALSE, collapse.singles = FALSE)
  tip.loci <- grep('locus', tree$tip.label, value = T)
}
tree$tip.label <- str_remove(tree$tip.label, "cell_")
tree <- ape::compute.brlen(tree, 1)

tree <- remove_fingers(tree)

cells <- intersect(unique(cn$cell_id), tree$tip.label)
cn <- filter(cn, cell_id %in% cells)
tree <- ape::keep.tip(tree, cells)

chroms <- unique(cn$chr)
chroms <- chroms[!chroms %in% c("14", "16", "18", "19", "21", "22")]
mycl <- data.frame(cell_id = tree$tip.label, clone_id = "0")

sample_label_idx <- 1

cn <- filter(cn, !is.na(copy))

message("Generating plots")
h1chr <- plotHeatmap(cn, 
                     chrlabels = chroms,
                     reorderclusters = FALSE, 
                     tree = tree, 
                     clusters = mycl,
                     sample_label_idx = sample_label_idx,
                     show_library_label = TRUE,
                     spacer_cols = 15, 
                     show_legend = TRUE, 
                     show_clone_label = TRUE, 
                     plotfrequency = TRUE, 
                     plottree = TRUE)
h2chr <- plotHeatmap(cn, 
                     chrlabels = chroms,
                     plotcol = "state_phase",
                     reorderclusters = TRUE, 
                     tree = tree, 
                     clusters = mycl,
                     sample_label_idx = sample_label_idx,
                     show_library_label = FALSE,
                     spacer_cols = 15, 
                     show_legend = TRUE, 
                     show_clone_label = FALSE, 
                     plotfrequency = TRUE, 
                     plottree = FALSE)


ncells <- length(unique(cn$cell_id))
mytitle <- paste0(snakemake@wildcards$sample, " (", ncells, " cells)")

message("Saving pdf")
pdf(snakemake@output$plot, width = 17)
ComplexHeatmap::draw(h1chr + h2chr, ht_gap = unit(0.6, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom", column_title = mytitle, column_title_gp = grid::gpar(fontsize = 10))
dev.off()

message("Saving png")
png(snakemake@output$plotpng, width = 1152)
ComplexHeatmap::draw(h1chr + h2chr, ht_gap = unit(0.6, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom", column_title = mytitle, column_title_gp = grid::gpar(fontsize = 10))
dev.off()
