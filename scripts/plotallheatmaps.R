library(tidyverse)
library(data.table)
library(cowplot)
library(signals)
library(ape)

hscn <- fread(snakemake@input$hscn)[chr != "Y"]

format_tree <- function(tree){
  tip.loci <- grep('locus', tree$tip.label, value = T)
  while (length(tip.loci) > 0) {
    tree <- ape::drop.tip(tree, tip.loci, trim.internal = FALSE, collapse.singles = FALSE)
    tip.loci <- grep('locus', tree$tip.label, value = T)
  }
  tree$tip.label <- str_remove(tree$tip.label, "cell_")
  tree <- ape::compute.brlen(tree, 1)
  return(tree)
}

hmap <- function(hscn, tree, mytitle = NULL, textsize = 15){
  
  cells <- intersect(hscn$cell_id, tree$tip.label)
  hscn_for_plot <- hscn[cell_id %in% cells]
  tree_for_plot <- keep.tip(tree, cells)
  
  mytitle <- paste0(mytitle, paste0("  ", length(cells), " cells"))
  
  h1 <- plotHeatmap(hscn_for_plot, raster_device = "png",
                    tree = tree_for_plot %>% compute.brlen(.,1), 
                    annofontsize = textsize,
                    clusters = data.frame(cell_id = cells, clone_id = "0"),
                    reorderclusters = T, 
                    linkheight = 2,
                    show_clone_label = F)
  ComplexHeatmap::draw(h1,
                       ht_gap = unit(0.6, "cm"),
                       column_title = mytitle,
                       column_title_gp = grid::gpar(fontsize = textsize + 10),
                       heatmap_legend_side = "bottom",
                       annotation_legend_side = "bottom",
                       show_heatmap_legend = TRUE)
}


pdf(snakemake@output$pdf, width = 15, height = 7)
for (tr in snakemake@input$trees_filt){
  message(tr)
  treetype <- strsplit(basename(tr), "-")[[1]]
  treetype <- treetype[length(treetype) - 2]
  tree <- ape::read.tree(tr)
  hmap(hscn, tree, mytitle = paste0("Filtered ", treetype))
}

for (tr in snakemake@input$trees_unfilt){
  message(tr)
  treetype <- strsplit(basename(tr), "-")[[1]]
  treetype <- treetype[length(treetype) - 1]
  tree <- ape::read.tree(tr) %>% format_tree()
  hmap(hscn, tree, mytitle = paste0("Unfiltered ", treetype))
}

dev.off()
