library(data.table)
library(tidyverse)
library(devtools)
library(signals)

snakemake@source("utils.R")

message("Read in data")
hscn <- data.table::fread(snakemake@input$hscn)
ncells_original <- length(unique(hscn$cell_id))
message(paste0("Number of cells (original): ", ncells_original))
message("Filter data")
hscn <- subset(hscn, select = c(chr, start, end, cell_id, state))
hscnnew <- hscn[chr != "Y"]

message("Create segments")
options("scipen"=20)

cn_jitter <- jitter_fix_breakpoints(hscnnew, nextend = snakemake@config$window) %>%
  as.data.table() %>% 
  .[, loci := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>% 
  .[, row_num := .I] %>% # add row numbers
  .[, remove_row_num := .I[.N], by=.(cell_id, chr)]# find last row in each cell_id - chr group

segs <- signals::create_segments(cn_jitter)
fwrite(segs, file = snakemake@output[["sitka_segs"]])

cn_transitions <- signals::create_cntransitions(cn_jitter) %>% as.data.table()
fwrite(cn_transitions, file = snakemake@output[["sitka_transitions"]])

cn_transitions <- cn_transitions %>%
  .[, tipInclusionProbabilities := 1] %>%
  .[, loci := paste(chr, start, end, sep = "_")] %>% 
  select(cell_id, loci, tipInclusionProbabilities) %>%
  rename(cells = cell_id)
cn_transitions2 <- cn_transitions
message(paste0("Number of transitions: ", dim(cn_transitions)[1]))
message(paste0("Number of unique transitions: ", length(unique(cn_transitions$loci))))
message(paste0("Number of cells: ", length(unique(cn_transitions$cells))))

singleton_segs <- table(cn_transitions$loci)
singleton_segs <- singleton_segs[singleton_segs == 1]
cn_transitions <- cn_transitions[!loci %in% names(singleton_segs)]

message(paste0("Number of transitions: ", dim(cn_transitions)[1]))
message(paste0("Number of unique transitions (removing singletons): ", length(unique(cn_transitions$loci))))
message(paste0("Number of cells: ", length(unique(cn_transitions$cells))))

freq_segs <- table(cn_transitions$loci)
freq_segs <- freq_segs / ncells_original
freq_segs <- freq_segs[freq_segs > snakemake@config$frac_loci]
cn_transitions <- cn_transitions[loci %in% names(freq_segs)]

message(paste0("Number of transitions: ", dim(cn_transitions)[1]))
message(paste0("Number of unique transitions (removing f < 0.001): ", length(unique(cn_transitions$loci))))
message(paste0("Number of cells: ", length(unique(cn_transitions$cells))))


message("Find all segments")
dt <- CJ(unique(hscnnew$cell_id), unique(cn_transitions$loci))
names(dt) <- c("cells", "loci")
sitka_df <- merge(dt, cn_transitions, all = T) %>%
  .[, tipInclusionProbabilities := fifelse(is.na(tipInclusionProbabilities), 0, 1)] %>%
  .[, loci := paste0("locus_", loci)]
message(paste0("Number of cells: ", length(unique(sitka_df$cells))))
message(paste0("Number of cells (original): ", ncells_original))  
# 
# #remove breakpoints associated with replicating cells
# 
# #rank cells by number of breakpoints per cell
# bp_rank <- sitka_df[, list(nbps = sum(tipInclusionProbabilities)), by = "cells"] %>% 
#   .[, r := frank(nbps, ties.method = "random")]
# loci_n <- sitka_df[, list(ncells = sum(tipInclusionProbabilities)), by = "loci"] 
# sitka_df_rank <- sitka_df[bp_rank, on = "cells"]
# sitka_df_rank <- sitka_df_rank[loci_n, on = "loci"]
# 
# sitka_df_rank[tipInclusionProbabilities == 1] %>% 
#   group_by(loci, ncells) %>% 
#   summarize(meanr = mean(r), minr = min(r), maxr = max(r),
#             meanbp = mean(nbps),medianbp = median(nbps), minbp = min(nbps), maxbp = max(nbps)) %>% 
#   ungroup() %>% 
#   mutate(intree = loci %in% tr$node.label) %>% 
#   ggplot(aes(x = meanbp, y = minbp, fill = intree, col = intree)) +
#   geom_point()
# 
# #identify breakpoints enriched in cells with many breakpoints
# sitka_df_rank[tipInclusionProbabilities == 1][, list(x = median(nbps)), by = "loci"] %>%
#   mutate(x2 = loci %in% tr$node.label) %>% 
#   ggplot(aes(x = x, fill = x2)) + geom_histogram()

message("Write file")
fwrite(sitka_df, file = snakemake@output[["sitka_input"]])
