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
hscn <- subset(hscn, select = c(chr, start, end, cell_id, Min, Maj))
hscnnew <- hscn[chr != "Y"]

message("Create segments")
options("scipen"=20)

cn_jitterA <- jitter_fix_breakpoints(hscnnew %>% dplyr::rename(state = Maj), #hack rename Maj to state
                                     nextend = snakemake@config$window) %>%
  as.data.table() %>% 
  .[, loci := paste(paste0(chr, "A"), end - 0.5e6 + 1, end, sep = "_")] %>% 
  .[, row_num := .I] %>% # add row numbers
  .[, remove_row_num := .I[.N], by=.(cell_id, chr)]# find last row in each cell_id - chr group

cn_jitterB <- jitter_fix_breakpoints(hscnnew %>% dplyr::rename(state = Min), 
                                     nextend = snakemake@config$window) %>%
  as.data.table() %>% 
  .[, loci := paste(paste0(chr, "B"), end - 0.5e6 + 1, end, sep = "_")] %>% 
  .[, row_num := .I] %>% # add row numbers
  .[, remove_row_num := .I[.N], by=.(cell_id, chr)]# find last row in each cell_id - chr group

cn_jitter <- rbindlist(list(cn_jitterA, cn_jitterB))

segsA <- signals::create_segments(cn_jitterA) %>% mutate(allele := "A")
segsB <- signals::create_segments(cn_jitterB) %>% mutate(allele := "B")
segs <- rbindlist(list(segsA, segsB))
fwrite(segs, file = snakemake@output[["sitka_segs"]])

cn_transitionsA <- signals::create_cntransitions(cn_jitterA) %>% as.data.table() %>% mutate(allele := "A")
cn_transitionsB <- signals::create_cntransitions(cn_jitterB) %>% as.data.table() %>% mutate(allele := "B")
cn_transitions <- rbindlist(list(cn_transitionsA, cn_transitionsB))
fwrite(cn_transitions, file = snakemake@output[["sitka_transitions"]])

cn_transitions <- cn_transitions %>%
  .[, tipInclusionProbabilities := 1] %>%
  .[, loci := paste(paste0(chr, allele), start, end, sep = "_")] %>% 
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

message("Write file")
fwrite(sitka_df, file = snakemake@output[["sitka_input"]])
