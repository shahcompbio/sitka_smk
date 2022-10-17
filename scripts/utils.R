fixjitter <- function(bps, nextend = 2){
  x <- as.data.frame(colSums(bps))
  names(x) <- "frequency"
  x$loci <- row.names(x)
  
  frequency_order <- order(x$frequency, decreasing = TRUE)
  nloci <- dim(bps)[2]
  
  while(length(frequency_order) > 0){
    idx <- frequency_order[1]
    minidx <- max(1, idx - nextend)
    maxidx <- min(nloci, idx + nextend)
    temp_mat <-  bps[, minidx:maxidx] #extract matrix nextend either side of locus of interest
    bps[, minidx:maxidx] <- 0 #set matrix to 0
    bps[, idx] <- as.numeric(rowSums(temp_mat) > 0)
    frequency_order <- setdiff(frequency_order, minidx:maxidx)
  }
  
  message(paste0("Original number of loci: ", nloci))
  x <- colSums(bps)
  x <- x[x>0]
  bps <- bps[,names(x)]
  message(paste0("New number of loci: ", dim(bps)[2]))
  
  return(as.data.frame(bps))
}

createbreakpointmatrix_startend <- function(segs,
                                            transpose = FALSE,
                                            internalonly = TRUE,
                                            use_state = FALSE,
                                            state_remove = 2){
  
  options("scipen"=20)
  
  if (use_state == FALSE){
    if (internalonly == TRUE){
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci_start := paste(chr, start, start + 0.5e6 - 1, sep = "_")] %>%
        .[, loci_end := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, row_num := .I] %>% # add row numbers
        .[, remove_row_num := .I[.N], by=.(cell_id, chr)]  %>% # find last row in each cell_id - chr group
        .[row_num != remove_row_num] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci_start, loci_end, chr, start, end, tipInclusionProbabilities)
    } else {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[state != state_remove] %>%
        .[, loci_start := paste(chr, start, start + 0.5e6 - 1, sep = "_")] %>%
        .[, loci_end := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci_start, loci_end, chr, start, end, tipInclusionProbabilities)
    }
  } else{
    if (internalonly == TRUE){
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci_start := paste(chr, start, start + 0.5e6 - 1, sep = "_")] %>%
        .[, loci_end := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, row_num := .I] %>% # add row numbers
        .[, remove_row_num := .I[.N], by=.(cell_id, chr)]  %>% # find last row in each cell_id - chr group
        .[row_num != remove_row_num] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci_start, loci_end, chr, start, end, tipInclusionProbabilities)
    } else {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[state != state_remove] %>%
        .[, loci_start := paste(chr, start, start + 0.5e6 - 1, sep = "_")] %>%
        .[, loci_end := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci_start, loci_end, chr, start, end, tipInclusionProbabilities)
    }
  }
  
  mapping <- dplyr::select(segs_bin, cell_id, chr, start, end, loci_start, loci_end)
  
  segs_mat_start <- segs_bin %>%
    dplyr::select(cell_id, loci_start, tipInclusionProbabilities) %>%
    data.table::dcast(., loci_start ~ cell_id, value.var = "tipInclusionProbabilities", fill = 0) %>%
    .[gtools::mixedorder(loci_start)]
  
  segs_mat_start <- as.data.frame(segs_mat_start)
  rownames(segs_mat_start) <- segs_mat_start$loci_start
  
  if (transpose == TRUE){
    segs_mat_start <- subset(segs_mat_start, select = -c(loci_start))
    segs_mat_start <- t(segs_mat_start)
  }
  
  segs_mat_end <- segs_bin %>%
    dplyr::select(cell_id, loci_end, tipInclusionProbabilities) %>%
    data.table::dcast(., loci_end ~ cell_id, value.var = "tipInclusionProbabilities", fill = 0) %>%
    .[gtools::mixedorder(loci_end)]
  
  segs_mat_end <- as.data.frame(segs_mat_end)
  rownames(segs_mat_end) <- segs_mat_end$loci_end
  
  if (transpose == TRUE){
    segs_mat_end <- subset(segs_mat_end, select = -c(loci_end))
    segs_mat_end <- t(segs_mat_end)
  }
  
  return(list(bps_start = segs_mat_start, bps_end = segs_mat_end, mapping = mapping))
}

jitter_fix_breakpoints <- function(CNbins,
                                   is_breakpoints = FALSE,
                                   nextend = 1,
                                   field = "state",
                                   bothends = FALSE){
  
  if (is_breakpoints == FALSE){
    segs <- signals::create_segments(CNbins, field = field)
  } else{
    segs <- CNbins
  }
  bps_ <- createbreakpointmatrix_startend(segs,
                                          transpose = TRUE,
                                          internalonly = FALSE,
                                          state_remove = 30)
  
  bps_start <- fixjitter(bps_$bps_start, nextend = nextend)
  bps_end <- fixjitter(bps_$bps_end, nextend = nextend)
  #bps <- as.data.frame(bps_$bps)
  
  bps_start$cell_id <- row.names(bps_start)
  bps_end$cell_id <- row.names(bps_end)
  
  newsegs_start <- bps_start %>%
    tidyr::pivot_longer(-cell_id, names_to = "loci", values_to = "is_segment") %>%
    dplyr::filter(is_segment == 1) %>%
    tidyr::separate(loci, c("chr", "start", "end"), convert = TRUE, remove = FALSE) %>%
    dplyr::mutate(chr = paste0(chr)) %>%
    dplyr::select(-end) %>%
    dplyr::rename(start_jitter = start, newloci_start = loci)
  
  newsegs_end <- bps_end %>%
    tidyr::pivot_longer(-cell_id, names_to = "loci", values_to = "is_segment") %>%
    dplyr::filter(is_segment == 1) %>%
    tidyr::separate(loci, c("chr", "start", "end"), convert = TRUE, remove = FALSE) %>%
    dplyr::mutate(chr = paste0(chr)) %>%
    dplyr::select(-start) %>%
    dplyr::rename(end_jitter = end, newloci_end = loci)
  
  if (bothends == FALSE){
    segs2 <- dplyr::left_join(bps_$mapping, newsegs_end) %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(dist_end = abs(end - end_jitter) / 1e6) %>%
      dplyr::arrange(cell_id, chr, start, end, dist_end) %>%
      dplyr::group_by(cell_id, chr, start, end) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::group_by(cell_id, chr) %>%
      dplyr::mutate(start_jitter = dplyr::lag(end_jitter) + 1) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(segs) %>%
      dplyr::rename(start_old = start, end_old = end, loci_old = loci_end, loci_jitter = newloci_end) %>%
      dplyr::rename(start = start_jitter, end = end_jitter) %>%
      dplyr::select(cell_id, chr, start, end, state, start_old, end_old, loci_old, loci_jitter) %>%
      dplyr::mutate(start = ifelse(is.na(start), start_old, start)) %>%
      dplyr::filter((end - start) > 0.5e6) %>%
      signals::create_segments()
  } else {
    segs2 <- dplyr::left_join(bps_$mapping, newsegs_end) %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::left_join(newsegs_start) %>%
      dplyr::mutate(dist_start = abs(start - start_jitter) / 1e6) %>% #calculate distance from jittered edge to old edge
      dplyr::mutate(dist_end = abs(end - end_jitter) / 1e6) %>%
      #dplyr::arrange(cell_id, chr, start, end, dist_end) %>%
      #dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::filter(end_jitter >= start_jitter) %>% #ensure end > start
      dplyr::group_by(cell_id, chr, start, end) %>%
      dplyr::mutate(start_jitter = start_jitter[which.min(dist_start)],
                    end_jitter = end_jitter[which.min(dist_end)]) %>% #find jittered edge closest to old edge
      dplyr::ungroup() %>%
      dplyr::filter(end_jitter >= start_jitter) %>%
      dplyr::distinct(cell_id, chr, start, end, start_jitter, end_jitter) %>%
      dplyr::left_join(segs) %>% #the following few lines make sure segments are overlapping
      dplyr::mutate(width_jitter = end_jitter - start_jitter) %>%
      dplyr::mutate(width = end - start) %>%
      dplyr::arrange(desc(width)) %>%
      dplyr::group_by(cell_id, chr, start_jitter, end_jitter) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      orderdf() %>%
      dplyr::mutate(start_lag = dplyr::lag(end_jitter)) %>%
      dplyr::mutate(start_jitter = ifelse(start_jitter > start_lag, start_lag + 1, start_jitter)) %>%
      dplyr::rename(start_old = start, end_old = end) %>%
      dplyr::rename(start = start_jitter, end = end_jitter) %>%
      dplyr::select(cell_id, chr, start, end, state, start_old, end_old) %>%
      dplyr::mutate(start = ifelse(is.na(start), start_old, start)) %>%
      dplyr::filter((end - start) > 0.5e6) %>% #remove segments that are only 1 bin long
      signals::create_segments()
  }
  
  return(segs2)
}

create_edge_table <- function(my_tree){
  node_labels_in_edge <- my_tree$node.label[my_tree$edge[,1]-ape::Ntip(my_tree)]
  tips_nodes <- my_tree$edge[,2]
  
  select.tip.or.node <- function(element, tree) {
    ifelse(element < ape::Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-ape::Ntip(tree)])
  }
  
  edge_table <- data.frame(
    "parent" = my_tree$edge[,1],
    "par.name" = sapply(my_tree$edge[,1], select.tip.or.node, tree = my_tree),
    "child" = my_tree$edge[,2],
    "chi.name" = sapply(my_tree$edge[,2], select.tip.or.node, tree = my_tree)
  )
}

mygetDescendants<-function(tree,node,curr=NULL){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  if(length(curr)==0&&node<=Ntip(tree)) curr<-node
  w<-which(daughters>Ntip(tree))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-mygetDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

add_desc_singletons <- function(et, tree){
  
  et <- left_join(et, data.frame(parent = 1:max(et$parent), nchild = tabulate(et$parent)), by = "parent")
  et$is_cell <- !str_detect(et$chi.name, "locus")
  
  #count number of descendents (only cells) and the number of singleton (only 1 child) descendents (loci + cells)
  #fingers cells tend to be composed of chains of singletons with a few cells at the end/interspersed
  #we therefore look for nodes that have a small number of cell descendents but large numbers of singleton descendents
  
  n_desc <- c()
  n_singletons <- c()
  for (i in 1:nrow(et)){
    desc <- et %>% filter(child %in% mygetDescendants(tree, et$parent[i]))
    
    n_desc <-  c(n_desc, desc %>% 
                   filter(!str_detect(chi.name, "locus")) %>% nrow())
    
    n_children <- filter(et, parent %in% desc$child) %>% pull(nchild)
    
    n_singletons <- c(n_singletons, sum(n_children == 1))
  }
  
  et$descendents <- n_desc
  et$singleton_descendents <- n_singletons
  
  et$descendents_to_singletons <- n_singletons / n_desc  
  
  return(et)
}


get_cells_to_remove <- function(et, tree, ratio = NULL, frac_cells = 0.03){
  
  if (max(et$descendents_to_singletons) < 0.2){
    return(c())
  }
  
  if (frac_cells == 0.0){
    return(c())
  }
           
  total_cells <- sum(et$is_cell)
  n_remove_cells <- round(frac_cells * total_cells)
  if (n_remove_cells < 1){
    return(c())
  }        
  
  if (is.null(ratio)){
    vals <- rev(sort(et$descendents_to_singletons))
    i <- 1
    ni <- 0
    while (ni < n_remove_cells){
      nodes_to_remove <- filter(et, descendents_to_singletons > vals[i]) %>% pull(parent)
      descendents_to_remove <- unique(unlist(lapply(nodes_to_remove, function(x) mygetDescendants(tree, x))))
      cells_to_remove <- et %>% filter(child %in% descendents_to_remove) %>% pull(chi.name)
      cells_to_remove <- cells_to_remove[!str_detect(cells_to_remove, "locus")] 
      i <- i + 1
      ni <- length(cells_to_remove)
    }
  } else {
    nodes_to_remove <- filter(et, descendents_to_singletons > ratio) %>% pull(parent)
    descendents_to_remove <- unique(unlist(lapply(nodes_to_remove, function(x) mygetDescendants(tree, x))))
    cells_to_remove <- et %>% filter(child %in% descendents_to_remove) %>% pull(chi.name)
    cells_to_remove <- cells_to_remove[!str_detect(cells_to_remove, "locus")] 
  }
  
  return(cells_to_remove)
}

remove_fingers <- function(tree, ratio = NULL, frac_cells = 0.03){
  et <- tree %>% create_edge_table(.)
  et <- add_desc_singletons(et, tree)
  
  cells_to_remove <- get_cells_to_remove(et, tree, ratio = ratio, frac_cells = frac_cells)
  
  if (round(length(cells_to_remove) / Ntip(tree), 3) > 0.5){
    cells_to_remove <- get_cells_to_remove(et, tree, ratio = ratio, frac_cells = 0.0)
  }
  
  message(paste0("Removing ", length(cells_to_remove), " (", round(length(cells_to_remove) / Ntip(tree), 3) * 100, "%) cells"))
  
  tree <- ape::drop.tip(tree, cells_to_remove)
  return(tree) 
}
