read_count_output <- function(dir, name = dir) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}
read_count_output <- Vectorize(read_count_output, SIMPLIFY = F)

#' Function to mark mitochondrial and ribosome genes with MT- and RB- labels
#' 
#' @param cells_x_genes path to cells_x_genes.genes.txt
#' @param mitoch - character vector with mitochondrial gene names
#' @param ribos - character vector with ribosomal gene names
#' 
#' @return Data.table of 1 column with ribosomal and mitochondrial gene names changed
#' from "MITOGENE" to "MT-MITOGENE", "RIBOGENE" to "RB-RIBOGENE"
#' 
#' @export mark_rb_mt_genes
#' @note Does not print result to console on some reason

mark_rb_mt_genes <- function(cells_x_genes, mitoch, ribos){
  genes <- fread(file = cells_x_genes, header = FALSE)
  if (nrow(genes[grepl("MT-", V1)])){
    warning(nrow(genes[grepl("MT-", V1)]), " mitochondrial genes seemed to be already replaced, skipping \n")
  }else{
    genes[V1 %in% mitoch, V1:=paste0("MT-", V1)]
    genes[V1 %in% mitoch_repeats$V1, V1:=paste0("MT-", V1)]
  }
  if (nrow(genes[grepl("RB-", V1)])){
    warning(nrow(genes[grepl("RB-", V1)]), " ribosomal genes seemed to be already replaced, skipping \n")
  }else{
    genes[V1 %in% ribos, V1:=paste0("RB-", V1)]
    genes[grepl("rRNA", V1), V1:=paste0("RB-", V1)]
  }
  return(genes)
}






#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions. Will be added to the next release
#' version of BUSpaRse.
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}


pca <- function(dgCMatrix){
  tmp <- prcomp_irlba(dgCMatrix, n = 2) # scales and centers by default
  tmp_df <- as.data.frame(tmp$x)
  return(tmp_df)
}
qc_pca <- Vectorize(pca, SIMPLIFY = F)


to_tibble <- function(dgCMatrix){
  tibl <- tibble(nCount = rowSums(dgCMatrix),
                 nGene = rowSums(dgCMatrix > 0),
                 id = deparse(substitute(dgCMatrix)))
  return(tibl)
}
qc_gene_num <- Vectorize(to_tibble, SIMPLIFY = F)

knee_filter <- function(dgCMatrix){
  bc_rank <- barcodeRanks(dgCMatrix, lower = 10)
  tot_counts <- Matrix::colSums(dgCMatrix)
  dgCMatrix_filtered <- dgCMatrix[, tot_counts > metadata(bc_rank)$inflection]
}


make_plots <- function(seu_object, path, top10_variable, good_markers){
  #Dotplots
  dp <- DotPlot(seu_object, features = top10_variable, cols = c("blue", "red")) +
    coord_flip()
  ggsave(file = paste(path, "top10_dot.png", sep = ""), plot = dp, dpi = 320)
  
  dp <- DotPlot(seu_object, features = unique(good_markers), cols = c("blue", "red")) +
    coord_flip()
  ggsave(file = paste(path, "good_marks_dot.png", sep = ""), plot = dp, dpi = 320)
  
  
  #Feature plots
  fp <- FeaturePlot(seu_object, top10_variable, cols = c("blue", "red"),
                    reduction="umap", order = T, ncol = 5, slot = "scale.data",
                    keep.scale = "all", repel = T,  label = T, 
                    interactive = F, combine = F, coord.fixed = T)
  for (i in fp){
    plot_name <- paste(path, "fp_top10_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  
  fp <- FeaturePlot(seu_object, good_markers, cols = c("blue", "red"),
                    reduction="umap", order = T, slot = "scale.data", 
                    keep.scale = "all", repel = T,  label = T, 
                    interactive = F, combine = F, coord.fixed = T)
  
  for (i in fp){
    plot_name <- paste(path, "fp_good_marks_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  #Violin plots
  
  vp <- VlnPlot(seu_object, slot = "scale.data", good_markers, pt.size = 0.01, combine = F)
  for (i in vp){
    plot_name <- paste(path, "vp_good_markers_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  vp <- VlnPlot(seu_morula, slot = "scale.data", top10_variable, pt.size = 0.01, combine = F)
  vp
  for (i in vp){
    plot_name <- paste(path, "vp_top10_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  VFP <- VariableFeaturePlot(seu_object, log = FALSE)
  VFP <- LabelPoints(plot = VFP, points = top10_variable, repel = TRUE)
  ggsave(paste(path, "VFP.png", sep = ""), plot = VFP, dpi = 320)
  #return (TRUE)
}