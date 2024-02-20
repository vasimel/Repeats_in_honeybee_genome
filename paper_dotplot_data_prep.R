#' Get data for Dot Plot
#'
#' @description 
#' `get_plot_data` returns the data.table with pct.exp and mean.exp
#' split by Seurat cluster numbers and caste feature. 
#'
#' @param seu_object The SeuratObject to extract data from. Should undergo 
#' the UMAP/TSNE reduction and contain clusters feature.
#' 
#' @param Markers_dt The 2 column data.table. First column - aliases for genes/gene groups
#' Second column - transcript names, delimited by semi-colons
#' 
#' @param name_caste character with the name of the caste. Should match the seu_object$caste feature
#' 
#' 
#' @returns Outputs 5-column data.table with following columns: 
#' - cluster - factor number of cluster from the seu_object (factor)
#' - transcript - character name of the gene group (matches the Markers_dt[,1])
#' - caste - character matches name_caste provided
#' - pct.exp - numeric percent of cells in cluster expressing given transcript
#' - mean.exp - numeric mean number of transcripts per cell in the given cluster
#' 
#' 
#' @examples
#'
#' ## don't run this in calls to 'example(get_plot_data)'
#' \dontrun{
#' #small sample dataset from Seurat
#' example_seu <- pbmc_small
#' #Randomly assign caste
#' example_seu$caste <- sample(c("worker", "queen"), ncol(example_seu), replace = TRUE)
#' Markers_dt <- data.table(V1 = c("MS4A1", "CD79B", "CD79A", "HLA-DRA", "TCL1A", "HLA-DQB1", "FAKE"),
#'                          V2 = c("MS4A1", "CD79B", "CD79A", "HLA-DRA", "TCL1A", "HLA-DQB1", "FAKE"))
#'
#' DP_data <- get_plot_data(example_seu, Markers_dt, "queen")}
#' str(DP_data)

library(Seurat)
library(data.table)
library(Matrix)
get_plot_data <- function(Seu_object, Markers_dt, name_caste){

  my_markers <- unlist(strsplit(Markers_dt$V2,split = ";"))
  #This line with subsetting seems to be the slowest. 
  #The possible workaround is to extract only necessary data to Matrix, for example
  subset_seu <- subset(Seu_object, subset = `caste` == name_caste) #If you want to change metafeature to select on, it's harcoded here
  plotting_data <- subset_seu@assays$RNA@data
  plotting_data <- plotting_data[(rownames(plotting_data) %in% my_markers),]
  absent_markers <- my_markers[!(my_markers %in% rownames(plotting_data))]
  #adding empty values for absent markers
  absent_matrix <- Matrix(ncol = ncol(plotting_data),
                          nrow = length(absent_markers), 
                          data = 0)
  colnames(absent_matrix) <- colnames(plotting_data)
  rownames(absent_matrix) <- absent_markers
  #Merging absent_matrix and plotting_data together
  plotting_data <- rbind(plotting_data, absent_matrix)
  rm(absent_matrix, absent_markers)
  
  
  #Creating DT for renaming
  rename_dt <- data.table(transcripts = strsplit(Markers_dt$V2, ";"),
                          group_names = Markers_dt$V1)
  
  rename_dt$transcripts <- sapply(rename_dt$transcripts, unlist)
  
  
  #Grouping transcripts 
  rename_transcripts <- function(transcripts, group_name, mtx_in = plotting_data){
    ifelse(length(transcripts) == 1,
           mtx_data <- plotting_data[rownames(plotting_data) %in% transcripts, ],  
           mtx_data <- apply(plotting_data[rownames(plotting_data) %in% transcripts, ], 2, sum))
    
    new_mtx <- Matrix(nrow = 1,
                      ncol = ncol(plotting_data), 
                      data = mtx_data, sparse = T)
    rownames(new_mtx) <- group_name
    colnames(new_mtx) <- colnames(plotting_data)
    return(new_mtx)
  }
  
  plotting_data <- mapply(rename_transcripts, rename_dt$transcripts, rename_dt$group_names,
                          SIMPLIFY = T, USE.NAMES = F)
  
  #rbinding together all tables yielded from plotting_data()
  plotting_data <- do.call(rbind, plotting_data)
  
  barcode <- colnames(plotting_data)
  transcripts <- rownames(plotting_data)
  plotting_data <- t(plotting_data)
  plotting_data <- as.data.table(plotting_data)
  names(plotting_data) <- transcripts
  plotting_data$barcode <- barcode
  rm(barcode, transcripts)
  plotting_data$cluster <- subset_seu$seurat_clusters
  plotting_data <- melt.data.table(plotting_data, id.vars= c("barcode", "cluster"),  variable.name = "transcript" , variable.factor = FALSE)
  plotting_data$caste <- name_caste #please rename meta.feature if needed
  
  
  plotting_data[, ncells := .N, by = cluster]
  
  #working_version
  plotting_data[value > 0, mean.exp := mean(value, na.rm=T), by = .(transcript, cluster)]
  
  PercentAbove <- function(x){
    return(length(x = x[x > 0]) / length(x = x))
  }
  
  #alternative version of calculating pct.exp
  #plotting_data[, pct.exp := (n.exp/ncells)*100] 
  
  plotting_data[, pct.exp := PercentAbove(value)*100, by = .(transcript, cluster)]
  
  plotting_data <- plotting_data[,.(cluster, transcript, caste, mean.exp, pct.exp)]
  plotting_data[is.na(plotting_data)] <- 0
  plotting_data[, max.pct := max(pct.exp), by = .(transcript, cluster)]
  plotting_data[, max.exp := max(mean.exp), by = .(transcript, cluster)]
  # 
  plotting_data$mean.exp <- NULL
  plotting_data$pct.exp <- NULL
  # 
  setnames(plotting_data, c("max.pct", "max.exp"), c("pct.exp", "mean.exp"))
  plotting_data <- unique(plotting_data)
  return(plotting_data)
}
