#!/home/l.s.adonin/miniconda3/bin/Rscript
library(Seurat)
library(ggplot2)
setwd("/mnt/r840_storage/gentech/l.s.adonin/")
load("anchors.RData")

#Create reductions object from anchors
tmp <- lapply(anchors@object.list, function(x) x@reductions$pca)
reducts <- merge(tmp[[1]], tmp[-1])
rm(tmp)

all_integrated <- IntegrateEmbeddings(anchorset = anchors, reductions = reducts, dims = 1:50)

rm(anchors)
rm(reducts)

save(all_integrated, file = "/mnt/r840_storage/gentech/l.s.adonin/all_integrated.RData")
#Running Scaling data
all_integrated <- FindVariableFeatures(all_integrated, nfeatures = 3000)
all_integrated <- ScaleData(all_integrated)

#Running PCA reductions
all_integrated <- RunPCA(all_integrated, verbose = FALSE, npcs = 20) # uses HVG by default
ElbowPlot = ElbowPlot(all_integrated, ndims = 20)
PCAplot = PCAPlot(all_integrated)
ggsave("Elbowplot.svg", ElbowPlot) 
ggsave("PCAPlot.svg", PCAplot)
save(all_integrated, file = "/mnt/r840_storage/gentech/l.s.adonin/all_integrated.RData")

#Finding clusters
library(data.table)
library(MAST)
all_integrated <- FindNeighbors(all_integrated, dims = 1:10)
all_integrated <- FindClusters(all_integrated, resolution = 0.5)
save(all_integrated, file = "/mnt/r840_storage/gentech/l.s.adonin/all_integrated.RData")

PCA_clusters = PCAPlot(all_integrated)
ggsave("PCA_clusters.svg", PCA_clusters)

#Running UMAP
print("Starting UMAP")
all_integrated <- RunUMAP(all_integrated, dims = 1:15)
print("Saving UMAP reduction")
save(all_integrated, file = "/mnt/r840_storage/gentech/l.s.adonin/all_integrated.RData")

UMAP_plot = UMAPPlot(all_integrated)
ggsave("UMAP_plot.svg", UMAP_plot)

#Finding markers
markers <- FindAllMarkers(all_integrated, test.use = "MAST", only.pos = FALSE, 
                               min.pct = 0.25, logfc.threshold = 0.25)

fwrite(markers, file = "all_markers.tsv", sep = "\t")


