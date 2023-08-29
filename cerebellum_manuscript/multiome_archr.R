'''
ArCHR
Most commands follow https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html.
'''
library(ArchR)
library(ggplot2)
library(cowplot)
addArchRGenome("hg38")
addArchRThreads(1)
library(parallel) 

# Process samples individually 
ArrowFiles <- c("/path/to/arrow")
proj <- ArchRProject(ArrowFiles)

seRNA <- import10xFeatureMatrix(
  input = c(".path/to/filtered_feature_bc_matrix.h5"),
  names = c("p014")
)

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)
proj <- proj[proj$TSSEnrichment > 2 & proj$nFrags > 2500 & !is.na(proj$Gex_nUMI) & proj$Gex_nUMI > 500] #filter

proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)

pdf("p014_combinedumap.pdf", width = 30, height = 10)
p1 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters3", embedding = "UMAP_Combined")
pp1 <- list(p1, p2, p3)
plot_grid(plotlist=pp1, ncol=3, align='h')
dev.off()

#rename clusters example 
labelNew =  c("gc_1", "gc_2", "gc_1", "gc_3", "gc_3", "gc_3", "gc_3", "gc_2", "gc_3", "gc_3" , "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19")
names = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19")
proj$Clusters_named <- mapLabels(proj$Clusters_res1.5, newLabels = labelNew, oldLabels = names)

# get marker genes list 
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
saveRDS(markerList, file = "/path/to/p014_markerList") 


# differential peaks for granule cells 
subset_df <- markerList[c("gc_1", "gc_2", "gc_3")]
write.csv(markerList$gc_3, file = "/path/to/markerlist.csv", row.names = FALSE)
write.csv(markerList$gc_1, file = "markerList_gc_1.csv", row.names = FALSE)
markerList <- readRDS("0.1a_3stagegranule_markerlist.rds")
markerList$gc_2$marker_id <- paste(markerList$gc_2$seqnames, "_", markerList$gc_2$start,"-",  markerList$gc_2$end, sep="")
markerList$gc_1$marker_id <- paste(markerList$gc_1$seqnames, "_", markerList$gc_1$start,"-",  markerList$gc_1$end, sep="")
markerList$gc_3$marker_id <- paste(markerList$gc_3$seqnames, "_", markerList$gc_3$start,"-",  markerList$gc_3$end, sep="")
saveRDS(markerList, file = "0.1a_3stagegranule_markerlist.rds")


#plot genes 
markerGenes  <- c("BOC", "FOXP2", "DDAH2", "CHRM3", "ERBB4", "RBFOX1")
p <- plotEmbedding(proj, colorBy="GeneExpressionMatrix", name=markerGenes, embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
p <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 0, color = "white") +
    theme(plot.margin = unit(c(0, 0, .5, .5), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

# plot liger cell type annotations onto archr 
library(rliger)
library(ArchR)

b <- readRDS(file = "~/research/multiome/liger/ligerobj.rds") 
proj <- loadArchRProject(path = "~/research/multiome/archr500/archr_0.7", force = FALSE, showLogo = TRUE)

#make the names the same
names(a@clusters) <- sub('_', '#', names(a@clusters)) 
names(a@clusters) <- sub('04', '0.4', names(a@clusters)) 
names(a@clusters) <- sub('01a', '0.1a', names(a@clusters)) 
names(a@clusters) <- sub('01b', '0.1b', names(a@clusters)) 
names(a@clusters) <- sub('02', '0.2', names(a@clusters)) 
names(a@clusters) <- sub('07', '0.7', names(a@clusters)) 
names(a@clusters) <- sub('376', '37.6', names(a@clusters)) 
names(a@clusters) <- sub('23', '2.3', names(a@clusters)) 

#create a dataframe with rownames = cells, col1 = cluster
x <- data.frame(row.names= as.vector(names(a@clusters)) , val= as.vector((a@clusters)))
common_names <- intersect(rownames(proj@cellColData), rownames(x))
clusters <- x[common_names, , drop = FALSE]
missing_names <- setdiff(rownames(cell_data), rownames(clusters))
missing_values <- rep("na", each = length(missing_names))
all_names <- c(rownames(clusters), missing_names)
all_values <- c(clusters[, 1], missing_values)
new_data <- data.frame(ClustersValue = all_values, row.names = all_names)
sorted_cell_data <- cell_data[order(row.names(cell_data)), ]
sorted_new_data <- new_data[order(row.names(new_data)), , drop = FALSE]
final_clusters <- sorted_new_data[, 1]
proj$Clusters3 <- final_clusters


# trajectory analysis for 0.1a
a <- loadArchRProject(path = "archr/archr_0.1a/")
plotEmbedding(a, colorBy="cellColData", name="Clusters2", embedding = "UMAP_Combined")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="SNAP25", embedding="UMAP_Combined", imputeWeights = NULL, plotAs="points")

# granule trajectory
a <- addTrajectory(a, name = "granule", groupBy = "Clusters2", trajectory = c("gc_1", "gc_2", "gc_3"), embedding = "UMAP_Combined", force = TRUE, spar = 2)
plotTrajectory(a, size=1, embedding = "UMAP_Combined", trajectory = "granule", colorBy = "cellColData", name = "granule", plotAs="points")[[1]]

trajGEM <- getTrajectory(a, name = "granule", useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)
plotTrajectoryHeatmap(trajGEM, pal = c(rgb(0.85,0.85,0.85),rgb(1,0,0)))
heatmapGEM <- plotTrajectoryHeatmap(trajGEM, returnMatrix = TRUE)
dim(heatmapGEM)
which(rownames(heatmapGEM)=="chr7:FOXP2")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="FOXP2", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="GRIA2", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="CHRM3", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="SYNE1", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")


trajPM <- getTrajectory(a, name = "granule", useMatrix = "PeakMatrix", log2Norm = FALSE)
plotTrajectoryHeatmap(trajPM, pal = c(rgb(0.85,0.85,0.85),rgb(1,0,0)))
heatmapPM <- plotTrajectoryHeatmap(trajPM, returnMatrix = TRUE)
dim(heatmapPM)
a@peakSet$name <- paste0(seqnames(a@peakSet),"_",ranges(a@peakSet))
a <- addPeakMatrix(a, force = TRUE)
plotEmbedding(a, colorBy="PeakMatrix", name="chr16_54930458-54930958", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="PeakMatrix", name="chr10_129964114-129964614", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="PeakMatrix", name="chr1_2050049-2050549", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="PeakMatrix", name="chr1_59814558-59815058", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")


devtools::install_github("GreenleafLab/chromVARmotifs")
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
a <- addMotifAnnotations(a, motifSet = "cisbp", name = "Motif")
a <- addBgdPeaks(a)
a <- addDeviationsMatrix(a, peakAnnotation = "Motif", force = TRUE)

trajMM <- getTrajectory(a, name = "granule", useMatrix = "MotifMatrix", log2Norm = FALSE)
plotTrajectoryHeatmap(trajMM, pal =c(rgb(0.85,0.85,0.85),rgb(1,0,0)))
heatmapMM <- plotTrajectoryHeatmap(trajMM, returnMatrix = TRUE)
dim(heatmapMM)
plotEmbedding(a, colorBy="MotifMatrix", name="z:NHLH2_80", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="MotifMatrix", name="z:RFX8_730", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="MotifMatrix", name="z:ZEB1_157", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="MotifMatrix", name="z:NEUROG1_87", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")

# MLI trajectory
a <- addTrajectory(a, name = "mli", groupBy = "Clusters2", trajectory = c("interneurons"), embedding = "UMAP_Combined", force = TRUE, spar = 2)
plotTrajectory(a, size=1, embedding = "UMAP_Combined", trajectory = "mli", colorBy = "cellColData", name = "mli", plotAs="points")[[1]]

trajGEM <- getTrajectory(a, name = "mli", useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)
plotTrajectoryHeatmap(trajGEM, pal = c(rgb(0.85,0.85,0.85),rgb(1,0,0)), labelTop=10)
heatmapGEM <- plotTrajectoryHeatmap(trajGEM, returnMatrix = TRUE)
dim(heatmapGEM)
which(rownames(heatmapGEM)=="chr7:FOXP2")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="NFIB", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="NRXN3", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="RORA", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="GRID2", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")


trajPM <- getTrajectory(a, name = "mli", useMatrix = "PeakMatrix", log2Norm = FALSE)
plotTrajectoryHeatmap(trajPM, pal = c(rgb(0.85,0.85,0.85),rgb(1,0,0)),labelTop=20)
heatmapPM <- plotTrajectoryHeatmap(trajPM, returnMatrix = TRUE)
dim(heatmapPM)
a@peakSet$name <- paste0(seqnames(a@peakSet),"_",ranges(a@peakSet))
a <- addPeakMatrix(a, force = TRUE)
plotEmbedding(a, colorBy="PeakMatrix", name="chr16_54930458-54930958", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="PeakMatrix", name="chr10_129964114-129964614", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="PeakMatrix", name="chr1_2050049-2050549", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="PeakMatrix", name="chr10_132537306-132537806", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")


devtools::install_github("GreenleafLab/chromVARmotifs")
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
a <- addMotifAnnotations(a, motifSet = "cisbp", name = "Motif")
a <- addBgdPeaks(a)
a <- addDeviationsMatrix(a, peakAnnotation = "Motif", force = TRUE)

trajMM <- getTrajectory(a, name = "mli", useMatrix = "MotifMatrix", log2Norm = FALSE)
plotTrajectoryHeatmap(trajMM, pal =c(rgb(0.85,0.85,0.85),rgb(1,0,0)))
heatmapMM <- plotTrajectoryHeatmap(trajMM, returnMatrix = TRUE)
dim(heatmapMM)
plotEmbedding(a, colorBy="MotifMatrix", name="z:SOX4_754", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="MotifMatrix", name="z:RORA_658", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="MotifMatrix", name="z:SNAI2_161", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="MotifMatrix", name="z:ESRRA_691", embedding="UMAP_Combined", size=1, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")

# plot markers for all ages
a <- loadArchRProject(path = "archr/archr_0.1a/")
a <- loadArchRProject(path = "archr/archr_0.1b/")
a <- loadArchRProject(path = "archr/archr_0.2/")
a <- loadArchRProject(path = "archr/archr_0.4/")
a <- loadArchRProject(path = "archr/archr_0.7/")
a <- loadArchRProject(path = "archr/archr_2.3/")
a <- loadArchRProject(path = "archr/archr_37.6/")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="FOXP2", embedding="UMAP_Combined", size=2, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="CHRM3", embedding="UMAP_Combined", size=2, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="RBFOX1", embedding="UMAP_Combined", size=2, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")
plotEmbedding(a, colorBy="GeneExpressionMatrix", name="SNAP25", embedding="UMAP_Combined", size=2, pal=c(rgb(0.85,0.85,0.85),rgb(1,0,0)), imputeWeights = NULL, plotAs="points")




