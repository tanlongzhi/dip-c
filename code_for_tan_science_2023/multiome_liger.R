'''
LIGER
We ran http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html
with k=100, k=50, and k=150. 
'''
library(rliger)
library(dplyr)

# create liger object, run liger and cluster
a <- readRDS(file = "/path/to/combined-seurat.rds")
a <- seuratToLiger(a, combined.seurat = TRUE, meta.var = 'age') 
a <- normalize(a)
a <- selectGenes(a)
a <- scaleNotCenter(a)
a <- optimizeALS(a, k = 100) # k=100, k=50, and k=150
a <- quantile_norm(a)
a <- louvainCluster(a, resolution = 1)

# plot dataset and cluster
a <- runUMAP(a, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
all.plots <- plotByDatasetAndCluster(a, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
all.plots[[1]] + all.plots[[2]]

# plot gene
plotGene(a,"NEUROD1",pt.size=1,zero.color=rgb(0.85,0.85,0.85),cols.use=c(rgb(0.85,0.85,0.85),rgb(1,0,0)),plot.by="none",return.plots = TRUE) 
  + coord_fixed() + theme(legend.position = "none")

#rename clusters example
a@clusters <- recode(a@clusters,"31" = "microglia", "62" = "microglia")

# generate csv for percentage of cell types 
total <- a@cell.data %>% group_by(dataset) %>% summarise(count = n())
table <- total
features =  c("0.1a", "0.2", "0.7", "2.3", "0.1b", "0.4", "37.6")

for (cluster in levels(a@clusters)) {
  asub <- subsetLiger(a,clusters.use = cluster, remove.missing = FALSE)
  temp <- asub@cell.data %>% group_by(dataset) %>% summarise(count = n())
  for (feature in features){
    if (!feature %in% temp$dataset){
      print(cat(feature, cluster))
      temp <-rbind(temp, c(feature, 0))
    }
  }
  temp <- temp %>% arrange(dataset)
  table[, cluster] <- as.numeric(temp$count) / total$count
}
table
write.csv(table, "/path/to/percents.csv")

# factor analysis 
factor_markers <- getFactorMarkers(a)
factor_markers_2 <- getFactorMarkers(a,dataset1="2.3",dataset2="37.6")
factor_plots <- plotFactors(a, plot.tsne = TRUE) + coord_fixed()

# get cell type markers
wilcox.results <- runWilcoxon(a, compare.method = "clusters")
write.csv(wilcox.results,'liger_marker_genes.csv')

# differential genes 
asub <- subsetLiger(a, clusters.use = c("gc_1", "gc_2"))
cluster.results <- runWilcoxon(asub, compare.method = "clusters")

table <- cluster.results %>%
  group_by(group) %>%
  slice_max(n = 10, order_by = logFC) 
write.csv(table, "degenes_granule.csv")

# get gene loadings 
BiocManager::install('org.Hs.eg.db')
BiocManager::install('reactome.db')
BiocManager::install('fgsea')
gsea_results <- runGSEA(k)

gene_loadings <- plotGeneLoadings(k, do.spec.plot = FALSE, return.plots = TRUE)

gene_loadings[[73]][[1]] + coord_fixed() # get the color bar

gene_loadings[[73]][[1]] + coord_fixed() + theme(legend.position = "none") # granule stage 1
cat(names(sort(k@W[73,],decreasing=TRUE)[1:10]),sep="\n")
sort(k@W[73,],decreasing=TRUE)[1:50]
cat(names(sort(k@W[73,],decreasing=TRUE)[1:50]),sep="\n")
gsea_results[[73]][1:10,]
write.table(sort(k@W[73,],decreasing=TRUE)[1:100],"module_a.txt",col.names=FALSE,sep="\t", quote = FALSE)


