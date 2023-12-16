library(data.table)
library(anndata)
library(Matrix)
library(dplyr)
library(Seurat)
library(data.table)
library(scAnnotatR)
library(tidyverse)
library(RCurl)
library(cowplot)
library(scPOP)

#SEURAT CLUSTERING


#convert anndata to single cell exp then seurat object
sce <- importAnnData(sampleDirs = "C:/Users/SarahSavvy/Documents/CMU PhD/Comp_Gen/Project/Data_files",
                     sampleNames = 'my_ann')
sce
        
rownames(sce) <- rowData(sce)$gene_id 
rowData(altExp(sce))
ann_seurat <- as.Seurat(sce, counts = "counts", data = "counts")
        
ann_seurat[["percent.mt"]] <- PercentageFeatureSet(ann_seurat, pattern = "^MT-")
        # Visualize QC metrics as a violin plot
VlnPlot(ann_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)
        #filtration
        
filtered_seurat <- subset(x = ann_seurat, 
                                  subset= (nFeature_originalexp >= 200))
        
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
        
# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 3
        
# Only keeping those genes expressed in more than 3 cells
filtered_counts <- counts[keep_genes, ]
        
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
        
sp <- NormalizeData(object = filtered_seurat)
        
        
sp <- FindVariableFeatures(object = sp)
sp <- ScaleData(object = sp)
sp <- RunPCA(object = sp)
sp <- FindNeighbors(object = sp)
        
sp <- FindClusters(object = sp)
sp <- RunTSNE(object = sp)
sp <- RunUMAP(object = sp, dims = 1:10, reduction = "pca")
DimPlot(object = sp, reduction = "tsne")
        
png("seurat_clusters.png")
DimPlot(object = sp, reduction = "umap")
        dev.off()
        
DefaultAssay(sp)
FeaturePlot(object = sp, features =  "Cd68")
row.names(sp@assays[["RNA"]]@data)
# Violin and Ridge plots
VlnPlot(object =sp, feature = c("Cd68","Aif1"))
        
VlnPlot(object =sp, feature = c("Slc39a12","Vstm2a"))
        
VlnPlot(object =sp, feature = c("Ccr2","Cd4"))
        
VlnPlot(object =sp, feature = c("Cd19","Cd22"))
        
VlnPlot(object =sp, feature = c("Plp1","Ermn"))
        
        
VlnPlot(object =sp, feature = c("Fabp7","Sept4"))
        
VlnPlot(object =sp, feature = c("Bsg","Cntn1"))
        
VlnPlot(object =sp, feature = c("Pdgfrb","Rgs5"))
        
VlnPlot(object =sp, feature = c("Cd34"))
        
        
VlnPlot(object =sp, feature = c("Cd2","Ncr1"))
        
VlnPlot(object =sp, feature = c("Mmp8","Cxcr2"))
        
VlnPlot(object =sp, feature = c("Itgb7","Cd69"))
        
VlnPlot(object =sp, feature = c("Cd34","Cd117"))
        
VlnPlot(object =sp, feature = c("Eps15","Fbxl7"))
        
VlnPlot(object =sp, feature = c("Itgam","Acap2"))
        
        
# Rename all identities
seurat_renamed <- RenameIdents(object = sp, 
                                       "0" = "microglial cell",
                                       "1" = "B cell",
                                       "2" = "oligodendrocyte",
                                       "3" = "natural killer cell",
                                       "4" = "Unknown",
                                       "5" = "mesenchymal stem cell",
                                       "6" = "neuron/ERCC genes",
                                       "7" = "neutrophil",
                                       "8" = "endothelial",
                                       "9" = "B cell",
                                       "10" = "T cell",
                                       "11" = "neuron/Bergmann glial cell",
                                       "12" = "mesenchymal stem cell/hematopoietic stem cell",
                                       "13" = "myeloid cell",
                                       "14" = "astrocyte",
                                       "15" = "B cell",
                                       "16" = "oligodendrocyte precursor cell",
                                       "17" = "neuron/Bergmann glial cell",
                                       "18" = "microglial cell",
                                       "19" = "granulocyte",
                                       "20" = "brain pericyte/ smooth muscle cell",
                                       "21" = "hematopoeitic stem cell",
                                       "22" = "oligodendrocyte")
png("annotated_seurat_clusters_no_subset.png")
        
        
seurat_renamed$CellType_new <- Idents(seurat_renamed)
annotation_seurat <- seurat_renamed@meta.data
        
DimPlot(object = seurat_renamed, 
                reduction = "umap", 
                label = TRUE,
                pt.size = 0.9,
                label.size = 5,
                repel = TRUE)
dev.off()
        

#benchmarking
#to run this you need the annotation from desc and seurat saved as a csv file or dataframe
DESC_annotated <- read.csv("annotations_DESC.csv")
        
ari(DESC_annotated$cell_ontology_class, DESC_annotated$Cell_type)
nmi(DESC_annotated$cell_ontology_class, DESC_annotated$Cell_type)
#seurat
ari(annotation_seurat$cell_ontology_class, annotation_seurat$CellType_new)
nmi(annotation_seurat$cell_ontology_class, annotation_seurat$CellType_new)


#bar plot
#create a dataframe of the scores 
df <- data.frame(Metric  = c("ARI", "NMI"),
                         Score = c("73.2", "77.9")
        )
        
df2 <- data.frame(Metric  = c("ARI", "NMI"),
                          Score = c("75.7", "75.6")
        )
        
# Basic barplot
ggplot(data=df2, aes(x=Metric, y=Score)) +
          geom_bar(stat="identity", fill = "steelblue")
        
#create data frame
df_merge <- data.frame(Algorithm= rep(c('Seurat', 'DESC'), each = 2),
                               Metric=rep(c('ARI', 'NMI'), times = 2),
                               Score =c(75.7, 75.6,73.2,77.9))
        
ggplot(df_merge, aes(x=Algorithm, y= Score, fill = Metric)) + 
          geom_bar(stat="identity", position=position_dodge(), width = 0.5)+
          theme_minimal()+
          labs(x='Algorithm', y='Score', title='ARI and NMI based performance of Seurat and DESC')+
          theme(plot.title = element_text(hjust=0.5, size=20, face='bold'), 
                axis.text = element_text(size = 15)) +
          scale_fill_manual('Position', values=c('coral2', 'steelblue'))
    
dev.print(file="scores.png", device=png, width=800)
        