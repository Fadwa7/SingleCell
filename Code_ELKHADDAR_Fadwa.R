#######################  PROJET TRANSCRIPTOMIQUE ##################
###################################################################
# Autor 	: Fadwa EL KHDDAR 

#### DATASET 1

library(dplyr)
library(readxl)
library(Seurat)
library(destiny)
library(slingshot)


# Importation des fichiers

TAM_data <- Read10X(data.dir='SingleCell/')

# Initialisation de variables
TAM <- CreateSeuratObject(counts = TAM_data, project = "tumor-associated macrophages", min.cells = 10, min.features = 800)

# Les variables de l'objet
str(TAM)

# Récupération des noms de gènes 
gene.names <- rownames(TAM)
length(gene.names)

# Récupération des gènes mitochondriaux pour les éliminer dans l'étape du contrôle qualité 

mt.genes <- gene.names[grep("^mt-", gene.names)]


names(TAM@meta.data)

TAM[["percent.mt"]] <- PercentageFeatureSet(TAM, pattern = "^mt")
TAM[["percent.mt"]]


# Visualisation des données par VlnPlot

VlnPlot(TAM, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)


# Filter les données pour éliminer que celles dont le pourcentage de gènes mitochondriaux est sup à  12
# et UMIs sup à 65000
TAM<- subset(TAM, subset = percent.mt < 12 & nCount_RNA < 65000)


# Normalisation
TAM <- NormalizeData(TAM)

Normalized <- NormalizeData(TAM)

# Les gènes variables
Gene_variable <- FindVariableFeatures(TAM, selection.method = "vst", nfeatures=2000)

# GO des gènes de cycles cellulaire
GO <- read_excel("GO_term_summary_20221217_042818.xlsx")
GO <- data.frame(GO)

# Récupération des noms de gènes dans le tableau des gènes GO de la souris 
symbol <- GO$Symbol

# Elimination des gènes GO.
TAM<- Gene_variable[!rownames(Gene_variable) %in% symbol,]

top25 <- head(VariableFeatures(TAM), 25)
plot1 <- VariableFeaturePlot(TAM)
LabelPoints(plot=plot1, points=top25, repel=TRUE, ynudge = 0, xnudge = 0)


# Application du PCA

all.genes <- rownames(TAM)
TAM <- ScaleData(TAM, features=all.genes)
TAM <- RunPCA(TAM, features=VariableFeatures(TAM), ndims.print=1:3, nfeatures.print=10)
DimPlot(TAM, reduction="pca")
ElbowPlot(TAM)


# t-SNE
TAM <- RunTSNE(TAM, dims=1:35)
DimPlot(TAM, reduction="tsne")

# UMAP
TAM <- RunUMAP(TAM, dims=1:35)
DimPlot(TAM, reduction="umap")

# Clustering

TAM <- FindNeighbors(TAM, dims = 1:35) ##Pour les 35 premiers PCA
TAM <- FindClusters(TAM, resolution=1.3) ## J'ai dû réduire la résolution pour pouvoir récupérer 8 populations 
DimPlot(TAM, reduction="pca", label= T)
DimPlot(TAM, reduction="tsne", label=T)
DimPlot(TAM, reduction="umap", group.by = "seurat_clusters",label=T)


# Les gènes marqueurs 
TAM_markers <- FindAllMarkers(TAM, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
TAM_markers<- TAM_markers[TAM_markers$p_val_adj<0.01,]
a <- TAM_markers %>% group_by(cluster) %>% slice_max(n = 5 , order_by = avg_log2FC )
View(a)
# Annotation des clusters

new.cluster.ids <- c("Macrophage_0","Macrophage_1","Macrophage_2","Macrophage_3","Macrophage_4","Macrophage_5","Macrophage_6","Macrophage_7")
names(new.cluster.ids) <- levels(TAM)
TAM <- RenameIdents(TAM, new.cluster.ids)
DimPlot(TAM, reduction = "pca", label = TRUE, pt.size = 0.5) + NoLegend()


# La Diffusion map  pour une projection 35D de données 

Emb <- Embeddings(TAM, reduction = "pca")
dif_map <- DiffusionMap(Emb)


#Trajectoires et Pseudotime 

names_clusters <- Idents(TAM)
sling <- slingshot(dif_map@eigenvectors[,1:15], clusterLabels = names_clusters, reducedDim = "diffusion", start.clus = "Macrophage_0")
sds <- as.SlingshotDataSet(sling)


plot(dif_map@eigenvectors[,1],  dif_map@eigenvectors[,2], col = names_clusters, legend_main = "Cell state")

lines(sds, type = "c", lwd=3 )

legend(x="topright",col=names_clusters,pch=16,legend=levels(names_clusters))
