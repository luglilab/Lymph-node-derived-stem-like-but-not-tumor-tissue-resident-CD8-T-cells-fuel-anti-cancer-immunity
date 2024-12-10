####################
####################
# FIG_15B script
####################
####################

R
set.seed(123)
library(Seurat)
library(dplyr)
library(qpcR)

cell_object <- readRDS("/GSE176022/ser.integrate.rds")
dim(cell_object)
#  22223 560920


# filter only patients with TCR info
cell_object_use <- subset(cell_object, subset = imid %in% c("MD01-004","MD01-005","MD043-011"))
dim(cell_object_use)
#  22223 247162

# filter only cells from LN and tumor
cell_object_use_2 <- subset(cell_object_use, subset = tissue %in% c("LN","tumor"))
dim(cell_object_use_2)
#  22223 159263


# filter only CD8 cells (cluster 0, 1, 3, 8, 9, 11, 13, 14)
CD8_cells <- subset(cell_object_use_2, subset = seurat_clusters %in% c("0","1","3","8","9","11","13","14"))
dim(CD8_cells)
# 22223 90450

# extract metadata
metadata <- CD8_cells@meta.data
dim(metadata)
# 90450    69

# extract counts
matrice <- CD8_cells@assays$RNA@counts
dim(matrice)
# 22223 90450

# re-create object to filter genes not expressed
CD8_cells <- CreateSeuratObject(counts = matrice, min.cells = 3, min.features = 200)
dim(CD8_cells)
# 19621 90450

CD8_cells <- AddMetaData(CD8_cells, metadata = metadata)


library(scRepertoire)
S1 <- read.csv("/MD01-004_LN_1.vdj/filtered_contig_annotations.csv")
S2 <- read.csv("/MD01-004_LN_2.vdj/filtered_contig_annotations.csv")
S3 <- read.csv("/MD01-004_LN_3.vdj/filtered_contig_annotations.csv")
S4 <- read.csv("/MD01-004_tumor_1.vdj/filtered_contig_annotations.csv")

S5 <- read.csv("/MD043-011_LN_1.vdj/filtered_contig_annotations.csv")
S6 <- read.csv("/MD043-011_LN_2.vdj/filtered_contig_annotations.csv")
S7 <- read.csv("/MD043-011_LN_3.vdj/filtered_contig_annotations.csv")
S8 <- read.csv("/MD043-011_tumor_2.vdj/filtered_contig_annotations.csv")
S9 <- read.csv("/MD043-011_tumor_3.vdj/filtered_contig_annotations.csv")
S10 <- read.csv("/MD043-011_tumor_4.vdj/filtered_contig_annotations.csv")
S11 <- read.csv("/MD043-011_tumor_5.vdj/filtered_contig_annotations.csv")

S12 <- read.csv("/MD01-005_LN_1.vdj/filtered_contig_annotations.csv")
S13 <- read.csv("/MD01-005_LN_2.vdj/filtered_contig_annotations.csv")
S14 <- read.csv("/MD01-005_LN_3.vdj/filtered_contig_annotations.csv")
S15 <- read.csv("/MD01-005_LN_4.vdj/filtered_contig_annotations.csv")
S16 <- read.csv("/MD01-005_tumor_2.vdj/filtered_contig_annotations.csv")
S17 <- read.csv("/MD01-005_tumor_3.vdj/filtered_contig_annotations.csv")
S18 <- read.csv("/MD01-005_tumor_4.vdj/filtered_contig_annotations.csv")
S19 <- read.csv("/MD01-005_tumor_5.vdj/filtered_contig_annotations.csv")
S20 <- read.csv("/MD01-005_tumor_6.vdj/filtered_contig_annotations.csv")
S21 <- read.csv("/MD01-005_tumor_7.vdj/filtered_contig_annotations.csv")
S22 <- read.csv("/MD01-005_tumor_8.vdj/filtered_contig_annotations.csv")
S23 <- read.csv("/MD01-005_tumor_9.vdj/filtered_contig_annotations.csv")



contig_list <- list(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23)


combined <- combineTCR(contig_list, 
                samples = c("MD01-004", "MD01-004", "MD01-004", "MD01-004", "MD043-011", "MD043-011", "MD043-011", "MD043-011", "MD043-011", "MD043-011", "MD043-011", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005", "MD01-005"), 
                ID = c(":LN-1", ":LN-2", ":LN-3", ":tumor-1", ":LN-1", ":LN-2", ":LN-3", ":tumor-2", ":tumor-3", ":tumor-4", ":tumor-5", ":LN-1", ":LN-2", ":LN-3", ":LN-4", ":tumor-2", ":tumor-3", ":tumor-4", ":tumor-5", ":tumor-6", ":tumor-7", ":tumor-8", ":tumor-9"),
                cells = "T-AB",
                removeNA = TRUE,
                removeMulti = TRUE)
  

# rename cells : con _:

head(x = colnames(x = CD8_cells))

a <-  colnames(CD8_cells)
a <- gsub(x = a, pattern = ":", replacement = "_:", fixed = TRUE)

CD8_cells <- RenameCells(object = CD8_cells, new.names = a)
head(x = colnames(x = CD8_cells))



seurat <- combineExpression(combined, CD8_cells, 
                  cloneCall="aa", 
                  proportion = FALSE, 
                  cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

tcr_meta <- seurat@meta.data[,79:84]


# mana's tcr
meta_mana <- subset(tcr_meta, grepl("CASSDLAGLNYGYTF|CASSLAGLAGPLAKNIQYF|CASSSISTGELFF|CASSLSGPVEDEQYF|CAWSEGRGNEKLFF|CASSYNNEQFF|CASSLPGDPYEQYF|CASSVDGGQSDEQYF|CASSRGEGGSNQPQHF|CASSSGPSYEQYF|CASSPLPSGGTYEQYF|CAWKEVRGQFF|CASSFWGNTEAFF|CATTGGQNTEAFF|CASSEVQGASNEKLFF|CASQSGILPWEQFF|CAISEWRAGSTDTQYF|CASSYLVGAGSNYGYTF|CASSFAGGGNQPQHF|CASSLAPGNEQFF|CASSLTGGSEAFF|CASSFQQNTEAFF|CAISGGTGANYGYTF|CASSLYTGELFF|CASSRSGSSYNEQFF|CATSDLRDRGNNEQFF|CASSVDWGAEAFF|CASRPYGTYGYTF|CASSEALGPGNTIYF|CASSLGANTGELFF|CASSWTGNQPQHF|CASSYPSAAAYNEQFF|CASSLETTGANVLTF|CASSLSGSSYNEQFF|CASSSGLRNIQYF|CASIPTGQNYGYTF|CASSIGTGSKPQHF|CASSLGDDSMNTEAFF|CASTPSAGANQPQHF|CASSEQGFWNGYTF|CASRPGQRYNSPLHF|CASGGTDTQYF|CAWESSRDIDDPEAFF|CASNKLGYQPQHF|CASSLTGGYTGELFF|CASSLLENQPQHF|CASNGEAETQYF|CASSDLAGLNYGYTF|CASSLAGLAGPLAKNIQYF|CASSLDPYEQYF", CTaa))
dim(meta_mana)
# 991   6


# filter mana cells
#create cell names as metadata colum
seurat[["CellName"]] <- colnames(seurat)
#subset the required cells
seurat_mana_cells <- subset(seurat, subset = CellName %in% rownames(meta_mana))
dim(seurat_mana_cells)
# 19621  991

table(seurat_mana_cells$imid)
#  MD01-004  MD01-005 MD043-011 
#        92        29       870 

matrice <- seurat_mana_cells@assays$RNA@counts

# re-create object to filter genes not expressed
mana_cells <- CreateSeuratObject(counts = matrice, min.cells = 3, min.features = 200)
dim(mana_cells)
# 11573   991

mana_cells <- AddMetaData(mana_cells, metadata = seurat_mana_cells@meta.data)

mana_cells <- NormalizeData(mana_cells)

mana_cells <- FindVariableFeatures(mana_cells, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(mana_cells)
mana_cells <- ScaleData(mana_cells, features = all.genes)

mana_cells <- RunPCA(mana_cells, features = VariableFeatures(object = mana_cells))


mana_cells <- JackStraw(mana_cells, num.replicate = 100)
mana_cells <- ScoreJackStraw(mana_cells, dims = 1:20)

mana_cells <- FindNeighbors(mana_cells, dims = 1:9)

####################
# Analyze different resolutions
####################

mana_cells <- FindClusters(mana_cells, resolution = seq(from=0, to=1, by=0.1), print.output = 0, save.SNN = T)

library(clustree)

clustree(mana_cells)

mana_cells <- FindClusters(mana_cells, resolution = 0.2)

mana_cells <- RunUMAP(mana_cells, dims = 1:9)

p1 <- DimPlot(mana_cells, reduction = "umap", group.by = "RNA_snn_res.0.2", label = TRUE)
p2 <- DimPlot(mana_cells, reduction = "umap", group.by = "tissue", label = TRUE)

pdf("Fig_15B_Right.pdf")
p1 
dev.off()


pdf("Fig_15B_Left.pdf", width = 23)
p2 
dev.off()


table(mana_cells$RNA_snn_res.0.2)
#   0   1   2   3 
# 493 340 124  34 

####### markers
 
cluster_markers <- FindAllMarkers(mana_cells, only.pos = TRUE, return.thresh = 0.05) 
cluster_markers_fdr001 <- cluster_markers[cluster_markers$p_val_adj <= 0.01,]
dim(cluster_markers_fdr001)
# 585    7