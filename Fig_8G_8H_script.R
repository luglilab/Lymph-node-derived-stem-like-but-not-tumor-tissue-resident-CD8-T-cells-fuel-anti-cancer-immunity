####################
####################
# FIG_8G script
####################
####################

R

set.seed(123)
library(edgeR)


# raw matrix
counts <- read.table(file="Dataset_1_raw_count_matrix.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t", row.names = 1)
# gene annotation
annotation <- read.table(file= "annotation.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t")

# "1" means Treated; "2" Controls
Treat <- factor(c(rep(1,4), rep(2,4)))
# DGE object
dge_data <- DGEList(counts, group = Treat)
# filter genes
keep <- filterByExpr(dge_data)
# lib size
dge_data <- dge_data[keep,,keep.lib.sizes=FALSE]
# normalize
dge_data < -calcNormFactors(dge_data)
dge_data$samples


# different donors
Sample <- factor(c("A","B","C","D","A","B","C","D"))

# design matrix
design <- model.matrix(~Sample+Treat)
dge_data <- estimateDisp(dge_data, design, robust = TRUE)

# CPM
matrix_cpm <- cpm(dge_data)
# log2CPM
matrix_log2cpm <- log2(matrice_cpm+1)

# DEGs
fit <- glmQLFit(dge_data,design)
qlf_contVStrat <- glmQLFTest(fit)
summary(decideTests(qlf_contVStrat))


topTags(qlf_contVStrat)
GeniDE_qlf_contVStrat <- topTags(qlf_contVStrat, n = dim(qlf_contVStrat)[1])
toptag <- GeniDE_qlf_contVStrat@.Data[[1]]

# DEGs FDR<0.05
DEG_FDR005 <- toptag[toptag$FDR<0.05,]
# add anno
DEG_FDR005 <- merge(DEG_FDR005, annotation, by.y="gene_ID", by.x="row.names", sort=FALSE)
DEG_FDR005 <- DEG_FDR005[, c("Row.names", "gene_Symbol", "gene_Function", "logFC", "logCPM", "F", "PValue", "FDR")]


# HEATMAP DEGs
library(gplots)

# DEGs FDR<0.05
DEGsign <- toptag[toptag$FDR<0.05,]
matrix_heat <- matrix_log2cpm[(rownames(DEGsign)),]
matrix_heat <- merge(matrix_heat, annotation, by.y="gene_ID", by.x="row.names", sort=FALSE)
row.names(matrix_heat) <- matrix_heat $gene_Symbol
column_to_delete <- c(1,10,11)
matrix_heat <- matrix_heat[,-column_to_delete]

# HEATMAP
heat.colors <- colorRampPalette(c("violet", "black", "gold"))(256)
# z-score
scale.data <- as.matrix((matrix_heat-apply(matrix_heat,1,mean))/apply(matrix_heat,1,sd))

pdf("FIG_8G.pdf", height=200, width=40)
heatmap.2(scale.data,Colv = T,
            dendrogram="both",scale="none",colsep=c(0),
            sepcolor="white",sepwidth=c(0.01),margins = c(13,13),
            col=heat.colors,trace="none",labRow = rownames(matrix_heat),labCol = colnames(matrix_heat),
            key=T,keysize=0.8,density.info="none",symkey=TRUE, key.xlab=c("log2ratio"),cexRow=0.2,ColSideColors=c(rep("blue",4),rep("red",4)))
            legend("topright",legend=c("Treated","Control"),col=c("blue", "red"),pch=c(15,15),cex=4,pt.cex=5,bty = "n")
par(cex.main=0.5)
dev.off()

####################
####################
# FIG_8H script
####################
####################

# rnk PreRanked GSEA; na_pos = Control, na_neg = Treated
all_degs_annotate <- merge(toptag, annotation, by.x = "row.names", by.y = "gene_ID", sort=FALSE)
input_GSEA <- cbind(all_degs_annotate$gene_Symbol, all_degs_annotate$logFC)
colnames(input_GSEA) <- c("GeneName","rank")
input_GSEA <- input_GSEA[order(as.numeric(input_GSEA[,2]),decreasing = TRUE),]
