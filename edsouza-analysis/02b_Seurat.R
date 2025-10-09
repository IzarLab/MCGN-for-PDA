#!/usr/bin/env Rscript

#### Title: Seurat analysis for CellBender output with sample name as argument
#### (includes TCR info integration and preliminary cell type annotation)
#### Author: Jana Biermann, PhD. Modified by Edridge D'Souza
#### This script takes the cellbender output plus the scrublet predictions 
#### and produces the processed Seurat object for each sample

print(paste('Start:',Sys.time()))

library(dplyr)
library(gplots)
library(ggplot2)
library(ggrastr)
library(DropletUtils)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(Matrix)
library(stringi)
library(Seurat)


minFeature<-200
maxFeature<- 12000
minCount<- 800
maxCount<- 70000
maxMT<-15


args <- commandArgs(trailingOnly=TRUE)
smpl <- args[1]
pat <- paste0('GM', smpl) 
message(paste(Sys.time(), 'Starting process on sample', pat))


### File definitions
file <- paste0('/home/ubuntu/data/cxcr4-pdac/cellbender_matrix/GM', smpl, '/GM', smpl, '_cellbender_filtered.h5')
doubletfile <- paste0('/home/ubuntu/data/cxcr4-pdac/scrublet/GM', smpl, '_doublets.txt')


outdir <- file.path('/home/ubuntu/data/cxcr4-pdac/seurat', pat)
if (!dir.exists(outdir)) {
       dir.create(outdir, FALSE, recursive=T)
}


seu_rds <- file.path(outdir, paste0(pat, '_cb.rds'))
outplot <- file.path(outdir, paste0(pat, '_plots.pdf'))


################################################################################


##### Loading, merging, QC, dimension reduction #####
### Load dataset
message(paste(Sys.time(), 'Loading datasets'))
seu.data <- Read10X_h5(file, 
                    use.names = TRUE, unique.features = TRUE)

### Initialize the Seurat object with the raw (non-normalized data)
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, 
                            min.cells = 1, min.features = 1)


################################################################################
message(paste(Sys.time(), 'Annotating MT genes and barcodes'))
# Annotate MT genes
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^MT-", col.name = "percent.mt")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPS", col.name = "percent.rps")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPL", col.name = "percent.rpl")
seu_raw$percent.rp=seu_raw$percent.rps + seu_raw$percent.rpl

# Annotate
seu_raw$sample <- pat
# seu_raw$patient<-stri_split_fixed(str = pat, pattern = "_", n = 2)[[1]][1]
# seu_raw$condition<-stri_split_fixed(str = pat, pattern = "_", n = 2)[[1]][2]
seu_raw$barcode <- rownames(seu_raw@meta.data)
# seu_raw$barcode_pat<-paste0(seu_raw$barcode_orig,'_',pat)

# Add clinical data. TODO: OBTAIN
# clin<-read.csv('data/clin.csv',na.strings = '')
# seu_raw@meta.data<-left_join(seu_raw@meta.data,clin,by='sample')
# rownames(seu_raw@meta.data)<-seu_raw$barcode_orig


################################################################################
# Identify doublets using scrublet
#doublet_rate<-doublet_rate[doublet_rate$sample==pat,2]
message(paste(Sys.time(), 'Adding doublet metadata'))
doublets <- read.table(doubletfile, header = T)
rownames(doublets) <- doublets$barcode
doublets$predicted_doublets <- as.logical(doublets$predicted_doublets)
seu_raw <- AddMetaData(seu_raw, doublets)

seu_raw <- NormalizeData(seu_raw)  # Need for cell cycle scoring
seu_raw <- CellCycleScoring(
			  object = seu_raw,
			  g2m.features = cc.genes$g2m.genes,
			  s.features = cc.genes$s.genes
			)
################################################################################
### subset 
message(paste(Sys.time(), 'Subsetting raw data'))
seu <- subset(seu_raw, 
            subset = nFeature_RNA > minFeature & 
                    nFeature_RNA < maxFeature & 
                    # nCount_RNA > minCount &   # Mori says this is redundant with nFeature
                    # nCount_RNA < maxCount & 
                    percent.mt < maxMT & 
                    (predicted_doublets == F |
                    Phase %in% c('G2M', 'S'))
              )


### Workflow RNA
message(paste(Sys.time(), 'Transforming and clustering'))


seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, vars.to.regress = "percent.mt") 


seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)



################################################################################
### Cell type identification with SingleR as first-pass
### These will be inherently inaccurate but allow us to compare w other methods
message(paste(Sys.time(), 'Celltype prediction'))
seu_sce <- as.SingleCellExperiment(seu)


# Blueprint Encode reference
bped<-BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels

pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels

# Human Primary Cell Atlas reference
hpca <- HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']]<-pred_hpca_main$pruned.labels

pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
pruneScores(pred_hpca_fine)
seu[['celltype_hpca_fine']]<-pred_hpca_fine$pruned.labels

################################################################################

### stats
message(paste(Sys.time(), 'Compiling summary stats'))
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_raw_features','n_raw_cells','n_predicted_doublets','n_features','n_cells','median_features','median_counts','cutoff_features','cutoff_counts','cutoff_mt')
rownames(stats)<-pat

stats$sample<-pat
stats$n_raw_features<-dim(seu_raw@assays$RNA@counts)[1]
stats$n_raw_cells<-dim(seu_raw@assays$RNA@counts)[2]
stats$n_predicted_doublets <-length(which(seu_raw@meta.data$predicted_doublets ==T))
stats$n_features<-dim(seu@assays$RNA@counts)[1]
stats$n_cells<-dim(seu@assays$RNA@counts)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))
stats$cutoff_features<-paste(minFeature,maxFeature)
stats$cutoff_counts<-paste(minCount,maxCount)
stats$cutoff_mt<-paste(maxMT)

### Save objects
saveRDS(seu, file = seu_rds)

################################################################################


### write pdf reports
message(paste(Sys.time(), 'Creating plots'))
# stats
pdf(file = outplot)
textplot(t(stats), cex=1.2 ,halign='left')

# plots raw data
ggplot(seu_raw@meta.data, aes(x=seu_raw$nCount_RNA,y = seu_raw$nFeature_RNA, col=seu_raw$percent.mt)) + 
rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="blue", high="green") + 
labs(color = "Percent MT") + theme_classic() + ggtitle('Raw object')

ggplot(seu_raw@meta.data, aes(x=seu_raw$nCount_RNA,y = seu_raw$nFeature_RNA, col=seu_raw$doublet_scores)) + 
rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="lightgrey", high="darkviolet") + 
labs(color = "doublet_scores") + theme_classic()+ ggtitle('Raw object')

print(VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','doublet_scores'),
            ncol = 3,group.by = 'sample',pt.size = 0))

# QC plot filtered
ggplot(seu@meta.data, aes(x=seu$nCount_RNA,y = seu$nFeature_RNA, col=seu$percent.mt)) + 
rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="blue", high="green") + 
labs(color = "Percent MT") + theme_classic()+ ggtitle('Filtered object')

ggplot(seu@meta.data, aes(x=seu$nCount_RNA,y = seu$nFeature_RNA, col=seu$doublet_scores)) + 
rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="lightgrey", high="darkviolet") + 
labs(color = "doublet_scores") + theme_classic()+ ggtitle('Filtered object')

print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','doublet_scores'),
            ncol = 3,
            # group.by = 'patient',
            pt.size = 0))

FeaturePlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','doublet_scores'), 
            min.cutoff = "q05", max.cutoff = 'q95',order=T, raster = T)

# PCA
print(ElbowPlot(seu))
DimPlot(seu, reduction = "pca",group.by = 'ident',raster = T,shuffle = T)

# UMAP
DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',raster = T,shuffle = T)



## bped
plotScoreHeatmap(pred_bped_fine, clusters=seu_sce@colData$ident, fontsize = 6, main='pred_bped_fine')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
theme(legend.text=element_text(size=6))


FeatureScatter(seu,feature1 ='nCount_RNA',feature2 = 'nFeature_RNA',shuffle = T,
            group.by = 'celltype_bped_main',raster = T)

VlnPlot(seu, features = c("nFeature_RNA"),group.by = 'celltype_bped_main',pt.size = 0)

## hpca
plotScoreHeatmap(pred_hpca_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_main')
DimPlot(seu, reduction = "umap",label = T, group.by = 'celltype_hpca_main',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
       guides(col = guide_legend(nrow = 30, override.aes = list(size=5))) +
       theme(legend.text=element_text(size=6))


###################################### Marker analysis
genes_t_cells <- c('BCL11B',   'FYN',   'IKZF1',   'CCL5',   'LYST',  'CD247',   'CBLB', 'IQGAP2',  'ETS1', 'SKAP1',  'ELMO1',  'ITGA4', 'DTDH1','AOAH',   'PTPRC',  'CD44',  'PYHIN1',  'CD8A',  'CD96', 'IFI6',  'NFAT5',  'THEMIS', 'PRKCH',  'ITK',  'PDE3B',  'BACH2','NCK2',  'CELF2',  'ITGA1', 'RBPJ',  'TNIK',   'LEF1',   'RIPOR2', 'AKT3', 'RORA', 'FOXP1',  'NCAM1',   'KLRC1',   'KLRC2','KLRD1',  'LYN',  'NCR1',   'KLRK1',  'RPLP1', 'RPL41','RPL13',  'ACTB',  'B2M', 'STAT4',   'RUNX2',   'ITGA4','ITGA1','CD69', 'ITGAE',   'CD44', 'CD58','IL2RA',  'TBC1D4',  'CTLA4','ICOS', 'CD4')
genes_myeloid <- c('CLEC9A','XCR1',  'CLNK', 'SLC38A1','RTN1','SIRPA','IDO1','LAMP3','CD200','MS4A2','KIT','MRC1','MERTK','FCGR1A','CD38','CD163L1','SELENOP','F13A1','DAB2','SIGLEC1','FTL','FTH1','NAV3','C3','P2RY12','VCAN','FCN1','LYZ')
genes_endothelial <- c('VEGFC', 'DLL4' ,'EFNB2' ,'RGCC', 'KIAA1217', 'ARHGAP18', 'MMRN1', 'FLT4', 'SEMA3D','RELN', 'IL1R1', 'NR2F2', 'CDH11', 'MKI67' )
genes_fibroblasts <- c('RGS5', 'CSPG4', 'ABCC9', 'PDGFRA', 'LUM', 'DCN', 'CTHRC1', 'COL1A1' ,'COL3A1', 'COL6A1' ,'LAMA2' )
genes_bcells <- c('BANK1','BLK','MKI67','TOP2A','SDC1','CD38','PRDM1', 'CD19', 'CD20')  # Added CD19 and CD20
genes_neuronal <- c('GFAP', 'AQP4', 'SLC1A2', 'GJA1', 'RELN', 'NRGN', 'STX1A', 'NEFH', 'SYP', 'SNRPN', 'SATB2', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'SLC6A1', 'MAG', 'MOG', 'MBP', 'MOBP', 'VCAN', 'CSPG4', 'SOX10', 'OLIG1' )
seu <- AddModuleScore(seu, features = list(genes_t_cells), name="T_cells")
seu <- AddModuleScore(seu, features = list(genes_myeloid), name="Myeloid")
seu <- AddModuleScore(seu, features = list(genes_endothelial), name="Endothelial")
seu <- AddModuleScore(seu, features = list(genes_fibroblasts), name="Fibroblasts")
seu <- AddModuleScore(seu, features = list(genes_bcells), name="B_cells")
seu <- AddModuleScore(seu, features = list(genes_neuronal), name="Neuronal")
FeaturePlot(seu,
            features = paste0(c('T_cells', 'B_cells', 'Myeloid'), 1), 
            order=T, raster = T, pt.size=5)
FeaturePlot(seu,
            features = paste0(c('Endothelial', 'Fibroblasts',  'Neuronal'), 1), 
            order=T, raster = T, pt.size=5)


dev.off()


message(paste(Sys.time(), 'Finished sample', pat))