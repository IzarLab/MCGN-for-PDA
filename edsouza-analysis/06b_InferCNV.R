#!/usr/bin/Rscript

#### Title: inferCNV with sample name as argument
#### Author: Jana Biermann, PhD. Modified by Edridge D'Souza

print(paste("Start:",Sys.time()))

library(dplyr)
library(Seurat)
library(infercnv)
library(stringr)
library(gplots)
library(ggplot2)
library(scales)
library(viridis)

pat <- commandArgs(trailingOnly=TRUE)[1] # 'GM18'

outdir <- file.path('/home/ubuntu/data/cxcr4-pdac/infercnv', pat)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=T)
}
message_ <- function(m) {
    message(paste0(Sys.time(), ': ', m))
}

seu_dir <- file.path('/home/ubuntu/data/cxcr4-pdac/seurat')
seu_file <- file.path(seu_dir, pat, paste0(pat,'_cb.rds'))
ref_file <- file.path(seu_dir, 'reference.RDS')
gene_order_file <- file.path('/home/ubuntu/InstallTemp/infercnv/refdata-gex-GRCh38-2020-A_gen_pos.txt')
seu_obj_out <- file.path(outdir, paste0(pat, '_infercnv.rds'))
pdf_file_ref <- file.path(outdir, paste0(pat, '_infercnv_withref.pdf'))
pdf_file_neverref <- file.path(outdir, paste0(pat, '_infercnv_noref.pdf'))

message_('Reading in data')
seu <- readRDS(seu_file)
ref <- readRDS(ref_file)
seu$is_ref <- 'not_ref'
ref$is_ref <- 'is_ref'

# counts matrix
seu <- merge(seu, y = ref)  
counts_matrix <- as.data.frame(GetAssayData(object = seu, slot = 'counts'))

# annotation
annot <- seu$is_ref %>% as.data.frame()
colnames(annot) <- c('cell')

# gene order
gene_order <- read.table(gene_order_file, header = F, row.names = 1)

# create the infercnv object
message_('Creating InferCNV object')
options(scipen = 100)
options(expressions=10000)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annot,
                                    gene_order_file=gene_order,
                                    ref_group_names = c('is_ref'))

# perform infercnv operations to reveal cnv signal
# Ran in about 4h on r5.4xlarge for sample GM36
message_('Running InferCNV')
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             denoise=T,
                             HMM=T,
                             analysis_mode='subclusters',
                             tumor_subcluster_partition_method='leiden',
                             leiden_resolution=0.01,
                             #tumor_subcluster_partition_method='random_trees',
                             output_format='pdf',
                             num_threads = parallel::detectCores()-4) # https://github.com/broadinstitute/infercnv/issues/446#issuecomment-1526827831
message_('Done running InferCNV')

# identify malignant cells
message_('Adding InferCNV metadata to seurat object')
seu <- add_to_seurat(seurat_obj = seu, 
                     assay_name = "RNA", 
                     infercnv_output_path = outdir)
cnv_cols<-grep('proportion_scaled_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$cnv_avg<-rowMeans(cnvs)
q10<- 0.1
seu$malignant<-ifelse(seu$cnv_avg > q10, 'malignant','non-malignant')

# add CNV metrics
cnv_cols<-grep('proportion_scaled_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$proportion_scaled_cnv_avg<-rowMeans(cnvs)

cnv_cols<-grep('proportion_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$proportion_cnv_avg<-rowMeans(cnvs)

cnv_cols<-grep('has_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$has_cnv_avg<-rowMeans(cnvs)


################################################################################
# Save JUST the non-reference cells
message_('Re-running PCA and UMAP before plotting')

seu_only <- readRDS(seu_file)
seu_only@meta.data <- seu@meta.data

seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu) 
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)

message_(paste0('Saving at', seu_obj_out))

saveRDS(seu_only, seu_obj_out)

################################################################################

message_(paste0('Plotting w/ never reference cells at', pdf_file_neverref))
pdf(pdf_file_neverref)
textplot(addmargins(table(seu_only$is_ref, seu_only$malignant)), cex=1.2,halign='left')
textplot(table(seu_only$celltype_bped_fine,seu_only$malignant), cex=0.9,halign='left')
hist(seu_only$cnv_avg,breaks=100, main='Average absolute CNV level; all cells',
     xlab = 'Average absolute CNV proportion')
abline(v=q10,col="red")
DimPlot(seu_only, reduction='pca',group.by='malignant')
DimPlot(seu_only, reduction='umap',group.by='malignant')
DimPlot(seu_only, reduction='umap',group.by='is_ref')

FeaturePlot(seu_only, features =  c("proportion_scaled_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)
FeaturePlot(seu_only, features =  c("proportion_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)

FeaturePlot(seu_only, features = c('MLANA',"MITF", 'PMEL', "MKI67"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)
FeaturePlot(seu_only, features = c('PTPRC',"CD8A", 'CD68','CD79A' ), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)
FeaturePlot(seu_only, features = c('ALB',"KRT15", 'COL1A1', "VWF"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

DimPlot(seu_only, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.position="none")

DimPlot(seu_only, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.position="none")
dev.off()


# For completeness, we'll include the reference cells
message_(paste0('Plotting with reference cells at', pdf_file_ref))
pdf(pdf_file_ref)
textplot(addmargins(table(seu$is_ref, seu$malignant)),cex=1.2,halign='left')
textplot(table(seu$celltype_bped_fine,seu$malignant),cex=0.9,halign='left')
hist(seu$cnv_avg,breaks=100, main='Average absolute CNV level; all cells',
     xlab = 'Average absolute CNV proportion')
abline(v=q10,col="red")
DimPlot(seu, reduction='pca',group.by='malignant')
DimPlot(seu, reduction='umap',group.by='malignant')
DimPlot(seu, reduction='umap',group.by='is_ref')

FeaturePlot(seu, features =  c("proportion_scaled_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)
FeaturePlot(seu, features =  c("proportion_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)

FeaturePlot(seu, features = c('MLANA',"MITF", 'PMEL', "MKI67"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)
FeaturePlot(seu, features = c('PTPRC',"CD8A", 'CD68','CD79A' ), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)
FeaturePlot(seu, features = c('ALB',"KRT15", 'COL1A1', "VWF"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.position="none")

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.position="none")
dev.off()

message_(paste('Finished with sample', pat))