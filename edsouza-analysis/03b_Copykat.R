# Code from https://github.com/navinlabcode/copykat/blob/master/vignettes/copycat-vignettes.Rmd
library(Seurat)
library(copykat)
library(dplyr)

outdir <- '/home/ubuntu/data/cxcr4-pdac/copykat'            

args <- commandArgs(trailingOnly=TRUE)
sample <- args[1] # 'GM18'

sample_outdir <- file.path(outdir, sample)
if (!dir.exists(sample_outdir)) {
    dir.create(sample_outdir, recursive=T)
}
samplepath <- paste0('/home/ubuntu/data/cxcr4-pdac/seurat/',
                     sample,'/', sample, '_cb.rds')
raw <- readRDS(samplepath)
exp.rawdata <- as.matrix(raw@assays$RNA@counts)
rm(raw); gc()

# Ran in 78 mins on sample of 3188 cells; some samples take ~4h (r5.large instance)
# Same sample ran in about 20 minutes on c5.4xlarge with 16 CPUs
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", 
                        cell.line="no", ngene.chr=5, 
                        win.size=25, KS.cut=0.2, 
                        sam.name=file.path(sample_outdir, sample), 
                        distance="euclidean", 
                        n.cores=parallel::detectCores())
saveRDS(copykat.test, file.path(outdir, sample, paste0(sample, '_copykat.rds')))

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)
print('Prediction test')
head(pred.test) %>% print()
print('CNA test')
head(CNA.test[ , 1:5]) %>% print()


### Heatmap 1
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]
pred <- pred[!is.na(pred)]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

message('Producing heatmap #1')
pdf(file.path(outdir, sample, paste0(sample, '_heatmap1.pdf')))
heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
dev.off()


# Heatmap 2
tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")] %>% sapply(function(x) gsub('-', '.', x))
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chrom <- cbind(CHR,CHR)
                                                                       
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = 'RdBu')))(n = 999)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
 

message('Producing heatmap 2')
pdf(file.path(outdir, sample, paste0(sample, '_heatmap2.pdf')))
heatmap.3(t(tumor.mat),
            dendrogram="r", 
            distfun = function(x) parallelDist::parDist(x, threads=4, method = "euclidean"), 
            hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chrom,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')

dev.off()
