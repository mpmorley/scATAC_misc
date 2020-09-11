devtools::install_github("timoast/signac", ref = "develop")
detach('package:Signac',unload = T)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(clustree)
#library(scExtras)
library(TFBSTools)
library(JASPAR2018)
################################################################################################
#    9/10/2020 MPM Revised this script using Signac version 1.0
#
scriptVersion=2
##################################################################################################





#### Species specific packages 
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


##################################
source('~/dsdata/projects/R_Packages/scExtras/R/markerGenes_human.R')



geneDB <- EnsDb.Hsapiens.v86
genomeName<-'hg38'
genome <- BSgenome.Hsapiens.UCSC.hg38


set.seed(1234)
######################################
#     Plot used for marker genes. 
#
########################################
genePlot<- function(g,cols=NULL,plotdir){
  try({
    p1 <- DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + NoLegend()
    p2 <- DimPlot(rna, label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + NoLegend()
    p3 <- FeaturePlot(atac,g, max.cutoff = 'q95',order = T)
    p4 <- FeaturePlot(rna,g,order=T)
    (p1|p2)/(p3|p4) + plot_annotation(title = g,theme = theme(plot.title = element_text(size = 18)))
    ggsave(filename = paste0(plotdir,'/',g,'_UMAP_Vln.png'), width = 8,height=8)
  })
}  




################################################################
#
# RNA-seq filename if it exits
#####################################

rnaseqfile <- '../scRNA/Seurat/AGIQ424_RUL_D_CD45neg_20190922.RDS'


#### Setup in and out dirs
indir<- 'cellranger/'
outdir<-'./Signac'
dir.create(outdir)
plotdir <- paste0(outdir,'/plots')
dir.create(plotdir)
qcdir <-paste0(outdir,'/qc')
dir.create(qcdir)

### Create the Seurat Object

counts <- Read10X_h5(paste0(indir,"/filtered_peak_bc_matrix.h5"))

metadata <- read.csv(paste0(indir,"/singlecell.csv"), header = TRUE, row.names = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seqinfo(genome),
  fragments = paste0(indir,'/fragments.tsv.gz'),
  min.cells = 10,
  min.features = 200
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)



annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atac) <- annotations


##########################################################################
#
# QC steps and plotting
#
##########################################################################
atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(object = atac, fast = FALSE)


# add blacklist ratio and fraction of reads in peaks
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$total * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

atac$high.tss <- ifelse(atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(atac, group.by = 'high.tss') + NoLegend()

atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = atac, group.by = 'nucleosome_group')



VlnPlot(
  object = atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

ggsave(filename = paste0(qcdir,'/QC_Plot.png'))



atac <- subset(
  x = atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
atac
##########################################################################
#                     
#  Dimension reduction
#
###########################################################################

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(object = atac)


#Here we see there is a very strong correlation between the first LSI component and the total number of counts for the cell, 
#so we will perform downstream steps without this component.
DepthCor(atac)

dims <- 2:25


atac <- RunUMAP(
  object = atac,
  reduction = 'lsi',
  dims = dims,n.neighbors = 30
)

atac <- FindNeighbors(
  object = atac,
  reduction = 'lsi',
  k.param = 30,
  dims = dims
)
#### Matching what was used for RNA-seq
atac <- FindClusters(
  object = atac,
  algorithm = 1,
  resolution = c(.2,.4,.6,.8)
)

clustree(atac, prefix = "peaks_snn_res.")

DimPlot(object = atac, label = TRUE) + NoLegend()
##########################################################################
#
#             Run Motif and ChromVar analysis 
#
##########################################################################

pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(collection='CORE',all_versions = FALSE,tax_group='vertebrates')
)

motif.matrix <- CreateMotifMatrix(
  features = granges(atac),
  pwm = pfm,
  genome = genome
)

motif <- CreateMotifObject(data = motif.matrix, 
                           pwm = pfm)


atac <- SetAssayData(
  object = atac,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)


atac <- RunChromVAR(
  object = atac,
  genome = genome
)




##########################################################################
#                     Compute Gene Act Data
#   
##########################################################################


gene.activities <- GeneActivity(atac)

atac[['RNA']] <- CreateAssayObject(counts = gene.activities)

atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)


saveRDS(atac,paste0(outdir,'/ATAC.RDS'))
atac <- readRDS(paste0(outdir,'/ATAC.RDS'))



############################################################## 
# Integrate with scRNA 
###############################################################


rna <- readRDS(rnaseqfile)
rna <- FindVariableFeatures(
  object = rna,
  nfeatures = 4000
)


DefaultAssay(atac) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  reduction = 'cca',
  dims = 1:50
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$var_cluster,
  weight.reduction = atac[['lsi']],
  dims = dims
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)



plot1 <- DimPlot(rna, group.by = 'var_cluster', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(atac, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2

##Plot Marker genes
purrr::map(markergenes,genePlot,plotdir = plotdir,cols=NULL)



####################################################
#   Can filter the ATAC data for poorly predicted 
#  cells. 
#
####################################################


filterStats <- table(atac$prediction.score.max > 0.5)
ggplot(atac@meta.data,aes(x=prediction.score.max)) + 
  geom_histogram() + 
  theme_bw() + 
  geom_vline(xintercept = 0.5,color='red',linetype='dashed') + ggtitle(table(atac$prediction.score.max > 0.5)) +
  ggtitle(paste0('Keep=',filterStats[2], ' Remove=',filterStats[1]))
ggsave(filename = paste0(plotdir,'/ATAC_RNA_PredictionScoreHistogram.png'))




atac.filtered <- subset(atac, subset = prediction.score.max > 0.5)
atac.filtered$predicted.id <- factor(atac.filtered$predicted.id, levels = levels(rna))  # to make the colors match




DefaultAssay(atac.filtered) <- 'RNA'

atac.filtered$predicted.id <- droplevels(atac.filtered$predicted.id)
atac.filtered$predicted.id <- as.numeric(levels(atac.filtered$predicted.id))[atac.filtered$predicted.id]

table(atac.filtered$predicted.id)
atacpallette <- cpallette[sort(unique(as.numeric(atac.filtered$predicted.id)+1))]


saveRDS(atac.filtered,paste0(outdir,'/ATAC_filtered.RDS'))



