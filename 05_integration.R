suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(ggplot2))
suppressWarnings(library(stringr))
suppressWarnings(library(ggrepel))
suppressWarnings(library(tidyverse))
suppressWarnings(library(cowplot))
suppressWarnings(library("optparse"))
suppressWarnings(library("harmony"))
})
options(future.globals.maxSize= 89128960000)



option_list = list(
make_option(c("--rds"), type="character", default=NA,
              help="RDS input, list of Seurat objects, eg after SCT and PCA", metavar="character"),
make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")
);

parser=OptionParser(option_list=option_list)
opt=parse_args(parser)

if(!is.na(opt$out)){
  outpath=(opt$out)
  if(str_sub(outpath,-1,-1) =='/'){
  outpath=substr(outpath,1,nchar(outpath)-1)
  }}


cd=getwd()
output=paste0('/', opt$out)
ifelse(!dir.exists(paste0(cd, output)), dir.create(paste0(cd, output)), FALSE)
ifelse(!dir.exists(paste0(cd, output, '/objects')), dir.create(paste0(cd, output, '/objects')), FALSE)

merged_seurat_sct=readRDS(opt$rds)
# For this huge dataset, clustering can take days
# so set only a minimal resolution
res <- c(0.3, 0.5)


integrated_per_pool_sct <- merged_seurat_sct %>% RunHarmony(group.by.vars =c("Pool"), plot_convergence = FALSE,reduction.use='pca', assay.use = "SCT")

rm(merged_seurat_sct)

print("Start Neighbors")
integrated_per_pool_sct <- FindNeighbors(integrated_per_pool_sct, reduction ='harmony', dims = 1:70,verbose = TRUE, k.param = 30)
print("FindClusters now")
integrated_per_pool_sct <- FindClusters(
  integrated_per_pool_sct, 
  resolution = res,
  verbose = TRUE, 
  algorithm=1, 
  method = "igraph",
  n.start = 5,
  n.iter = 5
  )
  #igraph to not have dense matrix, and 1 to use Louvain algorithm for speed
  # Reduce number of runs per resolution to n.start x n.iter = 25 instead of 100
print("UMAP Calculation now")
integrated_per_pool_sct <- RunUMAP(integrated_per_pool_sct, reduction = 'harmony', dims = 1:70,verbose = TRUE)
print("Saving rds now")

saveRDS(integrated_per_pool_sct, paste0(cd, output, '/objects/sct_integratedHarmony_k30_leiden.rds.gz'), compress = "gzip")
