suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(stringr))
suppressWarnings(library(SoupX))
suppressWarnings(library("optparse"))
suppressWarnings(library(scDblFinder))
})

option_list = list(

make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices. Available only if complete/typical has been set", metavar="character"),

make_option(c("--raw-matrix"), type="character", default=NA,dest='raw',
              help="Path with input raw matrices. Needed when --remove-ambient is set to True", metavar="character"),

make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")

);


parser=OptionParser(option_list=option_list)
opt=parse_args(parser)


NM='vst'
DR='UMAP'

path=(opt$matrix)
if(str_sub(path,-1,-1) =='/'){path=substr(path,1,nchar(path)-1)}

rawpath=(opt$raw)
if(str_sub(rawpath,-1,-1) =='/'){rawpath=substr(rawpath,1,nchar(rawpath)-1)}


output=paste0('/', opt$out)
cd=getwd()
dir.create(paste0(cd, output), showWarnings = FALSE)
dir.create(paste0(cd, output, '/objects'), showWarnings = FALSE)

gex_list=list()
dbl_list=list()

#run SoupX for ambient correction
runSoupx_h5<- function(filt, rawh5, resolution, name, norm, DR){

    filt.matrix=Read10X_h5(filt,use.names = T)
    raw.matrix=Read10X_h5(rawh5,use.names = T)
   
 soup.channel  <- SoupChannel(raw.matrix$`Gene Expression`, filt.matrix$`Gene Expression`)
    srat=CreateSeuratObject(counts = filt.matrix$`Gene Expression`)
        if(norm=='SCT'){srat<-SCTransform(srat, verbose = F)}

        else{
             	srat <- NormalizeData(srat)
                srat <- FindVariableFeatures(srat, selection.method = "vst",nfeatures = 2000, verbose=FALSE)
                srat <- ScaleData(srat)
            }

            srat    <- RunPCA(srat, verbose = F, npcs=30)
            if(DR=='UMAP'){srat    <- RunUMAP(srat, dims = 1:30, verbose = F)}
            if(DR=='TSNE'){srat    <- RunTSNE(srat, dims = 1:30, verbose = F,check_duplicates = FALSE)}
            srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
            srat    <- FindClusters(srat, verbose = T, resolution = resolution)

    soup.channel  <- setClusters(soup.channel, setNames(srat@meta.data$seurat_clusters, rownames(srat@meta.data)))
    if(DR=='UMAP'){soup.channel  <- setDR(soup.channel, srat@reductions$umap@cell.embeddings)}
    if(DR=='TSNE'){soup.channel  <- setDR(soup.channel, srat@reductions$tsne@cell.embeddings)}
    soup.channel  <- autoEstCont(soup.channel,forceAccept=TRUE)
    adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

    return(adj.matrix)

    }

print('REMOVE AMBIENT')
  cd=getwd()
  dir.create(paste0(cd,output,'/SoupX'), showWarnings = TRUE)

   n=0
   for(file in list.files(path, pattern='h5')){
     n=n+1
     name=paste(strsplit(file, '_')[[1]][1],strsplit(file, '_')[[1]][2], strsplit(file, '_')[[1]][3], sep='_' )
     print(paste0("Removing ambient RNA from ", name))

     rawp=paste0(cd,'/',rawpath,'/',name,'_raw.h5')
     filteredp=paste0(cd,'/',path,'/',name,'_filtered.h5')

	
     adj.matrix <-runSoupx_h5(filteredp, rawp, 0.5, name,NM,DR)
      setwd(cd)
      fobj=CreateSeuratObject(counts=adj.matrix, assay="RNA")
      gex_list[[name]]=fobj
    }
dev.off()


print('IDENTIFY DOUBLETS')
print('First mild QC is applied removing only cells with less than 200 UMIs')
min_umi=200

for(name in names(gex_list)){
	current=gex_list[[name]]
	print(current)
	print(paste0("Removing cells with <", min_umi," UMIs from ", name))
	print(paste0("Removing doublets from ", name))

	current=subset(current, subset = nCount_RNA>min_umi)
	print(current)

	sce <- scDblFinder(GetAssayData(current, slot="counts"), clusters=FALSE)
	current$scDblFinder.score <- sce$scDblFinder.score
	current$scDblFinder.class <- sce$scDblFinder.class
	dbl_list[[name]]<-current
}

saveRDS(dbl_list, paste0(cd, output, '/objects/preprocessed_list.rds.gz'), compress = "gzip")
