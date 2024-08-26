suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(ggplot2))
suppressWarnings(library(stringr))
suppressWarnings(library(ggrepel))
suppressWarnings(library(tidyverse))
suppressWarnings(library(cowplot))
suppressWarnings(library("optparse"))
})



 
option_list = list(

make_option(c("--matrix"), type="character", default=NA,
              help="Path with input matrices.", metavar="character"),
make_option(c("--rds"), type="character", default=NA,
              help="RDS input, list of Seurat objects, eg after removal of ambient RNA", metavar="character"),
make_option(c("--out-folder"), type="character", default=NA,dest='out',
              help="Name of folder where the output will be saved", metavar="character")
); 

parser=OptionParser(option_list=option_list) 
opt=parse_args(parser)

if(!is.na(opt$matrix)){
path=(opt$matrix)
if(str_sub(path,-1,-1) =='/'){
path=substr(path,1,nchar(path)-1)
}}

if(!is.na(opt$out)){
  outpath=(opt$out)
  if(str_sub(outpath,-1,-1) =='/'){
  outpath=substr(outpath,1,nchar(outpath)-1)
  }}

 
cd=getwd()
output=paste0('/', opt$out)
ifelse(!dir.exists(paste0(cd, output)), dir.create(paste0(cd, output)), FALSE)
ifelse(!dir.exists(paste0(cd, output, '/objects')), dir.create(paste0(cd, output, '/objects')), FALSE)



gex_list=list()

if(is.na(opt$matrix)){gex_list=readRDS(opt$rds)}

if(!is.na(opt$matrix)){
for(file in list.files(path)) {
        print(file)
        name=paste(strsplit(file, '_')[[1]][1],strsplit(file, '_')[[1]][2], strsplit(file, '_')[[1]][3], sep='_' )
        gex_list[[name]]<- Read10X_h5(paste0(path,'/',file))
				}
        gex_list= lapply(X=gex_list, FUN=function(x){
        x=CreateSeuratObject(counts=x$`Gene Expression`, assay="RNA", min.cells=3) })
			}

##RENAME CELLS
##Add Pool, Multiome (MO) and Condition based on file name, e.g. MOC15_P_1A
rename_list=list()
for(i in names(gex_list)){
  current=gex_list[[i]]
  current@meta.data['Pool']<- i
  current@meta.data['MO']<- strsplit(i,'_')[[1]][1]

  if(strsplit(i,'_')[[1]][2] =='P' ){
      current@meta.data['Condition']<- 'PD'}

  if(strsplit(i,'_')[[1]][2] =='C' ){
      current@meta.data['Condition']<- 'Control'}

  if( (strsplit(i,'_')[[1]][2] !='C') &  (strsplit(i,'_')[[1]][2] !='P') ){
      current@meta.data['Condition']<- 'Mixed'}

 current<- RenameCells(current, new.names=paste(i,rownames(current@meta.data), sep="_"))
  rename_list[[i]]<- current
			  }
#generate a merged file to check QC metrics separately
if(length(names(rename_list)) > 1){
rename_list[['merged']]<- merge(rename_list[[1]], y=c(rename_list[-1]))

}

saveRDS(rename_list, paste0(cd, output, '/objects/renamed_objects.rds.gz'), compress = "gzip")
