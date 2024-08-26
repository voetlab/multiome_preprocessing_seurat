suppressPackageStartupMessages({library(Seurat)
suppressWarnings(library(ggplot2))
suppressWarnings(library(stringr))
suppressWarnings(library(ggrepel))
suppressWarnings(library(tidyverse))
suppressWarnings(library(cowplot))
suppressWarnings(library("optparse"))
suppressWarnings(library("openxlsx"))
})



 
option_list = list(
make_option(c("--rds"), type="character", default=NA,
              help="RDS input, list of Seurat objects, eg after removal of ambient RNA/rename procedure", metavar="character"),
make_option(c("--QC"), type="character", default=NA, dest='QC',
              help="List of comma separated QC thresholds, in the order: min_feauture,max_feature,min_count, max_count, MT%, RB%. By default the MAD is estimated and everything that is larger or lower than 3x MAD is removed. The MT% is set at 5%.", metavar="character"),
make_option(c("--exclude"), type="character",dest='exclude', default=NA,
	      help="Multiomes to exclude", metavar="vector"),
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
ifelse(!dir.exists(paste0(cd, output, '/QC')), dir.create(paste0(cd, output, '/QC')), FALSE)
ifelse(!dir.exists(paste0(cd, output, '/objects')), dir.create(paste0(cd, output, '/objects')), FALSE)

gex_qc_list=c()

#Get QC metrics 
QC=as.numeric(strsplit(opt$QC, ',')[[1]])

rename_list=readRDS(opt$rds)


seurat_qc<- function(sobject,name,folder,nFeature_RNA_min,nFeature_RNA_max,percent_mito,percent_ribo,
                     nCount_RNA_min, nCount_RNA_max){

  sobject[["percent.mt"]] <- PercentageFeatureSet(sobject, pattern = "^MT-")
  sobject[["percent.rb"]] <- PercentageFeatureSet(sobject, pattern = "^RP[SL]")

  counts_df<- data.frame('Read_Counts'=sobject@meta.data['nCount_RNA'], 'Metric'='Read_Counts','Filter'='OG')
  colnames(counts_df)<- c('Read_Counts','Metric','Filter')
  genes_df<- data.frame('Gene_Counts'=sobject@meta.data['nFeature_RNA'],'Metric'='Gene_Counts','Filter'='OG')
  colnames(genes_df)<- c('Gene_Counts','Metric','Filter')
  mito<- data.frame("percent.mt"=sobject@meta.data['percent.mt'],'Metric'= 'MT%','Filter'='OG')


  sobject <- subset(sobject, subset = nFeature_RNA>nFeature_RNA_min & nFeature_RNA<nFeature_RNA_max & percent.mt<percent_mito & percent.rb< percent_ribo &
                      nCount_RNA < nCount_RNA_max & nCount_RNA > nCount_RNA_min )
  counts_df_qc<- data.frame('Read_Counts'=sobject@meta.data['nCount_RNA'], 'Metric'='Read_Counts','Filter'='QC')
  colnames(counts_df_qc)<- c('Read_Counts','Metric','Filter')
  genes_df_qc<- data.frame('Gene_Counts'=sobject@meta.data['nFeature_RNA'],'Metric'='Gene_Counts','Filter'='QC')
  colnames(genes_df_qc)<- c('Gene_Counts','Metric','Filter')
  mito_qc<- data.frame("percent.mt"=sobject@meta.data['percent.mt'],'Metric'= 'MT%','Filter'='QC')

  counts_df<-rbind(counts_df,counts_df_qc)
  genes_df<-rbind(genes_df,genes_df_qc)
  mito<- rbind(mito,mito_qc)

  p1<-ggplot(counts_df, aes(x=Metric, y=Read_Counts)) + geom_violin(fill='firebrick4') +
    theme_bw() + facet_grid(cols  = vars(Metric), rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5)

  p2<-ggplot(genes_df, aes(x=Metric, y=Gene_Counts)) + geom_violin(fill='turquoise4') +
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5 )

  p3<-ggplot(mito, aes(x=Metric, y=percent.mt)) + geom_violin(fill='forestgreen') +
    theme_bw() + facet_grid(cols  = vars(Metric),rows = vars(Filter), scales = "free")+
    geom_dotplot(binaxis='y',stackdir='center',alpha=0.2, color='gray81', dotsize=0.5)

  setEPS()
  postscript(paste0(folder,name,".eps"))
  plot(p1+p2+p3)
  dev.off()
  return(sobject)
}


#Check Values for QC and save
df<-data.frame()
	for(i in names(rename_list)){
		current<-rename_list[[i]]
		current[["percent.mt"]] <- PercentageFeatureSet(current, pattern = "^MT-")
		current[["percent.rb"]] <- PercentageFeatureSet(current, pattern = "^RP[SL]")

		quant_tbl <- data.frame(
		"nCount_RNA_10quantile" = quantile(x = current$nCount_RNA, probs = 10/100),
		"nCount_RNA_95quantile" = quantile(x = current$nCount_RNA, probs = 95/100),
		"nfeat_RNA_10quantile" = quantile(x = current$nFeature_RNA, probs = 10/100),
		"nfeat_RNA_95quantile" = quantile(x = current$nFeature_RNA, probs = 95/100),
		"percent.mt_95quantile" = quantile(x = current$percent.mt, probs = 95/100),
		"percent.rb_95quantile" = quantile(x = current$percent.rb, probs = 95/100)
		)
		rownames(quant_tbl)<-i
		df<-rbind(df,quant_tbl)
				}
write.xlsx(df, paste0(cd,output,'/QC/','QC_percentiles.xlsx'), sep='\t', row.names = T, col.names = T)


#Set exclude values
if(!is.na(opt$exclude)){exclude= as.vector(strsplit(opt$exclude, ',')[[1]])}
if(is.na(opt$exclude)){
        opt$exclude=c()
        exclude=opt$exclude}

print(paste0("MOs to exclude: ", exclude))

rename_list<- rename_list[which(!names(rename_list) %in% exclude)]

##Perform QC
nam=c()
lg=c()
hg=c()
lr=c()
hr=c()
hm=c()
hrib=c()

if(!is.na(opt$QC)){
for(i in names(rename_list)){
	 current<-rename_list[[i]]
	  current[["percent.mt"]] <- PercentageFeatureSet(current, pattern = "^MT-")
	  current[["percent.rb"]] <- PercentageFeatureSet(current, pattern = "^RP[SL]")


	nam=c(nam,i)
	low_genes=length(which(current@meta.data$nFeature_RNA<QC[1]))
	lg=c(lg,low_genes)
	high_genes=length(which(current@meta.data$nFeature_RNA>QC[2]))
	hg=c(hg,high_genes)
	low_rna=length(which(current@meta.data$nCount_RNA<QC[3]))
	lr=c(lr,low_rna)
	high_rna=length(which(current@meta.data$nCount_RNA>QC[4]))
	hr=c(hr,high_rna)
	high_mito=length(which(current@meta.data$percent.mt>QC[5]))
	hm=c(hm,high_mito)
	high_ribo=length(which(current@meta.data$percent.rb>QC[6]))
	hrib=c(hrib,high_ribo)

			}
qc_data=data.frame('Experiment'= nam, 'Low_Gene_Count'= lg, 'High_Gene_Count'=hg, 'Low_RNA_Count'=lr, 'High_RNA_Count'=hr, 'High_MT%'=hm, 'High_RB%'=hrib)
write.xlsx(qc_data, paste0(cd,output,'/QC/','QC_overview.xlsx'), sep='\t', row.names = T, col.names = T)
}



rown=c()
cbqc=c()
caqc=c()

if(!is.na(opt$QC)){
for(i in names(rename_list)){
	  print('QC performed of user-specified thresholds')
	  print(paste0(i))
	  current<-seurat_qc(rename_list[[i]],name=i,folder=paste0(cd, output, '/QC/'), nFeature_RNA_min = QC[1], nFeature_RNA_max = QC[2], nCount_RNA_min = QC[3], nCount_RNA_max = QC[4],percent_mito = QC[5],percent_ribo=QC[6])
	  rown=c(rown,i)
	  cbqc=c(cbqc,ncol(rename_list[[i]]))
	  caqc=c(caqc,ncol(current))   
	  print(paste0(i,' Cells before QC: ', ncol(rename_list[[i]]) ))
	  print(paste0(i,' Cells after QC: ', ncol(current) ))
          gex_qc_list[[i]]<- current
                      }
summary=data.frame('dataset'=rown, 'Cells_Before_QC'=cbqc, 'Cells_after_QC'=caqc)
write.xlsx(summary, paste0(cd,output,'/QC/','summary.xlsx'), sep='\t', row.names = T, col.names = T)
saveRDS(gex_qc_list, paste0(cd, output, '/objects/QC_list.rds.gz'), compress = "gzip")			
}
