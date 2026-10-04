suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(future))
plan('sequential')
args <- commandArgs(trailingOnly=TRUE); out <- args[1]
if (dir.exists(out)) stop('Output exists')
dir.create(out,recursive=TRUE)
x <- readMM('data/baseline_inputs_cellchat/expression.mtx'); genes <- readLines('data/baseline_inputs_cellchat/genes.txt')
meta <- read.delim('data/baseline_inputs/metadata.tsv',row.names=1,check.names=FALSE)
rownames(x)<-genes;colnames(x)<-rownames(meta)
db <- subsetDB(CellChatDB.human,search=c('Secreted Signaling','Cell-Cell Contact'),key='annotation')
profiles<-list();scores<-list()
for (s in sort(unique(meta$sample))) {
 t<-proc.time()[3];ii<-meta$sample==s;m<-meta[ii,,drop=FALSE];m$celltype<-factor(m$celltype)
 obj<-createCellChat(x[,ii,drop=FALSE],meta=m,group.by='celltype');obj@DB<-db;obj<-subsetData(obj)
 obj<-identifyOverExpressedInteractions(obj,features=rownames(obj@data.signaling))
 obj<-computeCommunProb(obj,type='truncatedMean',trim=0,nboot=100,seed.use=20260928)
 z<-as.data.frame(as.table(obj@net$prob),stringsAsFactors=FALSE)
 names(z)<-c('sender','receiver','interaction','score')
 z$feature<-paste(z$interaction,paste(z$sender,z$receiver,sep='|'),sep='::');z$sample<-s
 scores[[s]]<-z[,c('feature','score','sample')]
 saveRDS(obj@net,file=file.path(out,paste0(s,'_net.rds')))
 profiles[[s]]<-data.frame(method='CellChat',sample=s,seconds=proc.time()[3]-t,n_scores=nrow(z))
 write.table(do.call(rbind,profiles),file.path(out,'profile.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(do.call(rbind,scores),gzfile(file.path(out,'scores.tsv.gz')),sep='\t',row.names=FALSE,quote=FALSE)
 cat(s,'seconds',proc.time()[3]-t,'\n')
}
jsonlite::write_json(list(status='DONE',method='CellChat',samples=length(scores),nboot=100,CellChat_version=as.character(packageVersion('CellChat')),notes='Native score, full measured LR feature set; truncatedMean trim=0; original probability not identical to product statistic'),file.path(out,'receipt.json'),auto_unbox=TRUE,pretty=TRUE)
