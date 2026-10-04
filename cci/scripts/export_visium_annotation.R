suppressPackageStartupMessages(library(Seurat))
load('data/raw/visium_mouse_cortex_annotated.RData')
x <- visium.brain
write.table(data.frame(barcode=colnames(x),celltype=as.character(Idents(x))), 'data/resources/visium_annotation.tsv',sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c(paste('spots',ncol(x)),paste('genes',nrow(x)),capture.output(table(Idents(x)))), 'results/visium_annotation_summary.txt')
cat('annotations:',ncol(x),'spots\n')
