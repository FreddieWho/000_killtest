suppressPackageStartupMessages(library(CellChat))
dir.create('data/resources', recursive=TRUE, showWarnings=FALSE)
for (species in c('human','mouse')) {
 db <- get(paste0('CellChatDB.',species))
 resolve <- function(x) {
   if (x %in% rownames(db$complex)) {
     z <- as.character(unlist(db$complex[x,,drop=FALSE])); paste(unique(z[!is.na(z) & z!='']),collapse='|')
   } else x
 }
 z <- db$interaction
 z$ligand_subunits <- vapply(z$ligand,resolve,character(1))
 z$receptor_subunits <- vapply(z$receptor,resolve,character(1))
 write.table(z, paste0('data/resources/cellchat_',species,'.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
}
writeLines(c(paste('CellChat',packageVersion('CellChat')),capture.output(sessionInfo())), 'results/R_environment.txt')
cat('Exported human and mouse resources\n')

z<-unique(as.character(unlist(CellChatDB.human$cofactor)));writeLines(z[!is.na(z)&z!=""],"data/resources/cellchat_cofactor_genes.txt")

z<-unique(as.character(unlist(CellChatDB.mouse$cofactor)));writeLines(z[!is.na(z)&z!=""],"data/resources/cellchat_mouse_cofactor_genes.txt")
