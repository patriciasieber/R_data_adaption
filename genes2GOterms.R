## map gene ids to entrez id, map them to GO terms..

library("biomaRt")
library("readxl")

genes <- read.csv("/data/Postdoc/A_C3/deg_overlap_C_D3/deg_invitro_overlap_genes.csv",sep="\t")

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

goids <- getBM(attributes=c('external_gene_name','go_id', 'name_1006'),  
              filters='external_gene_name',values=genes$deg, 
              mart=ensembl)
goids <- goids[nchar(goids$go_id) > 0,]

gene2go <- lapply(as.character(genes$deg),function(x){
  current <- goids[goids$external_gene_name == x,]
  current_goid <- paste(current$go_id,collapse=";")
  current_godescr <- paste(current$name_1006,collapse=";")
  return(cbind(current_goid,current_godescr))
})

genes$go_id <- unlist(lapply(gene2go,function(x) return(x[1])))
genes$go_description <- unlist(lapply(gene2go,function(x) return(x[2])))


write.table(genes,"/data/Postdoc/A_C3/deg_overlap_C_D3/deg_invitro_overlap_genes2go.csv",sep="\t",row.names=F)


gene2go <- lapply(as.character(genes$deg),function(x){
  current <- goids[goids$external_gene_name == x,]
  ov <- as.character(genes[genes$deg == x,]$overlap)
  if(nrow(current) > 0){
    output <- data.frame(deg=rep(x,length(current$go_id)), overlap=ov, go_id=current$go_id, go_description=current$name_1006)
  }else{
    output <- data.frame(deg=x,overlap=ov,go_id="",go_description="")
  }
  return(output)
})
gene2go <- do.call("rbind",gene2go)

write.table(gene2go,"/data/Postdoc/A_C3/deg_overlap_C_D3/deg_invitro_overlap_genes2go.csv",sep="\t",row.names=F,quote=F)



## add information about innatedb 
genes <- read.csv("/data/Postdoc/A_C3/deg_overlap_C_D3/deg_invitro_overlap_genes.csv",sep="\t",stringsAsFactors=F)

innatedb <- read_excel("/data/Postdoc/A_C3/deg_overlap_C_D3/innatedb_curated_genes.xls")
innatedb_genes <- unique(innatedb$`Gene Symbol`)

in_innatedb <- sapply(genes$deg,function(x){ 
  current <- x %in% innatedb_genes
  if(current){
    return("x")
  }else{
    return(" ")
  } 
})

genes$involvedinInnateDB <- in_innatedb
write.table(genes,"/data/Postdoc/A_C3/deg_overlap_C_D3/deg_invitro_overlap_genes2InnateDB.csv",sep="\t",row.names=F,quote=F)
