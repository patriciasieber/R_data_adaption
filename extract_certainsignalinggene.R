## extract genes corresponding to a certain signaling pathway

library(ggplot2)

dir <- "/data/"
setwd(dir)

sig_genes <- read.csv("part1/treatment_2_VS_control_signaling.csv",sep="\t",header=T)
gene_ids <- as.character(sig_genes$id)

degs <- read.csv("treatment_2_VS_control.csv",sep="\t",header=T)
degs_id <- as.character(degs$id)

sig_pos <- unlist(lapply(gene_ids,function(x){
  return(grep(x,degs_id))
}))
degs_sig <- degs[sig_pos,]

degs_sig$external_gene_name <- degs_sig$external_gene_name
degs_sig$wikigene_description <- degs_sig$wikigene_description
degs_sig <- degs_sig[,c(1,20:21,2:19)]
degs_sig <- unique(degs_sig)

## calculate max p value over all comparisons
degs_pval <- degs_sig[,grep("_pval",colnames(degs_sig))]
maxpval <- as.numeric(apply(degs_pval,1,max))
degs_sig$max_pval <- maxpval

## sort by max p value
degs_sig <- degs_sig[order(degs_sig$max_pval),]

write.table(degs_sig,"treatment_2_VS_control_signaling_part2.csv",sep="\t",col.names=T,row.names=F,quote=F)


## hist of log2fc
pdf("log2_fc_mrn_hist_part2.pdf")
hist(degs_sig$log2_fc_mrn,breaks=100)
dev.off()


## expression of genes of interest over time
candidates <- NULL # TODO: fill with gene of interest!
degs_candidates <- unlist(lapply(candidates,function(x){as.character(degs_sig[degs_sig$external_gene_name == x,]$id)}))


mrn <- read.csv("mrn.csv",sep="\t",header=T) 
mrn_candidates_pos <- unlist(lapply(degs_candidates,function(x){grep(x,rownames(mrn))}))
mrn_candidates <- mrn[mrn_candidates_pos,]

## gene1
mrn_gene1 <- data.frame(name=colnames(mrn_candidates),donor=c(rep("patient1",3),rep("patient2",3),rep("patient3",2),rep("patient4",2)), 
                       condition=rep(c("t0","t10","t60"),4), 
                       value=as.numeric(mrn_candidates[1,]))

pdf("gene1.pdf")
ggplot(data=mrn_gene1,aes(x=donor,y=value,fill=condition)) +
  geom_bar(stat="identity",position=position_dodge())
dev.off()




# mrn_candidates_rel <- mrn_candidates
# for(i in 1:nrow(mrn_candidates)){
#   current_row <- as.numeric(mrn_candidates[i,])
#   row_max <- max(current_row)
#   
#   rel_list <- unlist(lapply(current_row,function(x){
#     return(x/row_max)
#   }))
#   mrn_candidates_rel[i,] <- rel_list
# }



## add biomart name and description ####
all_degs <- read.csv("treatment_2_VS_control.csv",sep="\t",header=T)
deg_ids <- as.character(all_degs$id)
len_degs <- length(deg_ids)

# library(biomaRt)
# 
# if(interactive()){
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description",
#                      "name_1006","definition_1006"),mart=mart)    ## check available options with: listAttributes(mart)
# select(mart,keys=as.character(all_degs$id),
#        columns=c("ensembl_gene_id","external_gene_name", "description","name_1006","definition_1006"),
#        keytype="ensembl_gene_id")
# } # --> did not work

mart <- read.csv("mart_export.txt",header=T)
genename <- vector(length=len_degs)
genedescr <- vector(length=len_degs)

for(i in 1:len_degs){
  current_id <- deg_ids[i]
  current_mart <- mart[mart$Gene.stable.ID == current_id,]
  if(nrow(current_mart) == 1){
    genename[i] <- as.character(current_mart$Gene.name)
    genedescr[i] <- as.character(current_mart$WikiGene.description)
  }else{
    genename[i] <- ""
    genedescr[i] <- ""
  }
}

all_degs$gene_name <- genename
all_degs$description <- genedescr

all_degs <- all_degs[,c(1,20:21,2:19)]

write.table(all_degs,"treatment_2_VS_control_descr.csv",sep="\t",col.names=T,row.names=F,quote=F)


## pca only for genes of interest ####
##
cd56_genes <- read.csv("treatment_2_VS_control_signaling_part2.csv",sep="\t",header=T)
cd56_ids <- as.character(sig_genes$id)

mrn <- read.csv("mrn.csv",sep="\t",header=T) 
mrn_genes <- rownames(mrn)
sig_pos <- unlist(lapply(sig_ids,function(x)grep(x,mrn_genes)))
mrn_sig <- mrn[sig_pos,]
write.table(mrn_sig,"mrn_signaling_part2.csv",sep="\t",col.names=NA,row.names=T,quote=F)

counts <- read.csv("fc.counts.tsv",sep="\t",header=T)
counts_genes <- as.character(counts$X)
sig_pos <- unlist(lapply(sig_ids,function(x)grep(x,counts_genes)))
counts_sig <- counts[sig_pos,]
write.table(counts_sig,"counts_signaling_part2.csv",sep="\t",col.names=T,row.names=F,quote=F)

