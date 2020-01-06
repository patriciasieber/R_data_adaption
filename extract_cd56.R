## extract genes corresponding to CD56 signaling pathway
library(ggplot2)

dir <- "/data/Postdoc/X_G/"
setwd(dir)

cd56_genes <- read.csv("part1/treatment_2_VS_control_CD56signaling.csv",sep="\t",header=T)
gene_ids <- as.character(cd56_genes$id)

degs <- read.csv("treatment_2_VS_control.csv",sep="\t",header=T)
degs_id <- as.character(degs$id)

cd56_pos <- unlist(lapply(gene_ids,function(x){
  return(grep(x,degs_id))
}))
degs_cd56 <- degs[cd56_pos,]

degs_cd56$external_gene_name <- cd56_genes$external_gene_name
degs_cd56$wikigene_description <- cd56_genes$wikigene_description
degs_cd56 <- degs_cd56[,c(1,20:21,2:19)]
degs_cd56 <- unique(degs_cd56)

## calculate max p value over all comparisons
degs_pval <- degs_cd56[,grep("_pval",colnames(degs_cd56))]
maxpval <- as.numeric(apply(degs_pval,1,max))
degs_cd56$max_pval <- maxpval

## sort by max p value
degs_cd56 <- degs_cd56[order(degs_cd56$max_pval),]

write.table(degs_cd56,"treatment_2_VS_control_CD56signaling_part2.csv",sep="\t",col.names=T,row.names=F,quote=F)


## hist of log2fc
pdf("log2_fc_mrn_hist_part2.pdf")
hist(degs_cd56$log2_fc_mrn,breaks=100)
dev.off()


## expression of CCL3, CCL4L2, CCL4 over time
candidates <- c("CCL3","CCL4L2","CCL4")
degs_candidates <- unlist(lapply(candidates,function(x){as.character(degs_cd56[degs_cd56$external_gene_name == x,]$id)}))


mrn <- read.csv("mrn.csv",sep="\t",header=T) 
mrn_candidates_pos <- unlist(lapply(degs_candidates,function(x){grep(x,rownames(mrn))}))
mrn_candidates <- mrn[mrn_candidates_pos,]

## CCL3
mrn_ccl3 <- data.frame(name=colnames(mrn_candidates),donor=c(rep("D4",3),rep("D6",3),rep("D2",2),rep("D5",2)), 
                       condition=c("NK-AF0","NK-AF10","NK-AF60","NK-AF0","NK-AF10","NK-AF60","NK-AF0","NK-AF60","NK-AF0","NK-AF60"), 
                       value=as.numeric(mrn_candidates[1,]))

pdf("CCL3.pdf")
ggplot(data=mrn_ccl3,aes(x=donor,y=value,fill=condition)) +
  geom_bar(stat="identity",position=position_dodge())
dev.off()

## CCL4L2
mrn_ccl4l2 <- data.frame(name=colnames(mrn_candidates),donor=c(rep("D4",3),rep("D6",3),rep("D2",2),rep("D5",2)), 
                       condition=c("NK-AF0","NK-AF10","NK-AF60","NK-AF0","NK-AF10","NK-AF60","NK-AF0","NK-AF60","NK-AF0","NK-AF60"), 
                       value=as.numeric(mrn_candidates[2,]))

pdf("CCL4L2.pdf")
ggplot(data=mrn_ccl4l2,aes(x=donor,y=value,fill=condition)) +
  geom_bar(stat="identity",position=position_dodge())
dev.off()

## CCL4
mrn_ccl4 <- data.frame(name=colnames(mrn_candidates),donor=c(rep("D4",3),rep("D6",3),rep("D2",2),rep("D5",2)), 
                         condition=c("NK-AF0","NK-AF10","NK-AF60","NK-AF0","NK-AF10","NK-AF60","NK-AF0","NK-AF60","NK-AF0","NK-AF60"), 
                         value=as.numeric(mrn_candidates[3,]))

pdf("CCL4.pdf")
ggplot(data=mrn_ccl4,aes(x=donor,y=value,fill=condition)) +
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


## pca only for cd65 genes ####
##
cd56_genes <- read.csv("treatment_2_VS_control_CD56signaling_part2.csv",sep="\t",header=T)
cd56_ids <- as.character(cd56_genes$id)

mrn <- read.csv("mrn.csv",sep="\t",header=T) 
mrn_genes <- rownames(mrn)
cd56_pos <- unlist(lapply(cd56_ids,function(x)grep(x,mrn_genes)))
mrn_cd56 <- mrn[cd56_pos,]
write.table(mrn_cd56,"mrn_CD56signaling_part2.csv",sep="\t",col.names=NA,row.names=T,quote=F)

counts <- read.csv("fc.counts.tsv",sep="\t",header=T)
counts_genes <- as.character(counts$X)
cd56_pos <- unlist(lapply(cd56_ids,function(x)grep(x,counts_genes)))
counts_cd56 <- counts[cd56_pos,]
write.table(counts_cd56,"counts_CD56signaling_part2.csv",sep="\t",col.names=T,row.names=F,quote=F)

