## calculate correlation 
## between genes within a condition (control vs S. aureus/ N. meningitidis/ A. fumigatus)
## Author: Patricia Sieber
## Date: 23.04.2019, edit: 14.05.2019

library(pheatmap)
library(RColorBrewer)


setwd("/home/psieber/correlation_X_B_genelist/")
R.utils::sourceDirectory("/sbidata/pidomics/X_B_eQTL/src/R/")

#load("R_workspace.RData")
#save(counts_bonn,get_mrn,file="correlation/correlation.RData")
#rm(list=ls())

load("correlation.RData")


# work with count data from Bonn, time point 6 hours
# each infection, including control, separatly
select <- c("K_6h","Af_6h","Sa_6h","Nm_6h")

samples <- lapply(select,function(x){grep(x, colnames(counts_bonn))})
counts_all_genes_samples_mrn <- lapply(samples,function(x){get_mrn(counts_bonn[,x])})

gene_list <- read.table("gene_list.txt")$V1

# counts_all_genes_control_samples_mrn = get_mrn(counts_bonn[, samples]) # apply MRN to all genes and only K_6h samples!
# # counts_all_genes_control_samples_mrn_scale = t(scale(t(counts_all_genes_control_samples_mrn))) # many NaNs because many 0s in mrn --> 0 var --> division by 0
# counts_all_genes_control_samples_mrn_plus1_log = log2(counts_all_genes_control_samples_mrn+1) # scale() not applicable here because of ^^^

## get deg lists
deg_dir <- "/sbidata/pidomics/X_B_eQTL/res/diff_exp_genes_bonn/"
#deg_files_6h <- c("Afum_6_VS_control_6.csv","Saur_6_VS_control_6.csv","Nmen_6_VS_control_6.csv",
#                  "Afum_6_VS_Saur_6.csv","Afum_6_VS_Nmen_6.csv","Nmen_6_VS_Saur_6.csv")
deg_files_6h <- c("Afum_6_VS_control_6.csv","Saur_6_VS_control_6.csv","Nmen_6_VS_control_6.csv")
degs_lists_6h <- lapply(deg_files_6h,function(x){
  temp_csv <- read.csv(paste0(deg_dir,x),sep="\t")
  temp_csv <- temp_csv[temp_csv$DESeq == T & temp_csv$DESeq2 == T & temp_csv$Limma == T & temp_csv$EdgeR == T,]
  return(as.character(temp_csv$id))
})

degs_K_6 <- unique(unlist(degs_lists_6h[grepl("control",deg_files_6h)]))
degs_Af_6 <- unique(unlist(degs_lists_6h[grepl("Afum",deg_files_6h)]))                    
degs_Sa_6 <- unique(unlist(degs_lists_6h[grepl("Saur",deg_files_6h)]))  
degs_Nm_6 <- unique(unlist(degs_lists_6h[grepl("Nmen",deg_files_6h)]))  

degs_concat <- c("degs_K_6","degs_Af_6","degs_Sa_6","degs_Nm_6")

## pairwise calculation of correlation and corresponding pvalue between degs

#ptm <- proc.time()
cor_lists <- NULL
for(s in 1:length(select)){
  print(select[s])
  current_selection <- counts_all_genes_samples_mrn[[s]]
  current_degs <- get(degs_concat[s])
  cor.list <- data.frame("gene1"=NULL,"gene2"=NULL,"correlation"=NULL,"pvalue"=NULL)
  
  len_genes <- length(current_degs)
  for(g1 in 1:length(gene_list)){
    for(g2 in 1:len_genes){  
      gene1 <- as.character(gene_list[g1])
      gene2 <- as.character(current_degs[g2])
      gene1_mrn <- as.numeric(current_selection[rownames(current_selection) == gene1,])
      gene2_mrn <- as.numeric(current_selection[rownames(current_selection) == gene2,])
      
      test.cor <- cor.test(gene1_mrn,gene2_mrn,method="pearson")
      current_correlation <- as.numeric(test.cor$estimate)
      current_pvalue <- test.cor$p.value
      
      sub_list <- data.frame("gene1"=as.character(gene1),"gene2"=as.character(gene2),
                             "correlation"=as.numeric(current_correlation),"pvalue"=as.numeric(current_pvalue))
      cor.list <- rbind(cor.list,sub_list)   
    }
  }
  cor_lists[[s]] <- cor.list

}
names(cor_lists) <- select
#proc.time() - ptm
gc()


## TODO: adjust pvalue
cor_lists <- lapply(cor_lists,function(x){
  x$fdr <- p.adjust(x$pvalue, method="fdr")
  return(x)
})



## print lists
lapply(1:length(select),function(s){
  current_cor_list <- cor_lists[[s]]
  current_name <- select[s]
  write.table(current_cor_list,paste0("correlation_",current_name,".csv"),sep="\t",col.names=T,row.names=F)
})



## select by p value thresholds
## p < 0.05
cor_lists_sign <- lapply(cor_lists,function(x){
  x <- x[!is.na(x$fdr),]
  return(x[x$fdr < 0.05,])
})

## save lists to file
lapply(1:length(select),function(s){
  current_cor_list <- cor_lists_sign[[s]]
  current_name <- select[s]
  write.table(current_cor_list,paste0("correlation_",current_name,"_significant.csv"),sep="\t",col.names=T,row.names=F)
})

## correlation distribution plot
for(i in 1:length(select)){
  pdf(paste0(select[i],"_correlation_sigificant_hist.pdf"))
  hist(as.numeric(cor_lists_sign[[i]]$correlation),xlab="Correlation",main=paste0("Pairwise correlation of genes for ",select[i]))
  dev.off()
}


cor_lists_sign_cor05 <- lapply(cor_lists_sign,function(x){
  #x <- x[!is.na(x$fdr),]
  x <- x[x$correlation < -0.5 | x$correlation > 0.5,]
  print(nrow(x))
  return(x)
})


## check which difference between treatment vs condition!!
control <- cor_lists_sign_cor05[[1]]
control_name <- select[1]
corr_diff_vscontrol <- NULL
for(i in 2:length(select)){
  treatment <- cor_lists_sign_cor05[[i]]
  treatment_name <- select[i]
  
  corr_comp <- NULL
  for(gene1 in as.character(gene_list)){
    print(gene1)
    treatment_g1 <- treatment[treatment$gene1 == gene1,]
    control_g1 <- control[control$gene1 == gene1,]
    
    compared_treat_boolean <- vector(length=nrow(treatment_g1))
    compared_contr_boolean <- vector(length=nrow(control_g1))
    
    ## compare treatment vs control
    if(nrow(treatment_g1) > 0){
      for(j in 1:nrow(treatment_g1)){
        treatment_g1g2 <- treatment_g1[j,]
        compared_treat_boolean[j] <- T
        
        status <- ""
        pos_contr_g1g2 <- which(as.character(control_g1$gene2) == as.character(treatment_g1g2$gene2))
        if(length(pos_contr_g1g2) == 1){
          control_g1g2 <- control_g1[pos_contr_g1g2,]
          compared_contr_boolean[pos_contr_g1g2] <- T  ## already compared, does not need to be compared again in controls
          
          ## unchanged correlation direction, pos vs neg. corr, neg vs pos. corr
          if((treatment_g1g2$correlation > 0.5 & control_g1g2$correlation > 0.5) | (treatment_g1g2$correlation < -0.5 & control_g1g2$correlation < -0.5) ){
            status <- "unchanged"
          }else{
            if(treatment_g1g2$correlation > 0.5 & control_g1g2$correlation < -0.5){
              status <- "pos.corrvsneg.corr"
            }else{
              if(treatment_g1g2$correlation < -0.5 & control_g1g2$correlation > 0.5){
                status <- "neg.corrvspos.corr"
              }else{
                status <- "undefined"
                print(j)
              }
            }
          }
          new_entry <- data.frame(gene1=as.character(treatment_g1g2$gene1), gene2=as.character(treatment_g1g2$gene2), 
                                  status_treatmentvscontrol=status, 
                                  correlation_tr=treatment_g1g2$correlation, fdr_tr=treatment_g1g2$fdr, 
                                  correlation_c=control_g1g2$correlation, fdr_c=control_g1g2$fdr,stringsAsFactors=FALSE)
        }else{
          ##status sigvsnon-sig: not significant in control
          status <- "sigvsnon-sig"
          new_entry <- data.frame(gene1=as.character(treatment_g1g2$gene1), gene2=as.character(treatment_g1g2$gene2), 
                                  status_treatmentvscontrol=status, 
                                  correlation_tr=treatment_g1g2$correlation, fdr_tr=treatment_g1g2$fdr, 
                                  correlation_c="", fdr_c="",stringsAsFactors=FALSE)
        } 
        corr_comp <- rbind(corr_comp,new_entry)
      }
      
    }

    
    ##status non-sigvssig --> in control but not in treatment
    if(nrow(control_g1) > 0){
      open_comparisons <- which(!compared_contr_boolean)
      for(j in open_comparisons){
        control_g1g2 <- control_g1[j,]
        
        ##status non-sigvssig: not significant in treatment
        status <- "non-sigvssig"
        new_entry <- data.frame(gene1=as.character(control_g1g2$gene1), gene2=as.character(control_g1g2$gene2), 
                                status_treatmentvscontrol=status, 
                                correlation_tr="", fdr_tr="", 
                                correlation_c=control_g1g2$correlation, fdr_c=control_g1g2$fdr,stringsAsFactors=FALSE)
        corr_comp <- rbind(corr_comp,new_entry)
      }
    }
  }
  #corr_comp
  ## print results
  write.table(corr_comp,paste0("correlation_comparison_",treatment_name,"vs",control_name,".csv"),sep="\t",col.names=T,row.names=F)
  
  changed_corr <- corr_comp[corr_comp$status_treatmentvscontrol != "unchanged",]
  write.table(changed_corr,paste0("changed_correlation_comparison_",treatment_name,"vs",control_name,".csv"),sep="\t",col.names=T,row.names=F)
  
  corr_diff_vscontrol[[i]] <- corr_comp
}

save.image("correlation_pearson.RData")
  

## read biomart to change gene ids into readable gene names
biomart <- read.csv("mart_export_id_map_20191001.txt",header=T)


gene_list_mart_index <- unlist(lapply(as.character(gene_list),function(x){which(as.character(biomart$Gene.stable.ID) == x) }))
gene_list_mart <- biomart[gene_list_mart_index,]

corr_diff_vscontrol_names <- NULL
for(i in 2:length(select)){
  current_comp <- corr_diff_vscontrol[[i]]
  treatment_name <- select[i]
  
  current_comp <- current_comp[current_comp$status_treatmentvscontrol != "unchanged",]
  
  current_comp$gene1 <- unlist(lapply(current_comp$gene1,function(x){
    output <- as.character(gene_list_mart[gene_list_mart$Gene.stable.ID == x,]$Gene.name)
    if(length(output) == 0) output <- x
    return(output[1])
  }))
  
  current_comp$gene2 <- unlist(lapply(current_comp$gene2,function(x){
    output <- as.character(biomart[biomart$Gene.stable.ID == x,]$Gene.name)
    if(length(output) == 0) output <- x
    return(output[1])
  }))
  
  write.table(current_comp,paste0("changed_correlation_comparison_",treatment_name,"vs",control_name,"_genenames.csv"),sep="\t",col.names=T,row.names=F)
  corr_diff_vscontrol_names[[i]] <- current_comp
}


## output for GO analyses
for(i in 2:length(select)){
  current_comp <- corr_diff_vscontrol_names[[i]]
  treatment_name <- select[i]
  
  genes_tr <- current_comp[current_comp$correlation_tr != "" & current_comp$correlation_c == "",] 
  genes_c <- current_comp[current_comp$correlation_c != "" & current_comp$correlation_tr == "",]
  
  genes_tr <- unique(c(genes_tr$gene1,genes_tr$gene2))
  genes_c <- unique(c(genes_c$gene1,genes_c$gene2))
  
  writeLines(genes_tr,paste0("genenames_GO_",treatment_name,"vs",control_name,"_treatment.txt"))
  writeLines(genes_c,paste0("genenames_GO_",treatment_name,"vs",control_name,"_control.txt"))
}


for(i in 2:length(select)){
  treatment_name <- select[i]
  current_counts <- unique(c(rownames(counts_all_genes_samples_mrn[[1]]),rownames(counts_all_genes_samples_mrn[[i]])))
  
  current_genenames <- unlist(lapply(current_counts,function(x){
    output <- as.character(biomart[biomart$Gene.stable.ID == x,]$Gene.name)
    if(length(output) == 0) output <- x
    return(output[1])
  }))
  writeLines(current_genenames,paste0("genenames_GO_",treatment_name,"vs",control_name,"_background.txt"))
}

