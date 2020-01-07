## compare GO/pathway categories and underlying genes between enrichr and innatedb
## comparisons based on human data control vs treatment with Aspergillus fumigatus

library(VennDiagram)
library(readxl)

setwd("/data/enrichments/results/")

paths_innatedb <- list.files(pattern = "innatedb")
paths_enrichr <- list.files(pattern = "enrichr")


extract_sheets <- function(path,timepoints){
  current_sheets <- excel_sheets(path)
  sheets_tps <- lapply(timepoints,function(tp){
    tp_sheets <- current_sheets[grep(tp,current_sheets)]
    if(length(tp_sheets) == 1){
      current_enrichment <- read_excel(path, sheet=tp_sheets)
      return(current_enrichment)
    }else{
      current_enrichment <- lapply(tp_sheets,function(x){read_excel(path, sheet=x)})
      current_enrichment <- do.call("rbind",current_enrichment)
      return(current_enrichment)
    }
  })
  return(sheets_tps)
}


system("mkdir comparison_enrichrvsinnatedb")
sapply(1:4,function(i){
  print(i)
  ## specify sheets per timepoint
  tps <- c("early","middle","late")
  data_innatedb <- extract_sheets(paths_innatedb[i],tps)
  data_enrichr <- extract_sheets(paths_enrichr[i],tps)
  
  name_comparison <- strsplit(paste0(strsplit(paths_innatedb[i],"_")[[1]][c(1:2,4)],collapse="_"),"[.]")[[1]][1]
  if(strsplit(paste0(strsplit(paths_enrichr[i],"_")[[1]][c(1:2,4)],collapse="_"),"[.]")[[1]][1] != name_comparison) print("comparison wrong!")
  
  sapply(1:length(tps),function(j){
    print(j)
    current_innatedb <- data_innatedb[[j]]
    current_enrichr <- data_enrichr[[j]]
    
    if(nrow(current_innatedb) > 0 & nrow(current_enrichr) > 0){
      ## venn diagram enrichment terms
      if(grepl("go",name_comparison)){
        terms_innatedb <- paste0(current_innatedb$`Pathway Name`," (",current_innatedb$`Pathway Id`,")")
        terms_enrichr <- current_enrichr$Term
      }else{
        terms_innatedb <- current_innatedb$`Pathway Name`
        terms_enrichr <- sapply(current_enrichr$Term,function(x) strsplit(x,"_")[[1]][1])
      }

      
      common_terms <- intersect(terms_enrichr,terms_innatedb)
      #grid.newpage();
      pdf(paste0("comparison_enrichrvsinnatedb/",name_comparison,"_",tps[j],"_comparison_enrichmentterms_venn.pdf"))
      draw.pairwise.venn(area1=length(terms_enrichr), area2=length(terms_innatedb), 
                         cross.area=length(common_terms),
                         category=c("enrichr","innateDB"))
      dev.off()
      
      ## print overlap into file
      writeLines(common_terms,paste0("comparison_enrichrvsinnatedb/",name_comparison,"_",tps[j],"_common_enrichmentterms_list.txt"))
      
      
      ## venn diagram underlying genes
      genes_innatedb <- unlist(lapply(current_innatedb$`Gene Symbols`,function(x) strsplit(x,";")[[1]]))
      genes_innatedb <- unique(gsub(" ","",genes_innatedb))
      genes_enrichr <- unlist(lapply(current_enrichr$Genes,function(x) strsplit(x,";")[[1]]))
      genes_enrichr <- unique(genes_enrichr)
      
      common_genes <- intersect(genes_enrichr,genes_innatedb)
      #grid.newpage();
      pdf(paste0("comparison_enrichrvsinnatedb/",name_comparison,"_",tps[j],"_comparison_enrichmentgenes_venn.pdf"))
      draw.pairwise.venn(area1=length(genes_enrichr), area2=length(genes_innatedb), 
                         cross.area=length(common_genes),
                         category=c("enrichr","innateDB"))
      dev.off()
      
      ## print overlap into file
      writeLines(common_genes,paste0("comparison_enrichrvsinnatedb/",name_comparison,"_",tps[j],"_common_enrichmentgenes_list.txt"))
    }
  })
})








setwd("/data/results/")

paths_innatedb <- list.files(pattern = "innatedb")
paths_innatedb <- paths_innatedb[grep("xlsx",paths_innatedb)]
paths_enrichr <- list.files(pattern = "enrichr")
paths_enrichr <- paths_enrichr[grep("DC_VS_AF",paths_enrichr)]


extract_sheets_C_D3 <- function(path){
  current_sheets <- excel_sheets(path)
  if(length(current_sheets) == 1){
    current_enrichment <- read_excel(path, sheet=current_sheets)
  }else{
    current_enrichment <- lapply(current_sheets,function(x){read_excel(path, sheet=x)})
    current_enrichment <- do.call("rbind",current_enrichment)
  }
  return(current_enrichment)
}


system("mkdir comparison_enrichrvsinnatedb")
sapply(1:4,function(i){
  print(i)
  ## specify sheets per timepoint
  data_innatedb <- extract_sheets_C_D3(paths_innatedb[i])
  data_enrichr <- extract_sheets_C_D3(paths_enrichr[i])
  
  name_comparison <- strsplit(paste0(strsplit(paths_innatedb[i],"_")[[1]][c(1:2,4)],collapse="_"),"[.]")[[1]][1]
  if(strsplit(paste0(strsplit(paths_enrichr[i],"_")[[1]][c(1:2,7)],collapse="_"),"[.]")[[1]][1] != name_comparison) print("comparison wrong!")
  
    current_innatedb <- data_innatedb
    current_enrichr <- data_enrichr
    
    if(nrow(current_innatedb) > 0 & nrow(current_enrichr) > 0){
      ## venn diagram enrichment terms
      if(grepl("go",name_comparison)){
        terms_innatedb <- paste0(current_innatedb$`Pathway Name`," (",current_innatedb$`Pathway Id`,")")
        terms_enrichr <- current_enrichr$Term
      }else{
        terms_innatedb <- current_innatedb$`Pathway Name`
        terms_enrichr <- sapply(current_enrichr$Term,function(x) strsplit(x,"_")[[1]][1])
      }
      
      common_terms <- intersect(terms_enrichr,terms_innatedb)
      #grid.newpage();
      pdf(paste0("comparison_enrichrvsinnatedb/",name_comparison,"_comparison_enrichmentterms_venn.pdf"))
      draw.pairwise.venn(area1=length(terms_enrichr), area2=length(terms_innatedb), 
                         cross.area=length(common_terms),
                         category=c("enrichr","innateDB"))
      dev.off()
      
      ## print overlap into file
      writeLines(common_terms,paste0("comparison_enrichrvsinnatedb/",name_comparison,"_common_enrichmentterms_list.txt"))
      
      
      ## venn diagram underlying genes
      genes_innatedb <- unlist(lapply(current_innatedb$`Gene Symbols`,function(x) strsplit(x,";")[[1]]))
      genes_innatedb <- unique(gsub(" ","",genes_innatedb))
      genes_enrichr <- unlist(lapply(current_enrichr$Genes,function(x) strsplit(x,";")[[1]]))
      genes_enrichr <- unique(genes_enrichr)
      
      common_genes <- intersect(genes_enrichr,genes_innatedb)
      #grid.newpage();
      pdf(paste0("comparison_enrichrvsinnatedb/",name_comparison,"_comparison_enrichmentgenes_venn.pdf"))
      draw.pairwise.venn(area1=length(genes_enrichr), area2=length(genes_innatedb), 
                         cross.area=length(common_genes),
                         category=c("enrichr","innateDB"))
      dev.off()
      
      ## print overlap into file
      writeLines(common_genes,paste0("comparison_enrichrvsinnatedb/",name_comparison,"_common_enrichmentgenes_list.txt"))
    }
})
