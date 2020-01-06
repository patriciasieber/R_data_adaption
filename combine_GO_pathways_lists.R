## combine lists of common GO terms and pathways between mouse and human

library("WriteXLS")


setwd("/data/Postdoc/Proteomics_SilkeMachata/results/")

lists_human <- list.files(pattern="regulated_significant.csv$", path="human/innatedb/",full.names=T)
lists_mouse <- list.files(pattern="regulated_significant.csv$", path="mouse/innatedb/",full.names=T)

comparisons <- c("go","pathways")
regulation <- c("up","down")

lapply(comparisons,function(comp){
  lists_comp <- lapply(list(lists_human,lists_mouse),function(current_list) current_list[grep(comp,current_list)])
  
  comp_regulation <- lapply(regulation,function(reg){
    lists_comp_reg <- lapply(lists_comp,function(current_list) current_list[grep(reg,current_list)])
    
    ## load data
    terms_comp_reg <- lapply(lists_comp_reg,function(current_file) read.csv(current_file,sep="\t",stringsAsFactors=F))
      
    ## overlap of mouse and human terms
    overlapping <- intersect(terms_comp_reg[[1]]$Pathway.Name, terms_comp_reg[[2]]$Pathway.Name)
    
    ## common lists of overlapping terms
    common_terms <- lapply(overlapping,function(ov){
      current_ov_human <- terms_comp_reg[[1]][terms_comp_reg[[1]]$Pathway.Name == ov,]
      current_ov_mouse <- terms_comp_reg[[2]][terms_comp_reg[[2]]$Pathway.Name == ov,]
      
      #if(current_ov_human$Pathway.Name != current_ov_mouse$Pathway.Id) print(paste0("different terms for ",ov))
      if(comp == "go"){
        current_ov <- current_ov_human[,c(1:3,5:6,8:9)]
        colnames(current_ov)[4:7] <- paste0("human_",colnames(current_ov)[4:7])
        current_ov <- cbind(current_ov,current_ov_mouse[,c(5:6,8:9)])
        colnames(current_ov)[8:11] <- paste0("mouse_",colnames(current_ov)[8:11])
      }else{
        current_ov <- current_ov_human[,c(1,3,2,5:6,8:9)]
        colnames(current_ov)[3:7] <- paste0("human_",colnames(current_ov)[3:7])
        current_ov <- cbind(current_ov,current_ov_mouse[,c(2,5:6,8:9)])
        colnames(current_ov)[8:12] <- paste0("mouse_",colnames(current_ov)[8:12])
      }
      return(current_ov)
    })
    common_terms <- do.call("rbind",common_terms)
    return(common_terms)
  })
  
  ## write to file
  WriteXLS::WriteXLS(
    comp_regulation,
    ExcelFileName = paste0("mail7/figures_final/common",comp,"terms.xlsx"),
    SheetNames = paste0(comp,"_",regulation,"regulation"),
    AdjWidth = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    row.names = F,
    col.names = T
  )
  
})


