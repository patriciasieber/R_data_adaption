
## part1: images of GO hierarchy using amigo ####
## filter GO terms by CC, BP, MF - all up- and down-regulated ones

## read files
path_to_human <- "/data/Postdoc/Proteomics_SilkeMachata/results/human/innatedb/"
path_to_mouse <- "/data/Postdoc/Proteomics_SilkeMachata/results/mouse/innatedb/"

## filter for main GO categories
categories <- c("cellular component","biological process","molecular function")
for(g in c(path_to_human,path_to_mouse)){
  go_files <- list.files(g,"_go_",full.names=T)
  go_files <- go_files[grep("regulated_significant",go_files)]
  
  g_lists <- lapply(go_files,function(x) read.csv(x,sep="\t"))
  g_list_all <- do.call("rbind",g_lists)
  
  go_per_category <- lapply(categories,function(c){
    current_terms <- unique(as.character(g_list_all[g_list_all$Source.Name == c,]$Pathway.Id))
    c <- sub(" ","",c)
    ## into file 
    write(current_terms,paste0(g,"all_go_",c,".txt"))
  })
}


## -> load into http://amigo.geneontology.org/visualize?mode=client_amigo
## adapt to see up and down regulated terms

## part2: plots of GO term abundance ####
path_to_human <- "/data/Postdoc/Proteomics_SilkeMachata/results/human/innatedb/revigo/"
path_to_mouse <- "/data/Postdoc/Proteomics_SilkeMachata/results/mouse/innatedb/revigo/"


regulation <- c("up","down")
for(g in c(path_to_human,path_to_mouse)){
  ## read revigo file
  revigo_files <- list.files(g,".csv$",full.names=T)
  go_files <- list.files(paste0(g,"../"),"_significant.csv$",full.names=T)
  
  for(reg in regulation){
    print(reg)
    revigo_files_reg <- revigo_files[grep(reg,revigo_files)]
    go_file_reg <- go_files[grep(paste0("go_",reg,"regulated"),go_files)]
    current_go_reg <- read.csv(go_file_reg,sep="\t")
    
    
    ## identify go terms of first (and second) level
    for(rev in revigo_files_reg){
      print(rev)
      current_revigo <- read.csv(rev)
      level1 <- current_revigo[current_revigo$eliminated == 0,]
      
      #go_level1_boolean <- unlist(lapply(as.character(current_go_reg$Pathway.Id),function(x){x %in% as.character(level1$term_ID)}))
      go_level1_index <- unlist(lapply(as.character(level1$term_ID),function(x){which(as.character(current_go_reg$Pathway.Id) == x)}))
      go_level1 <- current_go_reg[go_level1_index,]        

      go_level12_index <- unlist(lapply(as.character(current_revigo$term_ID),function(x){which(as.character(current_go_reg$Pathway.Id) == x)}))
      go_level12 <- current_go_reg[go_level12_index,]
      
      ## --> take size, term count, corr p value
      
      ## plot1: only count of GO terms
      new_xaxes <- seq(0,0.05,0.01)
      
      if(nrow(go_level1) > 0){
        #svg(paste0(rev,"_GOcount_level1.svg"),width=20,height=15)
        pdf(paste0(rev,"_GOcount_level1.pdf"))
        barplot(rev(go_level1$Pathway.uploaded.gene.count),horiz=TRUE,
                #names.arg=rev(as.character(go_level1$Pathway.Id)),
                names.arg=rev(as.character(go_level1$Pathway.Name)),
                las=2,cex.names=0.5,
                col="blue")
        dev.off()
        
        
        ## plot2: p-value, gene counts in header        
        ## sort results by pvalue
        go_level1 <- go_level1[order(go_level1$Pathway.p.value..corrected.,decreasing=F),]
        
        #svg(paste0(rev,"_GOpvalue_level1.svg"),width=20,height=15)
        pdf(paste0(rev,"_GOpvalue_level1.pdf"))
        barplot(rev(0.05-go_level1$Pathway.p.value..corrected.),horiz=TRUE,
                names.arg=rev(paste(as.character(go_level1$Pathway.Name),go_level1$Pathway.uploaded.gene.count,go_level1$Genes.in.InnateDB.for.this.entity,sep="  ")),
                #names.arg=rev(paste(as.character(go_level1$Pathway.Id),go_level1$Pathway.uploaded.gene.count,go_level1$Genes.in.InnateDB.for.this.entity,sep="  ")),
                las=2,cex.names=0.5,col="blue",
                xaxt="n")
        axis(1,at=new_xaxes,labels=rev(new_xaxes))
        dev.off()
        
        
        
        ## remove large terms (>= 1000 genes too broad)
        go_level1 <- go_level1[go_level1$Genes.in.InnateDB.for.this.entity < 1000,]
        #go_level12 <- go_level12[go_level12$Genes.in.InnateDB.for.this.entity < 1000,]
        
        #svg(paste0(rev,"_GOcount_level1_specific.svg"),width=20,height=15)
        pdf(paste0(rev,"_GOcount_level1_specific.pdf"))
        barplot(rev(go_level1$Pathway.uploaded.gene.count),horiz=TRUE,
                #names.arg=rev(as.character(go_level1$Pathway.Id)),las=2,cex.names=0.5,
                names.arg=rev(as.character(go_level1$Pathway.Name)),las=2,cex.names=0.5,
                col="blue")
        dev.off()
        
        
        ## plot2: p-value, gene counts in header
        #svg(paste0(rev,"_GOpvalue_level1_specific.svg"),width=20,height=15)
        pdf(paste0(rev,"_GOpvalue_level1_specific.pdf")) #,height=15
        barplot(rev(0.05-go_level1$Pathway.p.value..corrected.),horiz=TRUE,
                #names.arg=rev(paste(as.character(go_level1$Pathway.Id),go_level1$Pathway.uploaded.gene.count,go_level1$Genes.in.InnateDB.for.this.entity,sep="  ")),
                names.arg=rev(paste(as.character(go_level1$Pathway.Name),go_level1$Pathway.uploaded.gene.count,go_level1$Genes.in.InnateDB.for.this.entity,sep="  ")),
                las=2,cex.names=0.5,col="blue",
                xaxt="n")
        axis(1,at=new_xaxes,labels=rev(new_xaxes))
        dev.off()
        
        
        ## print table for go_level1 (Zahl der Go-terms; Genzahlen; korr. pvalues); order first
        print_go_level1 <- go_level1[,c(1:3,5:6,8:9)]
        print_go_level1 <- print_go_level1[order(print_go_level1$Pathway.p.value..corrected.,decreasing=F),]
        
        WriteXLS::WriteXLS(
          print_go_level1,
          ExcelFileName = paste0(rev,"_go_level1_specific.xls"),
          SheetNames = NULL,
          AdjWidth = T,
          BoldHeaderRow = T,
          FreezeRow = 1,
          FreezeCol = 1,
          row.names = F,
          col.names = T
        )
      }

      
      ## include all levels
      if(nrow(go_level12) > 0){
        go_level2_index <- which(current_revigo$eliminated == 1)
        col_level12 <- rep("blue",length(go_level12_index))
        col_level12[go_level2_index] <- "lightblue"
        
        ## plot1
        #svg(paste0(rev,"_GOcount_alllevels.svg"),width=20,height=15)
        pdf(paste0(rev,"_GOcount_alllevels.pdf"))
        barplot(rev(go_level12$Pathway.uploaded.gene.count),horiz=TRUE,
                #names.arg=rev(as.character(go_level12$Pathway.Id)),
                names.arg=rev(as.character(go_level12$Pathway.Name)),
                las=2,cex.names=0.5,
                col=rev(col_level12))
        dev.off()
        
        ## plot2
        go_level12 <- go_level12[order(go_level12$Pathway.p.value..corrected.,decreasing=F),]
        #svg(paste0(rev,"_GOpvalue_alllevels.svg"),width=20,height=15)
        pdf(paste0(rev,"_GOpvalue_alllevels.pdf")) 
        barplot(rev(go_level12$Pathway.p.value..corrected.),horiz=TRUE,
                #names.arg=rev(paste(as.character(go_level12$Pathway.Id),go_level12$Pathway.uploaded.gene.count,go_level12$Genes.in.InnateDB.for.this.entity,sep="  ")),
                names.arg=rev(paste(as.character(go_level12$Pathway.Name),go_level12$Pathway.uploaded.gene.count,go_level12$Genes.in.InnateDB.for.this.entity,sep="  ")),
                las=2,cex.names=0.5,col=rev(col_level12),
                xaxt="n")
        axis(1,at=new_xaxes,labels=rev(new_xaxes))
        dev.off()
        
      }

    }
    
  }
  
}





## top 10 for mouse and human, up- and down-regulated
path_to_human <- "/data/Postdoc/Proteomics_SilkeMachata/results/human/innatedb/revigo/"
path_to_mouse <- "/data/Postdoc/Proteomics_SilkeMachata/results/mouse/innatedb/revigo/"

regulation <- c("up","down")
for(g in c(path_to_human,path_to_mouse)){
  ## read revigo file
  revigo_files <- list.files(g,".csv$",full.names=T)
  go_files <- list.files(paste0(g,"../"),"_significant.csv$",full.names=T)
  
  lapply(regulation,function(reg){
    revigo_files_reg <- revigo_files[grep(reg,revigo_files)]
    go_file_reg <- go_files[grep(paste0("go_",reg,"regulated"),go_files)]
    current_go_reg <- read.csv(go_file_reg,sep="\t")
    
    ## identify go terms of first level, extract top 10
    go_level1_categories <- lapply(revigo_files_reg,function(rev){
      current_revigo <- read.csv(rev)
      level1 <- current_revigo[current_revigo$eliminated == 0,]
      go_level1_index <- unlist(lapply(as.character(level1$term_ID),function(x){which(as.character(current_go_reg$Pathway.Id) == x)}))
      go_level1 <- current_go_reg[go_level1_index,]      
      
      ## remove large terms (>= 500 genes too broad)
      go_level1 <- go_level1[go_level1$Genes.in.InnateDB.for.this.entity < 500,]   ## only this: _GOpvalue_top10_level1_specific.pdf
      go_level1 <- go_level1[go_level1$Pathway.uploaded.gene.count > 1,]            ## also this: _GOpvalue_top10_level1_wo1_specific.pdf
      ## sort by corrected p value
      go_level1 <- go_level1[order(go_level1$Pathway.p.value..corrected.),]
      return(go_level1[1:10,])
    })
    go_level1_categories <- do.call("rbind",go_level1_categories)
    
    
    new_xaxes <- seq(0,0.05,0.01)
    ## BP, CC, MF
    color <- c("#009900","#0099FF","#FF6600")
    
    pdf(paste0(g,reg,"_GOpvalue_top10_level1_wo1-500_specific.pdf")) #,height=15
    barplot(rev(0.05-go_level1_categories$Pathway.p.value..corrected.),horiz=TRUE,
            names.arg=rev(as.character(go_level1_categories$Pathway.Name)),
            las=2,cex.names=0.7,col=rev(rep(color,each=10)),border=NA,
            xaxt="n",xlim=c(0,0.05),
            legend.text=rev(as.character(unique(go_level1_categories$Source.Name))),
            args.legend=list(x="bottomright",fill=color,border=NA,bty = "n",cex = 0.75))
    axis(1,at=new_xaxes,labels=rev(new_xaxes),cex.axis=0.8)
    title(xlab="corrected p-value",cex.lab=0.9)
    dev.off()
    
    
    ## again with log scale
    log10_pvalue <- log10(go_level1_categories$Pathway.p.value..corrected.)
    new_xaxes <- seq(0,ceiling(abs(min(log10_pvalue)))+5,5)
    if(length(new_xaxes) <=2) new_xaxes <- seq(0,ceiling(abs(min(log10_pvalue)))+1,1)
    ## adapt to 1E+00, 1E-01, 1E-02, 1E-03 ..
    new_xaxes_labels <- sapply(1:length(new_xaxes),function(x){
      current_value <- new_xaxes[x]
      if(current_value == 0){
        return("1E+00") ## first value is zero
      }else{
        if(nchar(current_value) == 1){
          current_value <- paste0("0",current_value)
        }
        current_value <- paste0("1E-",current_value)
        return(current_value)
      }
    })
    
    pdf(paste0(g,reg,"_GOpvalue_top10_level1_wo1-500_log10_specific.pdf")) #,height=15
    barplot(rev(abs(log10_pvalue)),horiz=TRUE,
            names.arg=rev(as.character(go_level1_categories$Pathway.Name)),
            las=2,cex.names=0.7,col=rev(rep(color,each=10)),border=NA,
            xaxt="n",#xlim=c(0,0.05),
            legend.text=rev(as.character(unique(go_level1_categories$Source.Name))),
            args.legend=list(x="bottomright",fill=color,border=NA,bty = "n",cex = 0.75))
    axis(1,at=new_xaxes,labels=new_xaxes_labels,cex.axis=0.8)
    title(xlab="log10 (corrected p-value)",cex.lab=0.9)
    dev.off()
    
    ## identify go terms, extract top 10
    categories <- as.character(unique(current_go_reg$Source.Name))
    go_categories <- lapply(categories,function(cat){
      go_cat <- current_go_reg[current_go_reg$Source.Name == cat,]
      ## remove large terms (>= 1000 genes too broad)
      go_cat <- go_cat[go_cat$Genes.in.InnateDB.for.this.entity < 1000,]
      ## sort by corrected p value
      go_cat <- go_cat[order(go_cat$Pathway.p.value..corrected.),]
      return(go_cat[1:10,])
    })
    go_categories <- do.call("rbind",go_categories)
    
    
    new_xaxes <- seq(0,0.05,0.01)
    ## BP, CC, MF
    color <- c("#009900","#0099FF","#FF6600")
    
    pdf(paste0(g,reg,"_GOpvalue_top10_specific.pdf")) #,height=15
    barplot(rev(0.05-go_categories$Pathway.p.value..corrected.),horiz=TRUE,
            names.arg=rev(as.character(go_categories$Pathway.Name)),
            las=2,cex.names=0.7,col=rev(rep(color,each=10)),border=NA,
            xaxt="n",xlim=c(0,0.05),
            legend.text=rev(as.character(unique(go_categories$Source.Name))),
            args.legend=list(x="bottomright",fill=color,border=NA,bty = "n",cex = 0.75))
    axis(1,at=new_xaxes,labels=rev(new_xaxes),cex.axis=0.8)
    title(xlab="corrected p-value",cex.lab=0.9)
    dev.off()
    
  })
  
  
  
}  
  
