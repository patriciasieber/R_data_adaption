## enrichment analysis for human data using enrichR and globaltest
## different states of control vs infection (early, middle, late) are compared

# install.packages("enrichR")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("globaltest")
# install.packages("readxl")

library(enrichR)        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987924/
library(globaltest)     ## https://www.ncbi.nlm.nih.gov/pubmed/14693814
library("readxl")
library("WriteXLS")
library(Biobase)
library(annotate)

## read data
dem_degs <- read_excel("data/noP10--deg-dem-overlap.xlsx")

dem_alltargets <- read_excel("data/noP10--deg-all-targets-overlap.xlsx")


# overlap of dem_alltargets and all degs!  -- same results as for dem_degs
timepoints <- c("early","middle","late")
# degs <- c("/data/enrichments/data/noP10.T12.Wald.status.xls",
#           "/data/enrichments/data/noP10.T34.Wald.status.xls",
#           "/data/enrichments/data/noP10.T5.Wald.status.xls")
# demdeg_overlap <- lapply(seq(1,length(timepoints)),function(i){
#   x <- degs[i]
#   current <- read_excel(x)
#   current <- current[current$sig == T,]
#   colnames(current)[1] <- "mRNA"
#   
#   current_alltargets <- dem_alltargets[dem_alltargets$.id == timepoints[i],]
#   
#   overlap_boolean <- unlist(lapply(current_alltargets$mRNA,function(y){y %in% current$mRNA}))
#   overlap <- current_alltargets[overlap_boolean,]
#   return(overlap)
# })
# demdeg_overlap <- do.call("rbind",demdeg_overlap)
 


## enrichment based on enrichr ####
## relevant databases:  http://amp.pharm.mssm.edu/Enrichr/
dbs_go <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
dbs_pathways <- c("Reactome_2016","KEGG_2019_Human")

apply_enrichr <- function(gene_object, database, timepoints=c("early","middle","late"), output_file){
  ## run enrichr for relevant databases, separately for each timepoint
  enrichment <- lapply(timepoints,function(x){
    return(enrichr(gene_object[gene_object$.id == x,]$mRNA, databases = database))
  })
  
  ## filter for significant adjusted p value
  sign_enrichment <- lapply(enrichment,function(x){
    sign <- lapply(x,function(y){
      return(y[y$Adjusted.P.value < 0.05,])
    })
    return(sign)
  })
  
  ## print lists to file
  name_enrichments <- lapply(sign_enrichment,function(x){return(names(x))})
  name_tp <- lapply(seq(1,length(name_enrichments)),function(i){
    rep(timepoints[i],length(name_enrichments[[i]]))
  })
  output_names <- paste0(unlist(name_tp),"_",unlist(name_enrichments))
  output_names <- unlist(lapply(output_names,function(x){
    substr(x,1,30) ## limit for names of sheets in xls
  }))
  res <- do.call(c, sign_enrichment)
  
  WriteXLS::WriteXLS(
    res,
    ExcelFileName = output_file,
    SheetNames = output_names,
    AdjWidth = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    row.names = F,
    col.names = T
  )
  return(sign_enrichment)
}

  
dem_degs_enrichr_go <- apply_enrichr(dem_degs,dbs_go,output_file="results/dem_degs_enrichr_go.xls")
dem_degs_enrichr_pathways <- apply_enrichr(dem_degs,dbs_pathways,output_file="results/dem_degs_enrichr_pathways.xls")

dem_alltargets_enrichr_go <- apply_enrichr(dem_alltargets,dbs_go,output_file="results/dem_alltargets_enrichr_go.xls")
dem_alltargets_enrichr_pathways <- apply_enrichr(dem_alltargets,dbs_pathways,output_file="results/dem_alltargets_enrichr_pathways.xls")



## enrichment based on globaltest ####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GO.db")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("KEGG.db")


expr <- read.csv("/data/mrn.csv",sep="\t")
colnames(expr) <- gsub("[.]","-",colnames(expr))


## only consider samples of interest (early, middle, late)
## according to /sbidata/pidomics/A_C3_Biobank_mRNA/doc/meta_data_timed_interval.v2.unique.xlsx

samplesofinterest <- NULL # TODO: fill!
nrcols <- unlist(lapply(samplesofinterest,function(x){grep(x,colnames(expr))}))
expr <- expr[,nrcols]


## only dem_degs and dem_alltargets
dem_degs_expr_boolean <- sapply(rownames(expr),function(x){x %in% dem_degs$mRNA})
dem_degs_expr <- expr[dem_degs_expr_boolean,]

dem_alltargets_expr_boolean <- sapply(rownames(expr),function(x){x %in% dem_alltargets$mRNA})
dem_alltargets_expr <- expr[dem_alltargets_expr_boolean,]


## devide into early, middle, late
dem_degs_expr_tps <- lapply(timepoints,function(x){
  x_genes <- dem_degs[dem_degs$.id == x,]$mRNA
  dem_degs_expr_boolean <- sapply(rownames(dem_degs_expr),function(y){y %in% x_genes})
  dem_degs_expr_new <- dem_degs_expr[dem_degs_expr_boolean,]
  return(dem_degs_expr_new)
})

dem_alltargets_expr_tps <- lapply(timepoints,function(x){
  x_genes <- dem_alltargets[dem_alltargets$.id == x,]$mRNA
  dem_alltargets_expr_boolean <- sapply(rownames(dem_alltargets_expr),function(y){y %in% x_genes})
  dem_alltargets_expr_new <- dem_alltargets_expr[dem_alltargets_expr_boolean,]
  return(dem_alltargets_expr_new)
})




## design matrix for globaltest of remaining samples
patient <- unlist(lapply(colnames(expr),function(x) strsplit(x,"-")[[1]][1]))
design <- data.frame(patient=patient,
                    timepoint=rep(c("early","middle","late"),6),
                    condition=rep(c(rep("IA",3),rep("control",3)),3),
                    sex=c(rep("f",3*2),rep("m",3*4)))
rownames(design) <- colnames(expr)


eg <- as.list(org.Hs.egALIAS2EG)  ## to convert gene names into entrez ids, which are used internally by globaltest; adapt when different identifiers are used
# 
# lapply(seq(1,length(timepoints)),function(i){
#   current_dem_degs_exprset <- ExpressionSet(as.matrix(dem_degs_expr_tps[[i]]),phenoData=AnnotatedDataFrame(design))
#   #test <- gt(condition,current_dem_degs_exprset)
#   
#   dem_degs_globaltest_go <- gtGO(condition, current_dem_degs_exprset, annotation="org.Hs.eg", probe2entrez=eg, multtest="BH")
#   dem_degs_globaltest_pathways <- gtKEGG(condition, current_dem_degs_exprset, annotation="org.Hs.eg", probe2entrez=eg, multtest="BH") 
#   
#   ## läuft!
# })
 

apply_globaltest <- function(gene_object, timepoints=c("early","middle","late"), enrichment_type,
                             designmatrix, annot="org.Hs.eg", id2entrez=eg, test="BH", output_file){
  
  ## run globaltest for relevant databases, separately for each timepoint
  enrichment <- lapply(seq(1,length(timepoints)),function(i){
    current_exprset <- ExpressionSet(as.matrix(gene_object[[i]]),phenoData=AnnotatedDataFrame(designmatrix))
    
    if(enrichment_type == "GO") current_globaltest <- gtGO(condition, current_exprset, annotation=annot, probe2entrez=id2entrez, multtest=test)
    if(enrichment_type == "KEGG") current_globaltest <- gtKEGG(condition, current_exprset, annotation=annot, probe2entrez=id2entrez, multtest=test) 
    return(current_globaltest)
  })
  
  ## filter for significant adjusted p value
  sign_enrichment <- lapply(enrichment,function(x){
    #table of x@result, x@extra, x@subsets
    x_new <- cbind(x@result,x@extra)
    x_new$genes <- sapply(x@subsets,function(y){paste(y,collapse=";")})
    x_new <- x_new[x_new$BH < 0.05,]
    return(x_new)
  })
  
  
  ## print lists to file
  if(enrichment_type == "GO") output_names <- paste0(timepoints,"_GO")
  if(enrichment_type == "KEGG") output_names <- paste0(timepoints,"_KEGG")
    
  WriteXLS::WriteXLS(
    sign_enrichment,
    ExcelFileName = output_file,
    SheetNames = output_names,
    AdjWidth = T,
    BoldHeaderRow = T,
    FreezeRow = 1,
    FreezeCol = 1,
    row.names = T,
    col.names = T
  )
  return(sign_enrichment)
}


dem_degs_globaltest_go <- apply_globaltest(dem_degs_expr_tps, enrichment_type="GO", designmatrix=design, output_file="results/dem_degs_globaltest_go.xls")
dem_degs_globaltest_kegg <- apply_globaltest(dem_degs_expr_tps, enrichment_type="KEGG", designmatrix=design, output_file="results/dem_degs_globaltest_kegg.xls")

dem_alltargets_globaltest_go <- apply_globaltest(dem_alltargets_expr_tps, enrichment_type="GO", designmatrix=design, output_file="results/dem_alltargets_globaltest_go.xls")
dem_alltargetss_globaltest_kegg <- apply_globaltest(dem_alltargets_expr_tps, enrichment_type="KEGG", designmatrix=design, output_file="results/dem_alltargets_globaltest_kegg.xls")

