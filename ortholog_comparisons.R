## which proteins appear in both species to be differentially expressed?

setwd("/data/")

degs_human <- read.csv("data/human_signproteins.txt",header=T,sep="\t",row.names=NULL)
degs_human <- as.character(degs_human$Gene.names)

degs_mouse <- read.csv("data/mouse_signproteins.txt",header=T,sep="\t")
degs_mouse <- as.character(degs_mouse$Gene.names)

biomart <- read.csv("data/mart_export_orthologs.txt",header=T)

## compare gene names (protein IDs does not work with biomart data)
## translate mouse genes into human ones
degs_mouse2humans <- unlist(lapply(degs_mouse,function(x){
  current <- biomart[biomart$Mouse.gene.name == x,]
  current_name <- as.character(unique(current$Gene.name))
  if(length(current_name) == 1) return(current_name)
  else{
    if(length(current_name) > 1) return(current_name[1])#return(paste(current_name,collapse=";"))  ## --> both have the same results
    else return("")
  } 
}))

write(degs_mouse2humans,"mouse_deps_humanorthologs.txt")



## try again with new biomart data and both directions (human to mouse, mouse to human)


biomart_human <- read.csv("data/mart_export_orthologs_20191114_human.txt",header=T,sep="\t")
biomart_mouse <- read.csv("data/mart_export_orthologs_20191114_mouse.txt",header=T,sep="\t")

## compare gene names (protein IDs does not work with biomart data)
## translate mouse genes into human ones
degs_human2mouse <- unlist(lapply(degs_human,function(x){
  current <- biomart_human[biomart_human$Gene.name == x,]
  current_name <- as.character(unique(current$Mouse.gene.name))
  if(length(current_name) == 1) return(current_name)
  else{
    if(length(current_name) > 1) return(current_name[1])#return(paste(current_name,collapse=";"))  ## --> both have the same results
    else return("")
  } 
}))

degs_mouse2humans <- unlist(lapply(degs_mouse,function(x){
  current <- biomart_mouse[biomart_mouse$Gene.name == x,]
  current_name <- as.character(unique(current$Human.gene.name))
  if(length(current_name) == 1) return(current_name)
  else{
    if(length(current_name) > 1) return(current_name[1])#return(paste(current_name,collapse=";"))  ## --> both have the same results
    else return("")
  } 
}))


human2mouse <- cbind(degs_human,degs_human2mouse)
mouse2human <- cbind(degs_mouse,degs_mouse2humans)

# find overlap between human and mouse, look from both sites
ov_human2mouse <- intersect(human2mouse[,1],mouse2human[,2])
ov_mouse2human <- intersect(mouse2human[,1],human2mouse[,2])
#--> is acutally different!?!


