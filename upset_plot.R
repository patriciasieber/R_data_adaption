## create upset plot for overlapping lists
## https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html

library(UpSetR)

files_for_upset <- list.files(path="subsets/",full.names=T)

sets <- lapply(files_for_upset,readLines)

sets_names <- sapply(files_for_upset,function(x){
  splits <- strsplit(strsplit(x,"//")[[1]][2],"_")[[1]]
  return(paste0(substr(splits[1],1,1),substr(splits[2],1,2)))
})
names(sets) <- sets_names


svg(filename="upsetr.svg", 
    width=14, 
    height=12, 
    pointsize=12)
upset(fromList(sets), nsets = 20, nintersects = 40, sets = c("Ctr","Cmu","Csu","Cpn","Cpe","Cib","Cga","Cav","Cfe","Cca","Cps","Cab"), # sets=names(sets)
    mainbar.y.label = "No. of common genes", sets.x.label = "No. of annotated genes", 
    order.by = "freq", sets.bar.color = "#56B4E9", keep.order = T, 
    text.scale = 1.4, point.size = 2.6, line.size = 0.8,
    set_size.show = TRUE, empty.intersections = "off")
dev.off()
