## Patricia, 2019_06_11
## filter miRNA results interesting for us, sort conditions
library(gdata)
library(ggplot2)
library(gridExtra)

setwd("/data/Postdoc/miRNA_targets/ProjectSummary/Data/") # ran on local machine

results <- read.csv("normalized_dcq_values_ref_20653.csv",sep="\t",header=T)

sampleSelection <- read.xls("../../data/sampleSelection.xlsx",sheet=2)
sampleSelection <- sampleSelection[,1:7]

## sort for interesting candidate genes
mirnas <- c("hsa-miR-204-5p","hsa-miR-132-3p","hsa-miR-9-5p")
dems <- c(## DE miRNAs male
          "hsa-miR-194-5p","hsa-miR-376a-3p","hsa-miR-215-5p","hsa-miR-495-3p","hsa-miR-382-5p","hsa-miR-151a-5p",  
          "hsa-miR-182-5p","hsa-miR-151a-3p","hsa-miR-30d-5p","hsa-miR-23b-3p","hsa-miR-33a-5p","hsa-miR-223-3p",
          "hsa-miR-766-3p","hsa-miR-324-3p","hsa-miR-340-5p","hsa-miR-27b-3p","hsa-miR-339-5p","hsa-miR-409-3p",
          "hsa-miR-154-5p","hsa-miR-192-5p",
          ## DE miRNAs female
          "hsa-miR-21-5p","hsa-miR-342-3p","hsa-miR-484","hsa-miR-150-5p","hsa-miR-191-5p","ssa-miR-103a-3p","hsa-miR-181a-5p",
          "hsa-miR-222-3p","hsa-miR-17-5p","hsa-miR-328-3p","hsa-miR-106a-5p","hsa-miR-155-5p","hsa-miR-20a-5p","hsa-let-7d-5p",
          "hsa-miR-18b-5p","hsa-miR-15b-5p","hsa-miR-30e-5p","hsa-miR-107","hsa-miR-140-3p","hsa-miR-663a")

mirnas_pos <- unlist(lapply(mirnas,function(x){grep(x,results$dilution_id)}))
dems_pos <- unlist(lapply(dems,function(x){grep(x,results$dilution_id)}))
rest_pos <- 4:nrow(results)
rest_pos_rm <- unlist(lapply(c(mirnas_pos,dems_pos),function(x){which(x == rest_pos)}))
rest_pos <- rest_pos[-rest_pos_rm]

results <- results[c(1:3,mirnas_pos,dems_pos,rest_pos),]

dilution_name <- as.character(unlist(lapply(results[1,],as.character)))
patient_name <- dilution_name[-1]

association_table <- data.frame(patient=c("P043","P010","P055","P015"),associated_control=c("P018","None","P063","P053"))


## sort for patients in new data.frame:
## patient_id; timepoint; case; association; patient_name; miRNAs...
patient_id <- unlist(lapply(patient_name,function(x){
  strsplit(x,"-")[[1]][1]
}))

timepoint <- unlist(lapply(patient_name,function(x){
  sampleSelection[sampleSelection$Patient.Sample == x,]$postTX.d.
}))

case <- unlist(lapply(patient_name,function(x){
  temp <- as.character(sampleSelection[sampleSelection$Patient.Sample == x,]$Status)
  strsplit(temp," ")[[1]][2]
}))

association <- unlist(lapply(patient_name,function(x){
  temp <- strsplit(x,"-")[[1]][1]
  associate <- as.character(association_table[association_table$patient == temp,]$associated_control)
  if(length(associate) == 0) associate <- as.character(association_table[association_table$associated_control == temp,]$patient)
  return(associate)
}))

mirnas_patient_table <- data.frame(patient_id=patient_id,timepoint=timepoint,case=case,association=association,patient_name=patient_name)


## enter miRNA dcq values and sort by patients and timepoints
results_sort <- results[-c(1:3),]
colnames(results_sort) <- c("dilution_name",patient_name)
mir_names <- as.character(results_sort$dilution_name)
for(mi in 1:length(mir_names)){
  current_mirna <- mir_names[mi]
  current_values <- as.numeric(unlist(lapply(results_sort[mi,2:ncol(results_sort)],as.character)))
  mirnas_patient_table$current_mirna <- current_values
  colnames(mirnas_patient_table)[ncol(mirnas_patient_table)] <- current_mirna
}

mirnas_patient_table <- mirnas_patient_table[order(mirnas_patient_table$patient_id, mirnas_patient_table$timepoint), ]

write.table(mirnas_patient_table,"normalized_dcq_values_patients.csv",sep="\t",col.names=T,row.names=F,quote=F)
save(mirnas_patient_table,file="normalized_dcq_values_patients.RData")




## line plots
plot_mirnas <- c(mirnas,dems) ## c("hsa-miR-204-5p","hsa-miR-132-3p","hsa-miR-9-5p")
len_mir <- length(plot_mirnas)

plot_patient_pairs <- c("P043_P018","P055_P063","P015_P053","P010_None")
len_pat_pairs <- length(plot_patient_pairs)
yrange <- range(-8,6)

for(i in 1:len_mir){
  current_mirna <- plot_mirnas[i]
  mirnas_patient_table_currentmirna <- mirnas_patient_table[,c(1:5,grep(current_mirna,colnames(mirnas_patient_table)))]
  
  if(ncol(mirnas_patient_table_currentmirna) == 6){
    pltList <- NULL
    for(j in 1:len_pat_pairs){
      current_patient_pair <- plot_patient_pairs[j]
      pat <- strsplit(current_patient_pair,"_")[[1]]
      mirnas_patient_table_currentmirna_currentpatientpair <- mirnas_patient_table_currentmirna[mirnas_patient_table_currentmirna$patient_id == pat[1] | mirnas_patient_table_currentmirna$patient_id == pat[2],]
      mirnas_patient_table_currentmirna_currentpatientpair <- mirnas_patient_table_currentmirna_currentpatientpair[!is.na(mirnas_patient_table_currentmirna_currentpatientpair[,6]),]
      mirnas_patient_table_currentmirna_currentpatientpair$patient_id <- factor(mirnas_patient_table_currentmirna_currentpatientpair$patient_id, levels=c(pat[1],pat[2]))
      
      pltList[[j]] <- ggplot(data=mirnas_patient_table_currentmirna_currentpatientpair, 
                             aes(x=timepoint,y=get(current_mirna),group=patient_id)) +
        geom_line(aes(color=patient_id)) + 
        geom_point(aes(color=patient_id)) +
        scale_color_brewer(palette="Paired") +
        labs(title=paste0(current_mirna,": Patient=",pat[1],"; Control=",pat[2]),x="time points", y = "dcq values") +
        ylim(yrange) 
      
    }
    
    pdf(paste0("dcq_",current_mirna,".pdf"),width=20,height=4,onefile=F)
    print(marrangeGrob(pltList, nrow=1, ncol=len_pat_pairs))
    dev.off()
  }## else: not listed in csv
  
}


