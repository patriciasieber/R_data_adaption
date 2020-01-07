#! /home/psieber/bin/stringtie_stats/stringtie_annotation_edit+stats_geneswithintrons.R

args=commandArgs(trailingOnly=T)

## run: Rscript --vanilla /home/psieber/bin/stringtie_stats/stringtie_annotation_edit+stats_geneswithintrons.R path_to_gtf_file path_to_stringtie_annot

## read input (args)
if(length(args) == 2){
  file_name_original_annot <- args[1]
  file_name <- args[2]
}else stop("missing arguments",call.=F)

# file_name_original_annot <- "/home/psieber/rnaseq_analysis/A_fumigatus/Aspergillus_fumigatusa1163.CADRE.31.gtf"
# file_name <- "/sbidata/psieber/AlternativeSplicing/A_fumigatus/stringtie3_merge_g1.gtf"

## statistics of gene annotations

## sudo apt-get install libcurl4-openssl-dev
## sudo apt-get install libxml2-dev
#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")

library(rtracklayer)

gff <- import.gff(file_name, format="gtf")


original_annot <- import.gff(file_name_original_annot)
#original_annot <- as(original_annot, "GRanges")

types <- names(table(original_annot$type))
if(any(grepl("transcript",types)) | any(grepl("mRNA",types))){ original_mrna <- original_annot[original_annot$type == "transcript" | original_annot$type == "mRNA",]
}else stop("original annotation misses transcript or mRNA entries",call.=F)

mrna <- gff[gff$type == "transcript",]

## exclude transcripts overlapping with multiple genes ###
stringtie_transcript_id <- unique(mrna$transcript_id)
remaining_transcripts <- unlist(lapply(stringtie_transcript_id,function(x){
  current <- mrna[mrna$transcript_id == x,]
  ov <- findOverlaps(current,original_mrna)  
  ## no overlap --> new transcript predicted (not interesting for us)
  ## 1 overlap --> corresponds to exactly one gene (fits well)
  ## more than 1 overlaps --> too long transcript, discard it!
  if(length(ov) == 1) return(x)
  else return(NULL) 
}))

gff_remaining <- lapply(remaining_transcripts,function(x){
  return(gff[gff$transcript_id == x,])
})
gff_remaining <- do.call(c,gff_remaining) 

mrna_remaining <- gff_remaining[gff_remaining$type == "transcript",]


## check overlap of annotated transcripts with new transcripts, use reference gene_id
remaining_annotated_transcripts <- remaining_transcripts[!grepl("merged",remaining_transcripts)]
mrna_remaining_gene_id <- unlist(lapply(remaining_annotated_transcripts,function(x){ 
  tmp_tr <- mrna_remaining[mrna_remaining$transcript_id == x,]
  ov <- findOverlaps(tmp_tr,mrna_remaining)
  tmp_gene <- mrna_remaining[subjectHits(ov)]
  tmp_ref_gene_id <- unique(tmp_gene[!grepl("merged",tmp_gene$transcript_id)]$ref_gene_id)
  tmp_ref_gene_id <- tmp_ref_gene_id[!is.na(tmp_ref_gene_id)]
  if(length(tmp_ref_gene_id) == 1){
    tmp_gene$gene_id <- rep(tmp_ref_gene_id,length(tmp_gene))
    return(tmp_gene)
  } 
  else{
    if(length(tmp_ref_gene_id) == 0){
      return(NULL)
    }else{
      print(paste0("check ",x))
      newid <- tmp_gene[max(width(tmp_gene)) == width(tmp_gene),]$ref_gene_id
      tmp_gene$gene_id <- rep(newid,length(tmp_gene))
      return(tmp_gene)
    }
  } 
}))
mrna <- do.call(c,mrna_remaining_gene_id) 

remaining_transcript_ids <- unique(mrna$transcript_id)
gff_remaining_new <- unlist(lapply(remaining_transcript_ids,function(x){
  tmp_gff <- gff_remaining[gff_remaining$transcript_id == x,]
  new_gene_id <- unique(mrna[mrna$transcript_id == x,]$gene_id)
  if(length(new_gene_id) == 1) tmp_gff$gene_id <- rep(new_gene_id,length(tmp_gff))
  else tmp_gff$gene_id <- rep(tmp_gff$ref_gene_id[1],length(tmp_gff))
  return(tmp_gff)
}))
gff_remaining_new <- do.call(c,gff_remaining_new)

## the transcript should not just overlap with the reference range but they should share some exonic ranges 
## in order to exclude ncRNAs within the intronic regions
exons_remaining <- gff_remaining_new[gff_remaining_new$type == "exon",]
exons_remaining <- exons_remaining[!is.na(exons_remaining$gene_id),]
remaining_genes <- unique(exons_remaining$gene_id)
remaining_genes <- remaining_genes[!is.na(remaining_genes)]

remaining_transcript_ids <- unlist(lapply(remaining_genes,function(x){
  current_gene <- exons_remaining[exons_remaining$gene_id == x,]
  current_tr_ids <- unique(current_gene$transcript_id)
  current_ref_tr_id <- unique(current_tr_ids[!grepl("merged",current_tr_ids)])
  current_ref_tr_id <- current_ref_tr_id[!is.na(current_ref_tr_id)]
  if(length(current_ref_tr_id) != 1) print(paste0("there is no unique reference isoform for ",x))
  current_alt_tr_ids <- current_tr_ids[current_tr_ids != current_ref_tr_id]
  current_alt_tr_ids <- current_alt_tr_ids[!is.na(current_alt_tr_ids)]

  return_list <- list(current_ref_tr_id)
  if(length(current_alt_tr_ids) >= 1){
    for(i in 1:length(current_alt_tr_ids)){
      current_ref_tr <- current_gene[current_gene$transcript_id == current_ref_tr_id,]
      current_alt_tr <- current_gene[current_gene$transcript_id == current_alt_tr_ids[i],]
      ov <- findOverlaps(current_alt_tr,current_ref_tr)
      if(length(ov) > 0) return_list <- c(return_list,current_alt_tr_ids[i])
    }
  }
  return(return_list)
}))

gff_remaining_new <- unlist(lapply(remaining_transcript_ids,function(x){
  tmp_gff <- gff_remaining_new[gff_remaining_new$transcript_id == x,]
  return(tmp_gff)
}))
gff_remaining_new <- do.call(c,gff_remaining_new)


## consider just genes that are annotated (ignore completely new predicted genes)
export.gff(gff_remaining_new,paste0(file_name,"_remaining.gtf"),format="gtf")

## exclude exons of excluded transcripts
mrna <- gff_remaining_new[gff_remaining_new$type == "transcript",]
exons <- gff_remaining_new[gff_remaining_new$type == "exon",]


mrna_tr <- unique(mrna$transcript_id)
mrna_exons <- lapply(mrna_tr,function(x){
  return(exons[exons$transcript_id == x,])
})
mrna_exons <- do.call(c,mrna_exons) 


## annotation stats ####

ffile <- file(paste(file_name,"_stats_new",sep=""),"w")
writeLines("stats for all annotated transcripts",ffile)

genes <- unique(mrna_exons$gene_id)
writeLines(paste("number of genes:",length(genes),sep="\t"),ffile)

transcript_ids <- unique(mrna_exons$transcript_id)
writeLines(paste("number of transcripts:",length(transcript_ids),sep="\t"),ffile)

## define number of exons per transcript
writeLines(paste("number of exons:",length(mrna_exons),sep="\t"),ffile)

exons_per_transcript <- table(mrna_exons$transcript_id)
ept <- data.frame(exons_per_transcript)

## define exon length
exon_length <- width(mrna_exons)
writeLines(paste("exon length:","mean:",mean(exon_length),"standard deviation:",sd(exon_length),"min:",min(exon_length),"max:",max(exon_length),sep="\t"),ffile)

write(exon_length,paste0(file_name,"_exonLength"),ncolumns=1)

## define intron length (min, mean, max)
## exclude genes without introns
intron_table <- ept[ept$Freq > 1,]
tr_with_introns <- as.character(intron_table$Var1)
writeLines(paste("number of transcripts with more than one exon (intron-containing):",length(tr_with_introns),sep="\t"),ffile)
proportion <- length(tr_with_introns)/length(transcript_ids) * 100
writeLines(paste("proportion of intron-containing transcripts:",paste(proportion,"%",sep=" "),sep="\t"),ffile)

genes_with_introns <- unlist(lapply(tr_with_introns,function(x){
  return(mrna[mrna$transcript_id == as.character(x),]$gene_id)
}))
genes_with_introns <- unique(genes_with_introns)
proportion <- length(genes_with_introns)/length(genes) * 100
writeLines(paste("proportion of intron-containing genes:",paste(proportion,"%",sep=" "),sep="\t"),ffile)


## number of genes with more than one intron:
multiple_intron_table <- ept[ept$Freq > 2,]
mult_tr_with_introns <- as.character(multiple_intron_table$Var1)
writeLines(paste("number of transcripts with more than one intron:",length(mult_tr_with_introns),sep="\t"),ffile)
proportion <- length(mult_tr_with_introns)/length(transcript_ids) * 100
writeLines(paste("proportion of transcripts with multiple introns:",paste(proportion,"%",sep=" "),sep="\t"),ffile)

genes_with_mult_introns <- unlist(lapply(mult_tr_with_introns,function(x){
  return(mrna[mrna$transcript_id == as.character(x),]$gene_id)
}))
genes_with_mult_introns <- unique(genes_with_mult_introns)
proportion <- length(genes_with_mult_introns)/length(genes) * 100
writeLines(paste("proportion of genes with multiple introns:",paste(proportion,"%",sep=" "),sep="\t"),ffile)



intron_numbers <- intron_table$Freq
number_of_introns <- sum(intron_numbers) - length(intron_numbers)  # for each transcript one intron less than exons
writeLines(paste("number of introns:",number_of_introns,sep="\t"),ffile)

exons_with_introns <- lapply(genes_with_introns,function(x){
  return(mrna_exons[mrna_exons$gene_id == x,])
})
exons_with_introns <- do.call(c,exons_with_introns)

## define range between exons of each transcript --> introns
all_introns <- lapply(tr_with_introns,function(x){
  current_exons <- exons_with_introns[exons_with_introns$transcript_id == as.character(x),]
  return(gaps(current_exons,start=min(start(current_exons))))
})
all_introns <- do.call(c, all_introns)

intron_length <- width(all_introns)
writeLines(paste("intron length:","mean:",mean(intron_length),"standard deviation:",sd(intron_length),"min:",min(intron_length),"max:",max(intron_length),sep="\t"),ffile)
#length(intron_length)

write(intron_length,paste0(file_name,"_intronLength"),ncolumns=1)

writeLines("------------------------------------------------------",ffile)

## define AS patterns ####
writeLines("stats for genes with more than one transcript",ffile)

gene_table <- table(mrna$gene_id)
ttable <- data.frame(gene_table)

multiple_tr <- ttable[ttable$Freq > 1,]
writeLines(paste("number of genes with multiple transcripts:",nrow(multiple_tr),sep="\t"),ffile)

proportion <- nrow(multiple_tr)/length(genes) * 100
writeLines(paste("proportion of genes with multiple transcripts:",paste(proportion,"%",sep=" "),sep="\t"),ffile)

multiple_exons <- lapply(multiple_tr$Var1,function(x){
  return(mrna_exons[mrna_exons$gene_id == as.character(x),])
})
multiple_exons <- do.call(c,multiple_exons)

writeLines("------------------------------------------------------",ffile)
if(length(multiple_exons) == 0){
  writeLines("no alternative splicing patterns",ffile)
}else{
  writeLines("alternative splicing patterns",ffile)
  writeLines(paste("gene","reference transcript","alternative transcript","splicing pattern",sep="\t"),ffile)
  
  multiple_transcripts <- unique(multiple_tr$Var1)  ## genes with multiple transcripts, used for AS pattern identification
  len_multtr <- length(multiple_transcripts)
  
  altFirstExon <- 0
  altLastExon <- 0
  alt5End <- 0
  alt3End <- 0
  mxe <- 0
  intrRetention <- 0
  exonSkipping <- 0

  for(i in 1:len_multtr){
    current_gene <- as.character(multiple_transcripts[i])
    current_all <- multiple_exons[multiple_exons$gene_id == as.character(current_gene),]
    
    all_tr <- table(current_all$transcript_id)
    all_tr <- data.frame(all_tr)
    all_tr <- as.character(all_tr$Var1)
    
    if(grepl("merged",all_tr[1])){
      if(any(!grepl("merged",all_tr))){
        new_first <- all_tr[!grepl("merged",all_tr)]
        new_last <- all_tr[grepl("merged",all_tr)]
        all_tr <- c(new_first,new_last)
      }else print(paste0("annotated isoform is not the first for ",current_gene))
    }
    
    current_transcripts <- lapply(all_tr,function(x){
      output <- sort(current_all[current_all$transcript_id == x,])
      ## if two exons have the identical range, remove one
      disjoined <- disjoin(output)
      if(length(disjoined) != length(output)){
        ov <- findOverlaps(output,disjoined)
        hits <- data.frame(table(subjectHits(ov)))  
        ## we can have just multiple maps when a full exon is present multiple times
        if(any(hits$Freq > 1)){
          mult_hit <- hits[hits$Freq > 1,]$Var1
          toremove <- NULL
          for(j in 1:length(mult_hit)){
            current_hit <- mult_hit[j]
            hitsQuery <- queryHits(ov[subjectHits(ov) == current_hit,])
            toremove <- c(toremove,hitsQuery[2:length(hitsQuery)])
          }
          output <- output[-c(toremove),]
        }
        ## if different exons overlap, ? necessary ?
      }
      return(output)
    })
    strand <- as.character(strand(current_transcripts[[1]])[1])
    
    ## unique transcripts
    k <- 1
    l <- 2
    while(k < length(current_transcripts) & length(current_transcripts) > 1){
      if(length(current_transcripts[[k]]) == length(current_transcripts[[l]])){
        current_tr1 <- current_transcripts[[k]]
        current_tr2 <- current_transcripts[[l]]
        if(all(start(current_tr1) == start(current_tr2) & end(current_tr1) == end(current_tr2))){
          ## then transcript is duplicated, remove one
          current_transcripts[[k]] <- NULL
        }else{
          if(l == length(current_transcripts)){
            k <- k+1
            l <- k+1
          }else l <- l+1
        }
      }else{
        if(l == length(current_transcripts)){
          k <- k+1
          l <- k+1
        }else l <- l+1
      }
    }
    
    ## compare transcripts against each other, reference transcript against all the others
    ## we use the annotated transcript as referende
    reference_transcript <- current_transcripts[[1]]
    if(grepl("merged",reference_transcript$transcript_id[1])){
      print(paste0("annotated isoform is not the first for ",reference_transcript$transcript_id[1]))
    }
    
    if(length(current_transcripts) >= 2){
      for(j in 2:length(current_transcripts)){
        alternative_transcript <- current_transcripts[[j]]
        ## was every alternative exon checked? (in order to miss skipped exons)
        isExonChecked <- logical(length=length(alternative_transcript))

        ## Alternative First Exon
        if(strand == "+"){
          if(min(start(reference_transcript)) != min(start(alternative_transcript))){
            altFirstExon <- altFirstExon + 1
            writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative First Exon",sep="\t"),ffile)
          }
        }
        if(strand == "-"){
          if(max(end(reference_transcript)) != max(end(alternative_transcript))){
            altFirstExon <- altFirstExon + 1
            writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative First Exon",sep="\t"),ffile)
          }
        }
        ## Alternative Last Exon
        ## different last exon
        if(strand == "+"){
          if(max(end(reference_transcript)) != max(end(alternative_transcript))){
            altLastExon <- altLastExon + 1
            writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative Last Exon",sep="\t"),ffile)
          }
        }
        if(strand == "-"){
          if(min(start(reference_transcript)) != min(start(alternative_transcript))){
            altLastExon <- altLastExon + 1
            writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative Last Exon",sep="\t"),ffile)
          }
        }
      
        ## consider overlapping exons, seperately for each exon
        number_exons <- length(reference_transcript)
        for(l in 1:number_exons){
          current_exon <- reference_transcript[l]
          ov <- findOverlaps(current_exon,alternative_transcript)
          current_alternative_exons <- alternative_transcript[subjectHits(ov)]
          current_alternative_exons <- current_alternative_exons[current_alternative_exons$gene_id == as.character(current_exon$gene_id),]
          isExonChecked[subjectHits(ov)] <- TRUE
          if(length(current_alternative_exons) > 1){ ## intron retention takes place, still more?
            intrRetention <- intrRetention + 1
            writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Intron Retention",sep="\t"),ffile)
          }else{
            if(length(current_alternative_exons) == 0){
              range_to_neighbor_exons <- current_exon
              if(l != 1 && l != number_exons){
                start(range_to_neighbor_exons) <- end(reference_transcript[l-1])+1
                end(range_to_neighbor_exons) <- start(reference_transcript[l+1])-1
                ov <- findOverlaps(range_to_neighbor_exons,alternative_transcript)
                mxe_candidate <- alternative_transcript[subjectHits(ov)]
                mxe_candidate <- mxe_candidate[mxe_candidate$gene_id == as.character(current_exon$gene_id),]
                if(length(mxe_candidate) > 0){
                  ## check for no overlap back to reference transcript
                  ov <- findOverlaps(mxe_candidate,reference_transcript)
                  len_ov <- length(reference_transcript[subjectHits(ov)])
                  if(len_ov == 0){ ## mxe takes place
                    mxe <- mxe+1
                    writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Mutually Exclusive Exon",sep="\t"),ffile)
                  }else{ 
                    ## check whether the non-overlap is outside the range of the other transcript, then it's already considered as alternative start/end
                    if(end(current_exon) > start(alternative_transcript) && start(current_exon) < end(alternative_transcript)){
                      ## it's exon skipping with additional alternative SS (detected in later or earlier step)
                      exonSkipping <- exonSkipping+1
                      writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Exon Skipping",sep="\t"),ffile)
                      print(paste0("SE?!","\t",current_gene,"\t",reference_transcript$transcript_id[1],"\t",alternative_transcript$transcript_id[1],"; i=",i,"; j=",j))
                    }

                  }
                }## else: area is not in the range of the alternative transcript, already considered by alternative start / end
              }## else: alternative start or end, and then it's already considered 
            }else{  ## length(current_alternative_exons) == 1
              if(range(current_exon) != range(current_alternative_exons)){ ## else go to the next exon
                if(start(current_exon) != start(current_alternative_exons)){
                  if(l != 1){ ## otherwise already considered as first/last exon
                    if(start(current_alternative_exons) < end(reference_transcript[l-1])){
                      ## IR takes place, already counted in the previous step
                    }else{
                      ## Alternative 3'SS
                      if(strand == "+"){
                        alt3End <- alt3End + 1
                        writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative 3'End",sep="\t"),ffile)
                      }
                      if(strand == "-"){
                        alt5End <- alt5End + 1
                        writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative 5'End",sep="\t"),ffile)
                      }
                    }
                  }
                }
                if(end(current_exon) != end(current_alternative_exons)){
                  if(l != number_exons){ ## otherwise already considered as first/last exon
                    ## check overlap of alternative exon with other exons, to determine whether A5SS or IR is present
                    ov <- findOverlaps(current_alternative_exons,reference_transcript)
                    ir_test <- reference_transcript[subjectHits(ov)]
                    if(length(ir_test) > 1){ ## IR takes place
                      intrRetention <- intrRetention + 1
                      writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Intron Retention",sep="\t"),ffile)
                    }else{
                      if(strand == "+"){
                        alt5End <- alt5End + 1
                        writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative 5'End",sep="\t"),ffile)
                      }
                      if(strand == "-"){
                        alt3End <- alt3End + 1
                        writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Alternative 3'End",sep="\t"),ffile)
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if(sum(isExonChecked) < length(alternative_transcript)){
          unchecked_position <- which(isExonChecked==F)
          ## check if missing exon overlaps with any reference exon
          ## if not --> skipped exon
          for(m in unchecked_position){
            current_alternative_exon <- alternative_transcript[m]
            ov <- findOverlaps(current_alternative_exon,reference_transcript)
            if(length(ov) == 0){ ## SE or overlapping start/end; MXE was checked already
              if(end(current_alternative_exon) > min(start(reference_transcript)) && start(current_alternative_exon) < max(end(reference_transcript))){  
                ## else would be overlapping start/end, already considered
                ## now just internal exons without overlap considered
                if(m != 1 && m != length(isExonChecked)){ ## otherwise it is alternative start/end and already considered
                  exonSkipping <- exonSkipping+1
                  writeLines(paste(current_gene,reference_transcript$transcript_id[1],alternative_transcript$transcript_id[1],"Exon Skipping",sep="\t"),ffile)
                }
              }
            }else ## else: TODO check!
              print(paste0("SE?!","\t",current_gene,"\t",reference_transcript$transcript_id[1],"\t",alternative_transcript$transcript_id[1],"; i=",i,"; j=",j,"; m=",m))
          }
        }
      }
    }
  }
  writeLines("------------------------------------------------------",ffile)
  writeLines("AS summary:",ffile)
  writeLines(paste("Alternative First Exon:",altFirstExon,sep="\t"),ffile)
  writeLines(paste("Alternative Last Exon:",altLastExon,sep="\t"),ffile)
  writeLines(paste("Alternative 5'End:",alt5End,sep="\t"),ffile)
  writeLines(paste("Alternative 3'End:",alt3End,sep="\t"),ffile)
  writeLines(paste("Mutually Exclusive Exon:",mxe,sep="\t"),ffile)
  writeLines(paste("Exon Skipping:",exonSkipping,sep="\t"),ffile)
  writeLines(paste("Intron Retention:",intrRetention,sep="\t"),ffile)
}

close(ffile)

