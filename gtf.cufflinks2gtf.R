library(rtracklayer)

dir <- "/home/psieber/rnaseq_analysis/C_albicans"
setwd(dir)
gff <- import("transcripts_cufflinks_C.albicans.gtf")

strand <- strand(gff)

new_strand <- lapply(as.character(strand),function(x){
  if(x == "*") return("+")
  else return(x)
})

strand(gff) <- unlist(new_strand)

ffile <- file("transcripts_cufflinks_C.albicans.plus.gtf")
writeLines(paste(seqnames(gff),gff$source,gff$type,start(gff),end(gff),gff$score,strand(gff),".",
                 paste0("gene_id \"",gff$gene_id,"\"; transcript_id \"",gff$transcript_id,"\";"),sep="\t"),ffile)
close(ffile)