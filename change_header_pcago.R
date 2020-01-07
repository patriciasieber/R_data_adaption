## change header for pca plot

setwd("/data/")


## add metadata to header (patient information)
metadata <- read.csv("PatientDatabasehBALED.csv",sep="\t")

header <- read.table("header.txt",sep="\t")
n <- ncol(header)
new_header <- vector(length=n)

for(i in 1:n){
  input <- as.character(header[1,i])
  pat_number <- strsplit(strsplit(input," ")[[1]][1],"_")[[1]][2]
  if(substr(pat_number,1,1) == 0) pat_number <- substr(pat_number,2,nchar(pat_number))
  
  pat_metadata <- metadata[metadata$PATIENT.ID == pat_number,]
  
  ## add age, sex, disease, steroids as additional metadata into the header
  disease <- gsub(" ","",as.character(pat_metadata$Disease))
  steroids <- gsub(" ","",as.character(pat_metadata$steroids.IS))
  output <- paste(input,as.character(pat_metadata$Age),as.character(pat_metadata$Sex),disease,steroids,sep=" ")
  
  new_header[i] <- output
}


write.table(new_header,"header_new.txt",sep="\t",quote = F,row.names = F)
