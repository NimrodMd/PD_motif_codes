# First we want to merge the summary files

# Exclusive Alignment summary
setwd("Y:/nimrod.madrer/Folder_name/MINTmap_output") 

library("data.table") ; library(stringr)
file_list <- list.files(pattern = "exclusive-tRFs.countsmeta.txt")

# Read the files in, assuming comma separator
txt_files_df <- lapply(file_list, function(x) {read.table(file = x, header = T, sep ="\t")})
# Combine them
combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
rownames(combined_df) <- file_list
ExclusiveSummary <- combined_df

write.csv(ExclusiveSummary, "ExclusiveSummaryFile.csv")

# Ambiguous Alignment Summary
file_list <- list.files(pattern = "ambiguous-tRFs.countsmeta.txt")
# Read the files in, assuming comma separator
txt_files_df <- lapply(file_list, function(x) {read.table(file = x, header = T, sep ="\t")})
# Combine them
combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
rownames(combined_df) <- file_list
AmbiguousSummary <- combined_df

write.csv(AmbiguousSummary, "AmbiguousSummaryFile.csv")


# output files from MINTmap are .txt files with tons of information
# we want to only extract the information on the raw reads to process the files using DeSeq2
# use freads for reading multiple files
# use a for loop for grapping the files and putting them together

file_list <- list.files(pattern = "exclusive-tRFs.expression.txt")
txt_files_df <- lapply(file_list, function(x) {fread(file=x, select=c(1,4), data.table=F, col.names=c("MINTbase Unique ID",x))})
combined_df <-Reduce(function(x,y) {merge(x,y, by = "MINTbase Unique ID", all = TRUE)}, txt_files_df)
rownames(combined_df) <- combined_df$`MINTbase Unique ID`  ;  combined_df <- combined_df[,-1]
colnames(combined_df) <- gsub("-MINTmap_v1-exclusive-tRFs.expression.txt", "", colnames(combined_df))
combined_df[is.na(combined_df)] <- 0  ;  combined_df[1:5,1:5]
write.csv(combined_df, "tRNA_Exclusive_Combined_data.csv")

#RPM file:
file_list <- list.files(pattern = "exclusive-tRFs.expression.txt")
txt_files_df <- lapply(file_list, function(x) {fread(file=x, select=c(1,5), data.table=F, col.names=c("MINTbase Unique ID",x))})
combined_df <-Reduce(function(x,y) {merge(x,y, by = "MINTbase Unique ID", all = TRUE)}, txt_files_df)
rownames(combined_df) <- combined_df$`MINTbase Unique ID`  ;  combined_df <- combined_df[,-1]
colnames(combined_df) <- gsub("-MINTmap_v1-exclusive-tRFs.expression.txt", "", colnames(combined_df))
combined_df[is.na(combined_df)] <- 0  ;  combined_df[1:5,1:5]
write.csv(combined_df, "tRNA_Exclusive_Combined_data_RPM.csv")


#tRF meta information: sequence, type, amino acid
file_list <- list.files(pattern = "exclusive-tRFs.expression.txt")
txt_files_df <- lapply(file_list, function(x) {fread(file = x, select = c(1:3,8), data.table = FALSE)})
combined_df <- plyr::ldply(txt_files_df)
combined_df <- combined_df[order(combined_df$`MINTbase Unique ID`),]
combined_df <- combined_df[!duplicated(combined_df$`MINTbase Unique ID`),]
colnames(combined_df)<-gsub('Sequence.locations.in.tRNA.space..comma.deliminated.','details',colnames(combined_df))
combined_df$len<-unlist(lapply(combined_df$`MINTbase Unique ID`, function(x) strsplit(as.character(x),'-')[[1]][2]))
colnames(combined_df)<-gsub('MINTbase.Unique.ID','trf',colnames(combined_df))
combined_df$trna<-unlist(lapply(combined_df$details, function(x) substr(strsplit(as.character(x),'_')[[1]][2],1,3)))
combined_df$nuclear<-factor(unlist(lapply(combined_df$details,function(x) sum(str_count(as.character(as.character(x)),as.character('trnaMT')))!=0)),
                          labels=c('Nuclear','Mitochondrial'))
rownames(combined_df) <- combined_df$trf
head(combined_df) ; dim(combined_df)
write.csv(combined_df, "tRF_meta.csv")

# merge ambiguous reads
file_list <- list.files(pattern = "ambiguous-tRFs.expression.txt")
txt_files_df <- lapply(file_list, function(x) {fread(file = x, select = c(1,4), data.table=F, col.names=c("MINTbase Unique ID",x))})
combined_df <-Reduce(function(x,y) {merge(x,y, by = "MINTbase Unique ID", all = TRUE)}, txt_files_df)
rownames(combined_df) <- combined_df$`MINTbase Unique ID`
combined_df <- combined_df[,-1]
colnames(combined_df) <- gsub("-MINTmap_v1-ambiguous-tRFs.expression.txt", "", file_list)
combined_df[is.na(combined_df)] <- 0  ;  combined_df[1:5,1:5]
write.csv(combined_df, "tRNA_ambiguous_Combined_data.csv")
  


