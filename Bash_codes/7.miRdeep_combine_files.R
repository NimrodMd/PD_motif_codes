#!/usr/bin/env Rscript
setwd("Y:/nimrod.madrer/Folder_name/miRDeep2")

# combine seperate miRDeep2 files to one matrix
folderList <- list.dirs(path = getwd(),full.names = F,recursive = F)
for(f in folderList){
  #f = "1-2_S1_full.trimmed"
  sname <- strsplit(f,"_R1")[[1]][1]
  sfile <- read.delim(paste(f,"miRNAs_expressed_all_samples_16_19.csv",sep = "/"),header = T,sep = "\t")
  
  if(f == folderList[1]){
    unnormcounts <- sfile[,c("X.miRNA","precursor","read_count")]
    colnames(unnormcounts)[which(colnames(unnormcounts) == "read_count")] <- sname
    RPM <- sfile[,c("X.miRNA","precursor","seq.norm.")]
    colnames(RPM)[which(colnames(RPM) == "seq.norm.")] <- sname
  } else {
    if(all(unnormcounts$X.miRNA == sfile$X.miRNA) & all(unnormcounts$precursor == sfile$precursor)){
      print("all equal, continuing")
      unnormcounts <- cbind(unnormcounts,sfile[,c("read_count")])
      colnames(unnormcounts)[which(colnames(unnormcounts) == "sfile[, c(\"read_count\")]")] <- sname
      
      RPM<- cbind(RPM,sfile[,c("seq.norm.")])
      colnames(RPM)[which(colnames(RPM) == "sfile[, c(\"seq.norm.\")]")] <- sname
    }
    
  }
  print(paste0(which(folderList == f), " ", f))
}

write.csv(RPM,"RPM.csv")
write.csv(unnormcounts,"unnormcounts.csv")


### filter repeating miRnames (from different precursors - taking the minimum count for each)
dup_unnormcounts <-unique(unnormcounts$X.miRNA[duplicated(unnormcounts$X.miRNA)])
unnormcounts_compiled <- unnormcounts[which(!(unnormcounts$X.miRNA %in% dup_unnormcounts)),]

mm=1
for(m in dup_unnormcounts){
  #m = dup_unnormcounts[1]
  temp_dup_miR <- unnormcounts[which(unnormcounts$X.miRNA == m),]
  temp_precursor <- NULL
  for(p in 1:length(temp_dup_miR$precursor)){
    if(p == 1){
      temp_precursor <- temp_dup_miR$precursor[p]
    }else{
      temp_precursor <- paste0(temp_precursor,",",temp_dup_miR$precursor[p])
    }
  }
  
  temp_counts <- apply(temp_dup_miR[3:length(temp_dup_miR[1,])],MARGIN = 2, function(x){min(as.numeric(x))})
  temp_line <- cbind(as.character(temp_dup_miR$X.miRNA[1]),temp_precursor,t(temp_counts))
  colnames(temp_line)[c(1,2)] <- c("X.miRNA","precursor")
  
  unnormcounts_compiled <- rbind(unnormcounts_compiled,temp_line)
  #dup_unnormcounts <- dup_unnormcounts[-c(which(dup_unnormcounts == m))]
  print(mm)
  mm = mm+1
}
write.csv(unnormcounts_compiled,"unnormcounts_compiled.csv")



