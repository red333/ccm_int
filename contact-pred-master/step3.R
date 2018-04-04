#!/usr/bin/Rscript
step3 <- function(alignment_l, alignment_r, interface_l, interface_r, gap_l, gap_r, queryindex_l, queryindex_r){
  alignmatrix_pair <- cbind.data.frame(alignment_l[1:min(nrow(alignment_l),nrow(alignment_r)),],alignment_r[1:min(nrow(alignment_l),nrow(alignment_r)),],deparse.level = 0, stringsAsFactors = F)
  gapbinary_pair <- cbind.data.frame(gap_l[1:min(nrow(gap_l),nrow(gap_r)),],gap_r[1:min(nrow(gap_l),nrow(gap_r)),],deparse.level = 0, stringsAsFactors = F)
  interfacebinary_pair <- cbind.data.frame(interface_l[1:min(nrow(interface_l),nrow(interface_r)),],interface_r[1:min(nrow(interface_l),nrow(interface_r)),],deparse.level = 0, stringsAsFactors = F)
  target_pair <- cbind.data.frame(queryindex_l[1:min(nrow(queryindex_l),nrow(queryindex_r)),],queryindex_r[1:min(nrow(queryindex_l),nrow(queryindex_r)),],deparse.level = 0, stringsAsFactors = F)
  queryindex <- target_pair
  
  #only keep known index
  
  # colnum <- ncol(queryindex)
  # for(i in 1:colnum){
  # if(queryindex[2,i] == "0"){
  #   alignmatrix_pair <- alignmatrix_pair[,-i]
  #   gapbinary_pair <- gapbinary_pair[,-i]
  #   interfacebinary_pair <- interfacebinary_pair[,-i]
  #   queryindex <- queryindex[,-i]
  #   colnum <- colnum - 1
  #   i <- i - 1
  # }
  # 
  # }
  
  #cutoff the percentage of interface for columns para
  cut_interfacepersentage <- .2
  
  #drop no_align positions
  drop1 <- which(queryindex[2,] == "0")
  
  alignmatrix_pair <- alignmatrix_pair[,-drop1]
  gapbinary_pair <- gapbinary_pair[,-drop1]
  interfacebinary_pair <- interfacebinary_pair[,-drop1]
  queryindex <- queryindex[,-drop1]
  
  #gap_cutoff
  gappercentage <- c()
  for(i in 1:ncol(gapbinary_pair)){
    gappercentage[i] <- length(which(gapbinary_pair[,i] == "0")) / nrow(gapbinary_pair)
  }
  drop2 <- which(gappercentage > .75)
  alignmatrix_pair <- alignmatrix_pair[,-drop2]
  gapbinary_pair <- gapbinary_pair[,-drop2]
  interfacebinary_pair <- interfacebinary_pair[,-drop2]
  queryindex <- queryindex[,-drop2]
  
  interfacepersentage <- c()
  for(i in 1:ncol(interfacebinary_pair)){
    interfacepersentage[i] <- length(which(interfacebinary_pair[,i] == "1")) / (length(which(interfacebinary_pair[,i] == "1")) + length(which(interfacebinary_pair[,i] == "0")))
  }
  drop3 <- which(interfacepersentage < cut_interfacepersentage)
  alignmatrix_pair <- alignmatrix_pair[,-drop3]
  gapbinary_pair <- gapbinary_pair[,-drop3]
  interfacebinary_pair <- interfacebinary_pair[,-drop3]
  queryindex <- queryindex[,-drop3]
  
  return(list(alignmatrix_pair, interfacebinary_pair, gapbinary_pair, queryindex))
}