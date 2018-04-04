alignmatrix_pair <- cbind.data.frame(alignmatrix[1:min(nrow(alignmatrix),nrow(alignmatrix_2)),],alignmatrix_2[1:min(nrow(alignmatrix),nrow(alignmatrix_2)),],deparse.level = 0, stringsAsFactors = F)
gapbinary_pair <- cbind.data.frame(gapbinary[1:min(nrow(gapbinary),nrow(gapbinary_2)),],gapbinary_2[1:min(nrow(gapbinary),nrow(gapbinary_2)),],deparse.level = 0, stringsAsFactors = F)
interfacebinary_pair <- cbind.data.frame(interfacebinary[1:min(nrow(interfacebinary),nrow(interfacebinary_2)),],interfacebinary_2[1:min(nrow(interfacebinary),nrow(interfacebinary_2)),],deparse.level = 0, stringsAsFactors = F)
target_pair <- cbind.data.frame(target[1:min(nrow(target),nrow(target2)),],target2[1:min(nrow(target),nrow(target2)),],deparse.level = 0, stringsAsFactors = F)
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

drop1 <- which(queryindex[2,] == "0")

alignmatrix_pair <- alignmatrix_pair[,-drop1]
gapbinary_pair <- gapbinary_pair[,-drop1]
interfacebinary_pair <- interfacebinary_pair[,-drop1]
queryindex <- queryindex[,-drop1]

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

