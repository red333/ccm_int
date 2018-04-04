library(bio3d)


s2f <- function(id, chain, move=1){
  id <- as.character(id)
  chain <- as.character(chain)
  read.pdb(id) -> x
  x <- x$atom
  cord <- x[, c("elety","chain","resno","x","y","z")]
  cordCA <- cord[cord$elety =='CA',]
  cordMatrix <- cordCA[cordCA$chain == chain,c("resno","x","y","z")]
  #chain name gonna be para from input
  M <- matrix(unlist(cordMatrix),ncol=4)
  colnames(M) <- c("R","X","Y","Z")
  
  #frag library Matrix
  N <- read.table("lib_10_z_5.txt",header = T)
  N <- matrix(unlist(N),ncol=3)
  colnames(N) <- c("X","Y","Z")
  if(nrow(M)-4 > 0){
  P = matrix(0,ncol = length(seq(from = 1, to = (nrow(M)-4), by = move)), nrow = 10)
  count = 1
  for(i in seq(from = 1, to = (nrow(M)-4), by = move)){
    resno <- M[i:(i+4),1]
    problist <- list()
    if(resno[5]-resno[1] == 4){
    #tmp5 <- M[i:(i+4),2:4]
    for(j in 1:10){
      #rot = rot.lsq(tmp5,N[(5*j-4):(5*j),])
      P[j,count] = rmsd(c(M[i:(i+4),2:4]), c(N[(5*j-4):(5*j),]), fit = TRUE)
      ##fit = T in rmsd to see if its different
    }
    }else{problist = append(problist,count)}
    count = count + 1
  }
  fragSeq <- apply(P,2,which.min)
  fragSeq <- fragSeq - 1
  if(length(problist) > 0){
    for(i in 1:length(problist)){
      fragSeq[as.numeric(problist[i])] = "L"
    }
  }
  fragSeq <- paste0(fragSeq)
  fragSeq <- chartr("0123456789", "ARNDCQEGHI", fragSeq)
  
  cat(paste0(">",id," | Red333\n"), file = paste0(id,chain,".fas"), sep="")
  cat(fragSeq, file = paste0(id,chain,".fas"),sep = "",append = TRUE)
  }
  if(nrow(M)-4 <= 0){
    cat(paste0(id,"_",chain), file = "lessfourchain.txt", sep="\n", append = T)
  }
  cat(paste0(id), file = "worklog.txt", sep="\n", append = T)	
}








  
