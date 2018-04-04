#!/usr/bin/Rscript
step2 <- function(targetstart, querystart, mapfasta, alignmatrix){
  linkmap <- read.fasta(mapfasta,rm.dup = T, to.upper = T, to.dash = T)
  # targetstart <- 3
  # querystart <- 18
  target <- alignmatrix[1,]
  linkmap <- as.data.frame(linkmap$ali, stringsAsFactors = F)
  linkindex <- as.data.frame(linkmap,stringsAsFactors = F)
  targettrack <- 1
  target[2,] <- "0"
  
  
  for(i in 1:ncol(target)){
    if(target[1,i] != "." && target[1,i] != "-"){
      if(linkmap[2,targettrack] != "." && linkmap[2,targettrack] != "-"){
        if(linkmap[1,targettrack] != "." && linkmap[1,targettrack] != "-"){
          target[2,i] <- querystart
          querystart <- querystart + 1
          targetstart <- targetstart + 1
          targettrack <- targettrack + 1
        }else{
          targetstart <- targetstart + 1
          targettrack <- targettrack + 1
        }
      }else{
        while(linkmap[2,targettrack] == "." | linkmap[2,targettrack] == "-"){
          if(linkmap[1,targettrack] == "." | linkmap[1,targettrack] == "-"){
            
          }else{
            querystart <- querystart + 1
          }
          
          targettrack <- targettrack + 1
        }
        if(linkmap[1,targettrack] != "." && linkmap[1,targettrack] != "-"){
          target[2,i] <- querystart
          querystart <- querystart + 1
          targetstart <- targetstart + 1
          targettrack <- targettrack + 1
        }else{
          targetstart <- targetstart + 1
          targettrack <- targettrack + 1
        }
        
      }
    }
  }
  
  # linkmap2 <- read.fasta("PYRI_4fyyB_w_1m0kA.link.fasta",rm.dup = T, to.upper = T, to.dash = T)
  # targetstart2 <- 14
  # querystart2 <- 7
  # target2 <- alignmatrix_2[3,]
  # linkmap2 <- as.data.frame(linkmap2$ali, stringsAsFactors = F)
  # linkindex2 <- as.data.frame(linkmap2,stringsAsFactors = F)
  # targettrack2 <- 1
  # target2[2,] <- "0"
  # 
  # 
  # for(i in 1:ncol(target2)){
  #   if(target2[1,i] != "." && target2[1,i] != "-"){
  #     if(linkmap2[2,targettrack2] != "." && linkmap2[2,targettrack2] != "-"){
  #       if(linkmap2[1,targettrack2] != "." && linkmap2[1,targettrack2] != "-"){
  #         target2[2,i] <- querystart2
  #         querystart2 <- querystart2 + 1
  #         targetstart2 <- targetstart2 + 1
  #         targettrack2 <- targettrack2 + 1
  #       }else{
  #         targetstart2 <- targetstart2 + 1
  #         targettrack2 <- targettrack2 + 1
  #       }
  #     }else{
  #       while(linkmap2[2,targettrack2] == "." | linkmap2[2,targettrack2] == "-"){
  #         if(linkmap2[1,targettrack2] == "." | linkmap2[1,targettrack2] == "-"){
  #           
  #         }else{
  #           querystart2 <- querystart2 + 1
  #         }
  #         
  #         targettrack2 <- targettrack2 + 1
  #       }
  #       if(linkmap2[1,targettrack2] != "." && linkmap2[1,targettrack2] != "-"){
  #         target2[2,i] <- querystart2
  #         querystart2 <- querystart2 + 1
  #         targetstart2 <- targetstart2 + 1
  #         targettrack2 <- targettrack2 + 1
  #       }else{
  #         targetstart2 <- targetstart2 + 1
  #         targettrack2 <- targettrack2 + 1
  #       }
  #       
  #     }
  #   }
  # }
  
  return (target)
}