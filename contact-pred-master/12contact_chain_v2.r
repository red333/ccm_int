#!/usr/bin/Rscript
library(bio3d)

sss <- function(p1,p2){
  if(file.exists(p1) && file.exists(p2)){
    p1p <- read.pdb(p1)
    p2p <- read.pdb(p2)
    p1p <- p1p$atom
    p2p <- p2p$atom
    p1p = p1p[p1p$type == "ATOM",]
    p2p = p2p[p2p$type == "ATOM",]
    p2p = p2p[,c("resno","x","y","z")]
    p1p = p1p[,c("resno","x","y","z")]
    xx <- matrix(data = 0, ncol = 2)
    for(j in 1:nrow(p1p)){
      for(k in 1:nrow(p2p)){
        if(dist.xyz(p1p[j,2:4],p2p[k,2:4])<6){
          xx = rbind(xx,c(p1p[j,1],p2p[k,1]))
        }
      }
    }
  }
  write.table(xx, file = paste0(p1,"_",p2,"_contact6.table"))
}

#sss("1pg5_A.pdb", "1pg5_B.pdb")

#sss("1DE4_r_u_A.pdb","1DE4_l_u_A.pdb")
#sss("1H1V_r_u_B.pdb","1H1V_l_u_B.pdb")

#sss("1KXP_r_u_B.pdb","1KXP_l_u_B.pdb")

