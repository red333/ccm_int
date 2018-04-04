#table from true contact, remove the first row "0 0"
#b = b[-1,]

# minor the i seq L from j column
#a$j <- a$j - 311
grem_table <- read.table("1a2k_grem.txt",header = T)
grem <- grem_table[grem_table$gene == "AB",]

predict_table <- read.table("pred_1A2K_lArB.20.table")
truecontact <- read.table("1A2K_l_u_A.pdb_1A2K_r_u_B.pdb_contact6.table")
truecontact <- unique(truecontact[-1,])
a <- predict_table[,c(2,1)] + 1
for(i in 1:nrow(a)){
  # a = a[-1,]
  #print(a[i,1])
  a[i,1] <- queryindex[2,as.numeric(a[i,1])]
  a[i,2] <- queryindex[2,as.numeric(a[i,2])]
}
#desend order the table for top L/K
#a[order(a$V7, decreasing = T),]
abs <- 2

b <- truecontact


#PYRBPYRL 311 153

l1 = 202
l2 = 121
L = l1 + l2
K = c(10,5,2,1)
#could be L/K later
lk = round(L/K)
count = 0
index = c()
TP = c()
FP = c()
FN = c()
for(k in 1:4){
  for(i in 1:lk[k]){
    for(j in 1:nrow(b)){
      #+/- 2 for more TP
      if(as.numeric(as.character(a$V2[i])) <= as.numeric(as.character(b$V1[j])) + abs & as.numeric(as.character(a$V2[i])) >= as.numeric(as.character(b$V1[j])) - abs){
        if(as.numeric(as.character(a$V1[i])) <= as.numeric(as.character(b$V2[j])) + abs & as.numeric(as.character(a$V1[i])) >= as.numeric(as.character(b$V2[j])) - abs){
          #count = count + 1
          index = append(index, i)
        }
      }
    }
  }
  count2 = unique(index)
  count = length(count2)
  
  print(count)
  print(count2)
  TP = c(TP,count)
  FP = c(FP,lk[k]-count)
  FN = c(FN, nrow(b) - count)
  P = TP / (TP + FP)
  R = TP / (TP + FN)
  print(paste0("Number of TP: ", TP))
  print(paste0("Number of FP: ", FP))
  print(paste0("Precsion = ", P))
  print(paste0("Recall = ", R))
}