library(bio3d)
interfacelib = read.table("vi.table")

a = read.fasta("PYRB_1ekxA.alignment.fasta",rm.dup = T, to.upper = T, to.dash = T)
id = a$id
ali = a$ali
alignmatrix = as.data.frame(ali, stringsAsFactors = F)

gapbinary <- alignmatrix
gapbinary[ gapbinary != "-"] <- "1"
gapbinary[ gapbinary == "-"] <- "0"

interfacebinary <- alignmatrix


for(j in 1:nrow(interfacebinary)) {
  tryCatch ({
  x <- rownames(interfacebinary[j,])
  protein <- substr(x,1,5)
  chain <- substr(x,6,6)
  index <- gsub(".*/","",x)
  startindex <- gsub("-.*","",index)
  endindex <- gsub("*.-","",index)
  
  #interfacepair <- "TBA"
  interfacepair <- interfacelib[which(substr(interfacelib[,1],1,6) == paste0(protein,chain)),1]
  if(length(interfacepair)!=0){
    
    interfacetable <- read.table(paste0("/home/lidaphd/lidaphd/InterEvolLib/superpdbs/out_" ,interfacepair[1], ".pdb.interfaceinfo/molecule_1.txt"))
    #interfacetable <- read.table("molecule_1.txt")
    interfacetable <- interfacetable[,2]
    interfacetable <- unique(interfacetable)
    
    indextrack = as.numeric(startindex)
    for(i in 1:ncol(interfacebinary)){
      if(interfacebinary[j,i] != "-"){
        if(indextrack %in% interfacetable){
          interfacebinary[j,i] <- "1"
        }else{
          interfacebinary[j,i] <- "0"
        }
        
        indextrack <- indextrack + 1
      }
    }
    
  }
  }, error = function(e) {})
}

for(j in 1:nrow(interfacebinary)){
  tryCatch({
  x <- rownames(interfacebinary[j,])
  protein <- substr(x,1,5)
  chain <- substr(x,6,6)
  index <- gsub(".*/","",x)
  startindex <- gsub("-.*","",index)
  endindex <- gsub("*.-","",index)
  
  #interfacepair <- "TBA"
  interfacepair <- interfacelib[which(paste0(substr(interfacelib[,1],1,5),substr(interfacelib[,1],7,7)) == paste0(protein,chain)),1]
  if(length(interfacepair)!=0){
    
    interfacetable <- read.table(paste0("/home/lidaphd/lidaphd/InterEvolLib/superpdbs/out_" ,interfacepair[1], ".pdb.interfaceinfo/molecule_2.txt"))
    #interfacetable <- read.table("molecule_1.txt")
    interfacetable <- interfacetable[,2]
    interfacetable <- unique(interfacetable)
    
    indextrack = as.numeric(startindex)
    for(i in 1:ncol(interfacebinary)){
      if(interfacebinary[j,i] != "-"){
        if(indextrack %in% interfacetable | interfacebinary[j,i] == "1"){
          interfacebinary[j,i] <- "1"
        }else{
          interfacebinary[j,i] <- "0"
        }
        
        indextrack <- indextrack + 1
      }
    }
    
  }
  }, error = function(e) {})
}

b = read.fasta("PYRI_4fyyB.alignment.fasta",rm.dup = T, to.upper = T, to.dash = T)
id_2 = b$id
ali_2 = b$ali
alignmatrix_2 = as.data.frame(ali_2, stringsAsFactors = F)

gapbinary_2 <- alignmatrix_2
gapbinary_2[ gapbinary_2 != "-"] <- "1"
gapbinary_2[ gapbinary_2 == "-"] <- "0"

interfacebinary_2 <- alignmatrix_2


for(j in 1:nrow(interfacebinary_2)) {
  tryCatch ({
    x <- rownames(interfacebinary_2[j,])
    protein <- substr(x,1,5)
    chain <- substr(x,6,6)
    index <- gsub(".*/","",x)
    startindex <- gsub("-.*","",index)
    endindex <- gsub("*.-","",index)
    
    #interfacepair <- "TBA"
    interfacepair_2 <- interfacelib[which(substr(interfacelib[,1],1,6) == paste0(protein,chain)),1]
    if(length(interfacepair_2)!=0){
      
      interfacetable_2 <- read.table(paste0("/home/lidaphd/lidaphd/InterEvolLib/superpdbs/out_" ,interfacepair_2[1], ".pdb.interfaceinfo/molecule_1.txt"))
      #interfacetable <- read.table("molecule_1.txt")
      interfacetable_2 <- interfacetable_2[,2]
      interfacetable_2 <- unique(interfacetable_2)
      
      indextrack = as.numeric(startindex)
      for(i in 1:ncol(interfacebinary_2)){
        if(interfacebinary_2[j,i] != "-"){
          if(indextrack %in% interfacetable_2){
            interfacebinary_2[j,i] <- "1"
          }else{
            interfacebinary_2[j,i] <- "0"
          }
          
          indextrack <- indextrack + 1
        }
      }
      
    }
  }, error = function(e) {})
}

for(j in 1:nrow(interfacebinary_2)){
  tryCatch({
    x <- rownames(interfacebinary_2[j,])
    protein <- substr(x,1,5)
    chain <- substr(x,6,6)
    index <- gsub(".*/","",x)
    startindex <- gsub("-.*","",index)
    endindex <- gsub("*.-","",index)
    
    #interfacepair <- "TBA"
    interfacepair_2 <- interfacelib[which(paste0(substr(interfacelib[,1],1,5),substr(interfacelib[,1],7,7)) == paste0(protein,chain)),1]
    if(length(interfacepair_2)!=0){
      
      interfacetable_2 <- read.table(paste0("/home/lidaphd/lidaphd/InterEvolLib/superpdbs/out_" ,interfacepair[1], ".pdb.interfaceinfo/molecule_2.txt"))
      #interfacetable <- read.table("molecule_1.txt")
      interfacetable_2 <- interfacetable_2[,2]
      interfacetable_2 <- unique(interfacetable_2)
      
      indextrack = as.numeric(startindex)
      for(i in 1:ncol(interfacebinary_2)){
        if(interfacebinary_2[j,i] != "-"){
          if(indextrack %in% interfacetable_2 | interfacebinary_2[j,i] == "1"){
            interfacebinary_2[j,i] <- "1"
          }else{
            interfacebinary_2[j,i] <- "0"
          }
          
          indextrack <- indextrack + 1
        }
      }
      
    }
  }, error = function(e) {})
}







