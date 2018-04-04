for(i in 2:nrow(a)){
  # a = a[-1,]
  print(a[i,1])
  a[i,1] <- queryindex[2,as.numeric(a[i,1]) + 1]
  a[i,2] <- queryindex[2,as.numeric(a[i,2]) + 1]
}
