#!/usr/bin/Rscript
library(bio3d)
source("step1.R")
source("step2.R")
source("step3.R")
source("step4.R")
source("12contact_chain_v2.r")

#names:

#input names:
#xxxx_l_u_X.pdbX.fas
#xxxx_r_u_X.pdbX.fas
queue <- read.table("queue.table")
for(i in 1:nrow(queue)){
  l <- queue[i,1]
  r <- queue[i,2]
  
  l_pdb <- substr(l,1,14)
  r_pdb <- substr(r,1,14)
  # PYRI_4fyyB.alignment
  # PYRI_4fyyB.out
  # PYRI_fseq_4fyy_BD_B.pdbB.fas
  # 
  # 1A2K_l.alignment
  # 1A2K_l.alignment.fasta
  # 
  # PYRB_1ekxA.alignment.fasta
  # PYRB_1ekxA_w_1vbgA.link.fasta
  # 
  # 1A2K_lArB.aln
  # 1A2K_lArB.aln
  # 1A2K_lArB.mat
  # 
  # 1A2K_lArB.mat
  
  
  #pre-step
  #the alignments and output file for align back to query index for later steps from HMMer - phmmer
  # & translate the HMMer stockholm format back to Fasta
  system(paste0("/home/lidaphd/Desktop/hmmer-3.1b2/src/phmmer --notextw --mxfile /home/lidaphd/Desktop/Desktop/test.mat -E 0.1 -A /home/lidaphd/data/autofrag/", l, ".alignment -o /home/lidaphd/data/autofrag/", l, ".out /home/lidaphd/data/b4chainnew/", l, " /home/lidaphd/Desktop/Desktop/library.fas.fasta"), wait = TRUE)
  system(paste0("/home/lidaphd/Desktop/hmmer-3.1b2/src/phmmer --notextw --mxfile /home/lidaphd/Desktop/Desktop/test.mat -E 0.1 -A /home/lidaphd/data/autofrag/", r, ".alignment -o /home/lidaphd/data/autofrag/", r, ".out /home/lidaphd/data/b4chainnew/", r, " /home/lidaphd/Desktop/Desktop/library.fas.fasta"), wait = TRUE)
  
  # system("./phmmer --mxfile /home/lidaphd/lidaphd/Desktop/test.mat -E 0.1 -A /home/lidaphd/lidaphd/Gremlin_test/Gremlin_test/PYRI_4fyyB.alignment --tblout /home/lidaphd/lidaphd/Gremlin_test/Gremlin_test/PYRI_4fyyB.tblout -o /home/lidaphd/lidaphd/Gremlin_test/Gremlin_test/PYRI_4fyyB.out /home/lidaphd/lidaphd/Gremlin_test/Gremlin_test/PYRI_fseq_4fyy_BD_B.pdbB.fas /home/lidaphd/lidaphd/Desktop/library.fas.fasta", wait = TRUE)
  
  system(paste0("grep -A 5 '== domain 1' ", l,".out | sed -n '2p;4p' > tt.txt"), wait = T)
  system(paste0("awk 'BEGIN {ORS = \"\"} {print \">\"} {print $2} {print \"-\"} {print $4} {printf \"\n\"} {print $3} {printf \"\n\"}' tt.txt > ", l, ".link.fasta", wait = T))
  system(paste0("grep -A 5 '== domain 1' ", r,".out | sed -n '2p;4p' > tt.txt"), wait = T)
  system(paste0("awk 'BEGIN {ORS = \"\"} {print \">\"} {print $2} {print \"-\"} {print $4} {printf \"\n\"} {print $3} {printf \"\n\"}' tt.txt > ", r, ".link.fasta", wait = T))
  
  system(paste0("/home/lidaphd/Desktop/HHsuite/hhsuite-2.0.16-linux-x86_64/scripts/reformat.pl sto fas /home/lidaphd/data/autofrag/", l, ".alignment /home/lidaphd/data/autofrag/", l, ".alignment.fasta"), wait = T)
  system(paste0("/home/lidaphd/Desktop/HHsuite/hhsuite-2.0.16-linux-x86_64/scripts/reformat.pl sto fas /home/lidaphd/data/autofrag/", r, ".alignment /home/lidaphd/data/autofrag/", r, ".alignment.fasta"), wait = T)
  
  
  # system("./reformat.pl sto fas /home/lidaphd/lidaphd/Gremlin_test/Gremlin_test/1A2K_l.alignment /home/lidaphd/lidaphd/Gremlin_test/Gremlin_test/1A2K_l.alignment.fasta", wait = T)
  
  #step1
  #protein1
  alignment_interface_gap_l <- step1(paste0(l, ".alignment.fasta"))
  alignment_l <- alignment_interface_gap_l[[1]]
  interface_l <- alignment_interface_gap_l[[2]]
  gap_l <- alignment_interface_gap_l[[3]]
  #protein2
  alignment_interface_gap_r <- step1(paste0(r, ".alignment.fasta"))
  alignment_r <- alignment_interface_gap_r[[1]]
  interface_r <- alignment_interface_gap_r[[2]]
  gap_r <- alignment_interface_gap_r[[3]]
  
  
  #step2
  #p1
  queryfile <- file(paste0(l, ".link.fasta"))
  firstline <- readLines(queryfile, n=1)
  splat <- strsplit(firstline,">")[[1]]
  splat2 <- splat[[2]]
  querystart <- strsplit(splat2,"-")[[1]]
  querystart2 <- querystart[[1]]
  queryend <- querystart[[2]]
  
  #l_length <- as.numeric(queryend) - as.numeric(querystart2)
  close(queryfile)
  queryindex_l <- step2(0, as.numeric(querystart2), paste0(l, ".link.fasta"), alignment_l)
  #p2
  queryfile <- file(paste0(r, ".link.fasta"))
  firstline <- readLines(queryfile, n=1)
  splat <- strsplit(firstline,">")[[1]]
  splat2 <- splat[[2]]
  querystart <- strsplit(splat2,"-")[[1]]
  querystart2 <- querystart[[1]]
  close(queryfile)
  queryindex_r <- step2(0, as.numeric(querystart2), paste0(r, ".link.fasta"), alignment_r)
  
  #step3
  pair_msa_interface_gap <- step3(alignment_l, alignment_r, interface_l, interface_r, gap_l, gap_r, queryindex_l, queryindex_r)
  queryindex <- pair_msa_interface_gap[[4]]
  alignment_pair <- pair_msa_interface_gap[[1]]
  
  l_aftercut <- step3(alignment_l, alignment_l, interface_l, interface_l, gap_l, gap_l, queryindex_l, queryindex_l)
  l_alignment <- l_aftercut[[1]]
  l_length <- (ncol(l_alignment) / 2)
  ##inter-step
  #output the paired matrix for CCMpred & CCMpred predict & transfer the CCMpred output back to our prediction tables
  write.table(alignmatrix_pair, file = paste0(l,'_',r,'.aln'), quote = F, sep = '', row.names = F, col.names = F)
  system(paste0("/home/lidaphd/Desktop/CCMpred/bin/ccmpred ", l, "_", r, ".aln ",l, "_", r, ".mat"), wait = T)
  
  system(paste0("/home/lidaphd/Desktop/CCMpred/scripts/top_couplings.py -s ", l_length, " -n 1000 ", l, "_", r, ".mat > pred_",l,"_",r,".table"), wait = T)
  
  #step4
  #print TP FP Presicion Recall etc.
  #prediction_table is the prediction
  #truecontact_table is the true contacts
  prediction_table <- paste0("pred_",l,"_",r,".table")
  
  sss(l_pdb, r_pdb)
  truecontact_table <- paste0(l,"_",r,"_contact6.table")
  
  
  l1 <- length(read.fasta(l)$ali)
  l2 <- length(read.fasta(r)$ali)
    
  step4(prediction_table, truecontact_table, l1, l2, queryindex)
  
}
