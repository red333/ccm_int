README FOR RELATED R SCRIPT FILES

1.s2f.R
This is the file we generate 1-d seqs for structure fragments data based on PDB files.
The length of the fragment seqs is N-4, where N is the length of orignal fasta AA seqs.
When the resisdue index is missed, the position we use a L for it. THe 11th type from BLOSUM62. for structure gaps



The structural fragments library for protein chains we choose, is the 10-5 fragments library from Kolodny[1]. 10-5 means we have overall ten fragments type and each fragment is decided by five Calpha atoms in the protein chain. The reason we choose the 10-5 version from the other library with different total types and size of each fragment is that it is the most accurate library after the improvement by Kolodny's research.
The original library from Kolodny's database marks the ten types as 0-9. In our case, we want to use the alphabets fragments data so we convert the 0-9 to A-J.
For the protein-protein interfaces, we use the interfaces library from InterEvol database, which is designed for exploring 3D structures of homologous interfaces of protein complexes. The version we choose from InterEvol database is the INTER70 table which lists all the non-redundant reference interfaces (first column) and those considered as redundant. Every line corresponds to a group of redundant interfaces. Interfaces are defined as dimers interacting through at least 10 contacts involving different amino acids. Contacts were defined using a 5 A full atom cut-off. For 2 dimers AB and A'B', if AA' and BB' share at least 70% sequence identity and 70% of coverage, and the position of residues at the interfaces AA' and BB' overlay by more than 40%, AB and A'B' are defined as redundant interface. The PDB in the first column generally corresponds to the structure of best resolution (cluster reference).
What we do is that we cut the first column and paste them into another table file and let the shell script download all PDB files for those protein IDs from RCSB.
Finally, we use R script to finish the process of convert interfaces library into a structural fragments library. Inside the R script, we use the Bio3d package, which provides the PDB file read function for us. And to calculate the coordinate root-mean-square deviation of the C?? atom to measure the structural similarity of any two fragments.
We set the sliding window size as a parameter in our R script and it is currently set to 2 in our case (the sliding window size means how many atoms you want to jump across each time you start a new 5 atoms fragment inside a chain). Because we are using a 10-5 fragment library and we want to improve our efficiency by decrease the data input size, which is an advantage of structural fragments compare to 3D structure data.

2.whole steps
3.step1
4.step2
5.step3
6.step4

wholestep.R & step 1,2,3,4.R
The wholestep is the overall automatic script to take the list of jobs and run them.
Step 1 takes the MSA from jackhmmer function from hhsuite, and our interface residue library. The output from step1 is the MSA with interface residue percentage and gap percentage cutoff.
Step 2 takes the output of step1 and take the output of alignment.o files for residue index assignment for query seqs.
Step 3 takes the percentage informations and do the deletion based on percentages. And output the final MSA to CCMpred and save the residue index for evaluation. 
Step 4 check the top xxx predictions and output the TP counts and other prediction stats.

Because we generate the MSA for AA and Frags separately, we face the residue index where we may not have the index available in both AA MSA and Frags MSA. In this case, we put a -1 for the index matrix and # for the whole column in MSA matrix for the MSA where the residue is unavailable.  
