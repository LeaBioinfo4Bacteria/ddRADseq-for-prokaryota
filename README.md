# ddRADseq-for-prokaryota
This project was written as a first step to apply the ddRADSeq genotyping method to prokaryota.
This code is made to dertemine, in silico, the couple of enzymes to apply on a given sequence. It is based on the work of Charlotte COUCHOUD (University of Besançon, France) who started the in vitro applications. 

Here, we aim to select the best enzymes, either from the ones we already had (NEB supply) or other databases, on the Pseudomonas aeruginosa PAO1 genome and pGFP and pMMB207 plasmids as controls. As such, the loop_souchier.sh file parse the files contained in the seq_test directory. 

If you have several sequences, you can :
- change the name and/or path of the repository seq_test in loop_souchier.sh to match yours
- put your sequences in seq_test

If you have one sequence, you can launch main_multi.py on it. 

Here are the arguments you can use : 

  Argument		Description				                        Default value
--ratio, -r		ratio of cutting frequencies of enzyme pairs, must be > 1	3
--fragment, -f 		minimum number of interest fragments				0
--lenmin, -lmin		minimal size in base pair of interest fragments			200
--lenmax, -lmax		maximal size in base pair of interest fragments			800
--linear, -lin		linear genome					False
--tab_supp, -t 		table with additional information on enzymes tested		nebuffer-performancechart-with-restrictionenzymes.csv


In the repository, you will find the workflow chart and an overview of the created files and interactions.

