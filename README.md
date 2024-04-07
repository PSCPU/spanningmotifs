A python code for identification of spanning motifs in RNA–protein complexes 
This code is written in python3. To run this code, go through the following steps:
1.	Place all five codes (‘1_aa_motifs_False.py’, ‘2_nt_motifs_False.py’, ‘3_Remove_fasle_spanning_motifs’, ‘4_topologies- amino acid centered spanning motifs’, and ‘5_topologies- nucleotide centered spanning 
   motifs’) into a folder along with ‘.hb2’ files obtained from HBPLUS and respective PDB files.
2.	Run ‘1_aa_motifs_False.py’,  and ‘2_nt_motifs_False.py’ these script will generate spanning motifs along with some extended network of hydrogen bonded amino acids and nucleotides. Name of their output file will 
    be ‘aa_motifs_False.txt’ and ‘nt_motifs_False.txt’ respectively.
3.	To remove false entries, run script ‘3_Remove_fasle_spanning_motifs’. It will separate spinning motifs and other extended network of amino acids and nucleotides. Here, maximum four output file will be 
    generated:‘AA_spanningmotifs.txt’, ‘NT_spanningmotifs.txt’, ‘Extended_networks_1.txt’, and ‘Extended_networks_2.txt’.
4.	Further identified spanning motifs can be classified as per their topology using ‘4_topologies- amino acid centered spanning motifs.py’ and ‘5_topologies- nucleotide centered spanning  motifs’. A separate file 
    for each topology will be generated. 
