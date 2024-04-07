###########################################################################################################################
#                                          Identification of spanning motifs - 3                                          #
###########################################################################################################################
# The previous two scripts search for amino acid- and nucleotide- mediated motifs. But there are many cases where the same amino
# acids or nucleotides appear in both amino acid- and nucleotide- mediated motifs. That means these a not spanning motifs
# but are part of some extended network of hydrogen-bonded amino acids and nucleotides. Our study is focused only on
# spanning motifs, so we have to remove such large networks which can be analyzed later.

# Importing the required modules
import os

# A. To separate false amino acid-mediated spanning motifs.
output1 = open('aa_motifs_False.txt', 'r') # Reads output of script "1_aa_motifs_False.py".
data = output1.readlines()
dict_of_pdbs_A1N2 = {}
start = []
end = []
s = -1
e = -1
for d in data:
    s += 1
    e += 1
    if '****************' in d: start.append(s)
    if '......................................................' in d: end.append(e)


for i in range(len(start)):
    pdb_name = data[start[i]].split('********************  ')[1].split('  ********************')[0]
    dict_of_pdbs_A1N2[pdb_name] = data[start[i]:end[i]]
    
    
output2 = open('nt_motifs_False.txt', 'r') # Reads output of script "2_nt_motifs_False.py".
data = output2.readlines()
dict_of_pdbs_A2N1 = {}
start = []
end = []
s = -1
e = -1
for d in data:
    s += 1
    e += 1
    if '****************' in d: start.append(s)
    if '......................................................' in d: end.append(e)


for i in range(len(start)):
    pdb_name = data[start[i]].split('********************  ')[1].split('  ********************')[0]
    dict_of_pdbs_A2N1[pdb_name] = data[start[i]:end[i]]

    
    
dict_removed_A1N2 = {} 
for p in dict_of_pdbs_A2N1:
    pdb_name = '********************  ' + p + '  ********************'
    val = []
    if p in dict_of_pdbs_A1N2:
        aa_nt = []
        for j in dict_of_pdbs_A2N1[p]:
            if "PDB:" not in j:
                splt = j.strip().split(' --- ')
                for s in splt:
                    if (s.strip()[0:9]) not in aa_nt:
                        aa_nt.append(s.strip()[0:9])
                    if (s.strip()[14:23]) not in aa_nt: 
                        aa_nt.append(s.strip()[14:23])
            
            for an in aa_nt:
                for i in dict_of_pdbs_A1N2[p]:
                    if an in i:
                        if i in dict_of_pdbs_A1N2[p]:
                            val.append(i)
                            dict_of_pdbs_A1N2[p].remove(i)
                            
    dict_removed_A1N2[pdb_name] = val
                            
new_dict1 = {key: val for key, val in dict_of_pdbs_A1N2.items() if len(val) != 1}
for i in new_dict1:
    for j in new_dict1[i]:
        print(j.strip(), file=open("AA_spanningmotifs.txt", "a"))
    print('...............................................', file=open("AA_spanningmotifs.txt", "a"))    
    
    

new_dict2 = {key: val for key, val in dict_removed_A1N2.items() if len(val) > 0}
for i in new_dict2:
    print(i.strip(), file=open("Extended_networks_1.txt", "a"))
    for j in new_dict2[i]:
        print(j.strip(), file=open("Extended_networks_1.txt", "a"))
    print('...............................................', file=open("Extended_networks_1.txt", "a"))                           

    
######################################################################################################
# B. To separate false nucleotide-mediated spanning motifs.

output1 = open('aa_motifs_False.txt', 'r') # Reads output of script "1_aa_motifs_False.py".
data = output1.readlines()
dict_of_pdbs_A1N2 = {}
start = []
end = []
s = -1
e = -1
for d in data:
    s += 1
    e += 1
    if '****************' in d: start.append(s)
    if '......................................................' in d: end.append(e)


for i in range(len(start)):
    pdb_name = data[start[i]].split('********************  ')[1].split('  ********************')[0]
    dict_of_pdbs_A1N2[pdb_name] = data[start[i]:end[i]]
        
output2 = open('nt_motifs_False.txt', 'r') # Reads output of script "2_nt_motifs_False.py".
data = output2.readlines()
dict_of_pdbs_A2N1 = {}
start = []
end = []
s = -1
e = -1
for d in data:
    s += 1
    e += 1
    if '****************' in d: start.append(s)
    if '......................................................' in d: end.append(e)


for i in range(len(start)):
    pdb_name = data[start[i]].split('********************  ')[1].split('  ********************')[0]
    dict_of_pdbs_A2N1[pdb_name] = data[start[i]:end[i]]

dict_removed_A2N1 = {} 
for p in dict_of_pdbs_A1N2:
    pdb_name = '********************  ' + p + '  ********************'
    val = []
    if p in dict_of_pdbs_A2N1:
        aa_nt = []
        for j in dict_of_pdbs_A1N2[p]:
            if "PDB:" not in j:
                splt = j.strip().split(' --- ')
                for s in splt:
                    if (s.strip()[0:9]) not in aa_nt:
                        aa_nt.append(s.strip()[0:9])
                    if (s.strip()[14:23]) not in aa_nt: 
                        aa_nt.append(s.strip()[14:23])

        for an in aa_nt:
            for i in dict_of_pdbs_A2N1[p]:
                if an in i:
                    if i in dict_of_pdbs_A2N1[p]:
                        val.append(i)
                        dict_of_pdbs_A2N1[p].remove(i)
                        
    dict_removed_A2N1[pdb_name] = val
    
new_dict2 = {key: val for key, val in dict_of_pdbs_A2N1.items() if len(val) != 1}
for i in new_dict2:
    for j in new_dict2[i]:
        print(j.strip(), file=open("NT_spanningmotifs.txt", "a"))
        
    print('...............................................', file=open("NT_spanningmotifs.txt", "a"))
    
    
new_dict3 = {key: val for key, val in dict_removed_A2N1.items() if len(val) > 0}
for i in new_dict3:
    print(i.strip(), file=open("Extended_networks_2.txt", "a"))
    for j in new_dict3[i]:
        print(j.strip(), file=open("Extended_networks_2.txt", "a"))
    print('...............................................', file=open("Extended_networks_2.txt", "a")) 
             

