###########################################################################################################################
#                                          Identification of spanning motifs - 2                                          #
###########################################################################################################################

# Importing the required modules
import os, re, itertools

t_c = 0
motif = {}

# A used defined function to remove ambiguity arising due to hydrogen bond donor as well as acceptor nature of 
# some atoms as HBPLUS prints the same interaction twice by exchanging the role (hydrogen bond donor as well as 
# acceptor) of interacting atoms.def uq_interactions(set1):
def uq_interactions(set1):
    extra = []
    for h in range(len(set1)):
        for i in range(h+1,len(set1)):
            if (set1[h][:27] == set1[i][:27]):
                extra.append(set1[h])
    for h in extra:
        if h in set1:
            set1.remove(h)

            
# List of nucleotides and amino acids.
nucleotides = ['-  A', '-  C', '-  G', '-  U' ]
all_amino_acids = ['-GLY', '-ALA', '-SER', '-THR', '-CYS', '-VAL', '-LEU', '-ILE', '-MET', '-PRO', '-PHE', '-TYR', '-TRP', '-ASP', '-GLU', '-ASN', '-GLN', '-HIS', '-LYS', '-ARG']

# Find and read output files of HBPLUS (.hb2) present in the current directory.
files_pdbs = []
path = os.getcwd()
directory = os.listdir(path)  # specify folder path.
for fis in directory:             # fis is files in directory.
    if fis.find('.hb2') >= 0:
        files_pdbs.append(fis)
files_pdbs.sort()

# Scripting each .hb2 file recursively to identify water bridges.
for f in files_pdbs:
    ff = 'PDB: ' + f[3:7].upper()
    motif_list = []
    t_c += 1
    head = '********************  ' + ff + '  ********************'
    print(head)
    #print(t_c, '  ', f )
    file = open(f, 'r')
    data = file.readlines()
    amino_acids_int_nb = [] # List containing hydrogen bonds of amino acid and nucloeotides.
    # Rearranging the format of writing hydrogen bonds of nucleotides with amino acids; 
    # information of nucleotide will be written first followed by amino acids. 
    for l in range(len(data)):
            if data[l][5:9] in nucleotides and data[l][19:23] in all_amino_acids:
                r = data[l][:33]
                amino_acids_int_nb.append(r.strip())
            if data[l][5:9] in all_amino_acids and data[l][19:23] in nucleotides:
                p1 = data[l][:14]
                p2 = data[l][14:28] 
                p3 = data[l][28:33]
                r = p2 + p1 + p3
                amino_acids_int_nb.append(r.strip())
    uq_interactions(amino_acids_int_nb)  # This will give unique interactions, i.e. will remove repetitions due to 
                                         # to ambigous donors and acceptors.
    aa_resi_uq = [] # List of nucleotides intercting with amino acids along with their residue numbers.
    for a in amino_acids_int_nb:
        if a[:9].strip() not in aa_resi_uq:
            aa_resi_uq.append(a[:9].strip())
    aa_resi_uq.sort()
    amino_acids_int_nb.sort()
    
    # Grouping all hydrogen bonds of a particular nucleotide, resulting in a spanning motif formed by nucleotide 
    # and amino acids.
    strings = amino_acids_int_nb
    result = []
    iterator = itertools.groupby(strings, lambda string: string[:9])
    
    
    
    # Here 'element' is information of nucleotide, and 'group' is all interaction of that nucleotide. 
    # Making a list by grouping all interacting amino acids with particualar nucleotide.
    for element, group in iterator:
       # appending the group by converting it into a list
       result.append(list(group))
    joined_result = []
    for i in result:
        joined_result.append(" --- ".join(i))
        
    triplets = []
    
    for i in joined_result:
        res = len(re.findall('GLY|ALA|SER|THR|CYS|VAL|LEU|ILE|MET|PRO|PHE|TYR|TRP|ASP|GLU|ASN|GLN|HIS|LYS|ARG', i))
        if res > 1 : triplets.append(i)
            
    ## Remvoe those false spanning motifs, which include only one amino acid, i.e. nucleotide interacts only 
    # with single amino acid, but multiple times.
    to_be_removed = []
    for t in triplets:
        if len(t) == 69 :
            if t[51:60] == t[14:23]: to_be_removed.append(t)
        if len(t) == 106  :
            if t[14:23] == t[51:60] == t[88:97] : to_be_removed.append(t)
        if len(t) == 143  :
            if t[51:60] == t[14:23] == t[88:97] == t[125:134] : to_be_removed.append(t)
        if len(t) == 180  :
            if t[51:60] == t[14:23] == t[88:97] == t[125:134] == t[162:171] : to_be_removed.append(t)
        if len(t) == 217  :
            if t[51:60] == t[14:23] == t[88:97] == t[125:134] == t[162:171] == t[199:208]: to_be_removed.append(t)
        if len(t) == 254  :
            if t[51:60] == t[14:23] == t[88:97] == t[125:134] == t[162:171] == t[199:208] == t[236:245]: to_be_removed.append(t)
    for z in to_be_removed: triplets.remove(z)
  
    
    # Searcing amino acid - amino acid interactions ... 
    # i) create a list containing all amino acid - amino acid  hydrogen bonds. 
    #    There repition of interactions due to atoms capable of acting as acceptor as well donors.
    # ii) Rearange r1 = p1 + p2 + p3 to r2 = p2 + p1 + p3. 
    # iii) append all r1 to AA_AA, and if r2 is in AA_AA, aAA_AA it to another list and the remove this extra list fron AA_AA.
    AA_AA = [] # List of amino acid - amino acid  hydrogen bonds. 
    repetitive_AA = []
    for l in range(len(data)):
        if data[l][5:9] in all_amino_acids and data[l][19:23] in all_amino_acids: 
            AA_AA.append(data[l][:33].strip())

    for l in range(len(data)):
        p1 = data[l][:14]
        p2 = data[l][14:28] 
        p3 = data[l][28:33]
        r = p2 + p1 + p3
        if r.strip() in AA_AA: repetitive_AA.append(r.strip())  
    repetitive_AA_1 = []       
    for l in range(1,len(repetitive_AA),2):
        repetitive_AA_1.append(repetitive_AA[l])
    for l in repetitive_AA_1:
        if l in AA_AA:AA_AA.remove(l)
            
    # Adding more information to the spanning motifs; i.e., whether cyclic or noncyclic.
    # For this, hydrogen bonds of amino acids involved in formation of spanning motifs will be searched            
    my_dict = {}
    key = []
    val = []
    for i in triplets:
        nb_nb = []
        for a in range(14,len(i),37):
            if i[a:a+9] not in nb_nb:
                nb_nb.append(i[a:a+9])
                
        for n in range(len(nb_nb)):
            for m in range(n+1, len(nb_nb)):
                for l in (data): 
                    if nb_nb[m] in l and nb_nb[n] in l:
                        key.append(i)
                        val.append(l[:33])   
        if i not in key:
            key.append(i)
            val.append('NIL')
        
    for i in range(len(key)):
        for j in range(len(val)):
            if i == j: 
                my_dict.setdefault(key[i], []).append(val[j])
    
    for a, b in my_dict.items():
        x = str(a) + ' ' + str(b)
        motif_list.append(x)
    
    motif_list.append('........................................................................')
    motif[head] = motif_list

new_dict = {key: val for key, val in motif.items() if len(val) != 1}


# Printing results.
for i in new_dict:
    print(i, file=open("nt_motifs_False.txt", "a"))
    for j in new_dict[i]:
        print(j, file=open("nt_motifs_False.txt", "a"))