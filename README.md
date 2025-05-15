<img src="head.png" alt="My Logo" width="800" />


# Spanning Motif Identifier for RNA–Protein Complexes

A Python-based tool to identify and classify **spanning motifs** in RNA–protein complexes using hydrogen bond data and corresponding PDB structures.

<!--- BADGES: START --->
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
<!--- BADGES: END --->

## Overview

Spanning triplets comprise a single amino acid hydrogen bonded with two ribonucleotides (not mutually involved in a Watson–Crick base pair), or vice versa. With their high specificity, these motifs can aid in a better understanding of RNA–protein recognition and associated phenomena. To achieve this, we developed a Python algorithm that provides a topology-based classification and nomenclature scheme to automatically identify spanning triplet motifs in the crystal structures of RNA–protein complexes. This technique can also be extended to include higher-order spanning motifs.

This pipeline processes `.hb2` files (from HBPLUS) along with PDB files.

### Scripts:

1. `1_aa_motifs_False.py`
2. `2_nt_motifs_False.py`
3. `3_Remove_false_spanning_motifs.py`
4. `4_topologies_amino_acid_centered_spanning_motifs.py`
5. `5_topologies_nucleotide_centered_spanning_motifs.py`

## Requirements

- Python 3.6
- `.hb2` files generated using **HBPLUS**
- Corresponding PDB files

## Usage
Place the scripts in a single folder along with:
- `.hb2` files
- Corresponding `.pdb` files


**Run**   
  -  python3 `1_aa_motifs_False.py`
  -  python3 `2_nt_motifs_False.py`  
     (These scripts will generate spanning motifs and some extended network of hydrogen-bonded amino acids and nucleotides (if found). The name of their output file will 
      be ‘aa_motifs_False.txt’ and ‘nt_motifs_False.txt’ respectively.)     
  -  python3 `3_Remove_fasle_spanning_motifs.py`     
     (It will separate spinning motifs and other extended networks of amino acids and nucleotides. Here, a maximum of four output files will be generated: ‘AA_spanningmotifs.txt’, ‘NT_spanningmotifs.txt’, ‘Extended_networks_1.txt’, and ‘Extended_networks_2.txt’.)
-   python3 `4_topologies- amino acid centered spanning motifs.py`   
-   python3 `5_topologies- nucleotide centered spanning  motifs.py`  
    (These scripts classify the previously identified spanning motifs and generate a separate output file for each topology.)

## Dataset


## Contact

For queries or support, please feel free to contact us:

- **Dr. Purshotam Sharma**: [psharmapuchd@gmail.com](mailto:psharmapuchd@gmail.com)

    


