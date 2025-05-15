# Spanning Motif Identifier for RNA–Protein Complexes

A Python-based pipeline to identify and classify **spanning motifs** in RNA–protein complexes using hydrogen bond data.

## Overview

This toolkit processes hydrogen bond data (`.hb2` files from HBPLUS) and PDB structures to extract amino acid and nucleotide-centered spanning motifs, filter out false positives, and classify the valid motifs by topology.

## Requirements

- Python 3.x
- `.hb2` files generated using [HBPLUS](http://www.ebi.ac.uk/thornton-srv/software/HBPLUS/)
- Corresponding PDB files

## File Structure

Place the following scripts in a single folder along with:
- `.hb2` files
- Corresponding `.pdb` files

### Scripts:

1. `1_aa_motifs_False.py`
2. `2_nt_motifs_False.py`
3. `3_Remove_false_spanning_motifs.py`
4. `4_topologies_amino_acid_centered_spanning_motifs.py`
5. `5_topologies_nucleotide_centered_spanning_motifs.py`

---

## Usage

### Step 1: Generate Motifs

Run the following scripts to extract potential spanning motifs:

```bash
python3 1_aa_motifs_False.py
python3 2_nt_motifs_False.py
