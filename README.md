# dASPR

[![Build Status](https://github.com/Drew-Thomson/dASPR/actions/workflows/build.yml/badge.svg)](https://github.com/Drew-Thomson/dASPR/actions/workflows/build.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/language-C++-blue.svg)](https://isocpp.org/)
[![GitHub stars](https://img.shields.io/github/stars/Drew-Thomson/dASPR.svg?style=social&label=Star)](https://github.com/Drew-Thomson/dASPR/stargazers)

**dASPR** is a fast and accurate program for protein side-chain packing that supports both L- and D-amino acids. It is a modified version/fork of the original **FASPR** algorithm.

By mirroring the original Dunbrack rotamer library (copying each entry and inverting the sign of the dihedral angles) and incorporating D-amino acid logic, dASPR enables the modeling of heterochiral peptides and mixed-chirality protein structures.

---

### Credits & Acknowledgments
dASPR was built upon the open-source **FASPR** tool developed by Xiaoqiang Huang, Robin Pearce, and Yang Zhang. All credit for the core energy-based packing algorithm goes to the original authors. 

If you use dASPR, please cite this repository and the original FASPR papers:
1. Huang X, Pearce R, Zhang Y. *FASPR: an open-source tool for fast and accurate protein side-chain packing.* Bioinformatics (2020) 36: 3758-3765.
2. Huang X, Pearce R, Zhang Y. *Toward the Accuracy and Speed of Protein Side-Chain Packing: A Systematic Study on Rotamer Libraries.* J. Chem. Inf. Model. 2020; 60:410-420.

Copyright (c) Xiaoqiang Huang. FASPR is free to academic users. For suggestions regarding the original FASPR algorithm, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com. Feedback on the heterochiral additions in dASPR is welcomed at drew.thomson@glasgow.ac.uk.

---

## What's New in dASPR

* **Mixed Chirality Packing:** Support for packing mixed L- and D-amino acid residues. D-amino acids are appropriately formatted with an `ATOM` header in the PDB output.
* **Sequence Input:** Lowercase letters in a provided sequence string represent D-amino acid residues (Note: the original FASPR ability to fix side chain conformations via lowercase letters has been replaced by this feature). 
* **Intelligent Rotamer Library Resolution:** dASPR now searches for the `dun2010_mirror.bin` rotamer library using an intelligent fallback order, eliminating the need to explicitly pass `-r` if it's placed in a standard path.
* **Expanded Sequence Flags:** Support for `-f` and `--fasta` aliases when supplying sequence files.
* **Output Improvements:** Enhanced PDB output with correctly indexed atoms per residue.

## Usage

```bash
$path/dASPR -i input.pdb -o output.pdb [-s sequence.txt] [-r rotlib.bin]
```

- `-i`: The input backbone in PDB format for packing. There should be no missing main-chain atoms (N, CA, C and O). If side-chain atoms are included in input.pdb, they are ignored.
- `-o`: The repacked model output in PDB format. The residue positions are kept identical to the input PDB file.
- `-s` / `-f` / `--fasta` (optional): An input sequence to be repacked on the input protein backbone. The sequence should be a single line of one-letter amino acid codes. 
  - **Upper-case letters** represent standard L-amino acids.
  - **Lower-case letters** represent D-amino acids.
  - *Example:* `AcDefG` packs a sequence containing D-cysteine, D-glutamic acid, and D-phenylalanine. When the sequence is different from the structure, mutations are introduced.
- `-r` (optional): Path to the binary Dunbrack rotamer library file. **Important:** You must use `dun2010_mirror.bin`. Standard rotamer libraries will not work for D-amino acids.

### Rotamer Library Resolution Order
If `-r` is not provided, dASPR searches for `dun2010_mirror.bin` in the following order:
1. Environment variable (`DASPR_ROTLIB` or `FASPR_ROTLIB`)
2. The directory containing the `dASPR` executable
3. The current working directory (`./dun2010_mirror.bin`)
4. Standard system directories (`/usr/local/share/daspr/`, `/usr/local/bin/`, `/usr/share/daspr/`, `/opt/daspr/`)

## Installation

We recommend users download the source code and build the `dASPR` executable locally.

After downloading and unzipping the package, change into the dASPR directory and use the provided build script (on UNIX/Linux):

```bash
./build.sh
```

Alternatively, you can manually compile it:
```bash
g++ -O3 --fast-math -o dASPR src/*.cpp
```

*Note: For Mac users, use `-fast-math` or omit it entirely. Windows users will need to install a `g++` compiler first.*