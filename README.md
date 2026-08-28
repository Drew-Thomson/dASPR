# UPDATE

dASPR is a modified version of FASPR (original readme below) that allows packing of mixed l- and d-amino acid residues in protein structures. This was achieved by mirroring the original Dunbrack rotamer library, simply by copying each entry and inverting the sign of the dihedral angles.

The installation and usage of the algorithm are basically the same (but use the name dASPR instead of FASPR). Hence installation uses:

"<b>g++ -O3 --fast-math -o dASPR src/*.cpp</b>" 

It is also important that the dun2010_mirror.bin rotamer library be used. Other rotamer libraries won't work.

<b> $path/dASPR -i input.pdb -o output.pdb [-s sequence.txt] </b>

The other main difference is that lowercase letters in an input string represent d-amino acid residues. The ability to fix side chain conformations in the input model has been lost. D-amino acids are output with an ATOM header in the pdb file, not HETATM.

dASPR has only been tested on a linux system. The basic functionality for heterochiral peptides works for the cases examined so far, however the accuracy of mirroring the Dunbrack rotamer library remains to be evaluated.

Feedback on performance is welcomed at drew.thomson@glasgow.ac.uk

A communication is being prepared, but in the meantime if you use this please reference this repo, and of course cite the original FASPR paper!


# INTRODUCTION

FASPR is a fast and accurate program for protein side-chain packing, which is an important step in conventional, energy-based protein structure prediction and protein design.

# USAGE

<b> $path/FASPR -i input.pdb -o output.pdb [-s sequence.txt] [-r rotlib.bin] </b>

-i: the input backbone in pdb format for packing. There should be no missing main-chain atoms (N, CA, C and O). If side-chain atoms are included in input.pdb, they are ignored by FASPR.

-o: the repacked model output in pdb format by FASPR. The residue positions are kept identical to the input pdb file.

-s (optional): an input sequence to be repacked on the input protein backbone. The sequence should be written as one-single line of one-letter alphabet of amino acid types, e.g., ACDEFGHIKLMNPQRSTVWYYWVTSRQPNMLKIHGFEDCA. Only 20 canonical amino-acid types are allowed in the sequence. When the input sequence is the same as the one extracted from the pdb structure, the amino acid side-chain conformations are repacked only. When the input sequence is different, mutations will be introduced. Thus, <b>FASPR can be used to build mutant models efficiently</b>. If you want to fix the conformation of some residues during packing, you can specify them in lower-case letters while keeping other residues in upper-case, e.g., acdefghiklmnpqrstvwyYWVTSRQPNMLKIHGFEDCA.

-r (optional): path to the binary Dunbrack rotamer library file (`dun2010_mirror.bin`).

### Rotamer Library Resolution Order
FASPR searches for the rotamer library (`dun2010_mirror.bin`) in the following order of precedence:
1. Command-line argument (`-r` or `--rotlib`)
2. Environment variable (`DASPR_ROTLIB` or `FASPR_ROTLIB`)
3. Directory containing the executable binary
4. Current working directory (`./dun2010_mirror.bin`)
5. Standard system directories (`/usr/local/share/daspr/`, `/usr/local/bin/`, `/usr/share/daspr/`, `/opt/daspr/`)

# INSTALLATION
We recommend users to download the FASPR source-code package to your computer and build the FASPR executable on your own. After downloading and unzipping the package, change into the $path/FASPR/ directory and run "<b>g++ -O3 --fast-math -o FASPR src/*.cpp</b>" if you are working on UNIX or Linux. For Mac users, use "-fast-math" or ignore it. If you are working on the Windows system, you need to install the g++ compiler first.

# COPYRIGHT & CONTACT
Copyright (c) Xiaoqiang Huang. FASPR is free to academic users. For suggestions, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com.

# REFERENCES
1. Huang x, Pearce R, Zhang Y, FASPR: an open-source tool for fast and accurate protein side-chain packing. Bioinformatics (2020) 36: 3758-3765.
2. Huang X, Pearce R, Zhang Y. Toward the Accuracy and Speed of Protein Side-Chain Packing: A Systematic Study on Rotamer Libraries. J. Chem. Inf. Model. 2020; 60:410-420.
