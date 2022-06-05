# Secondary Structure Identification in RNA Sequences
## About The Project
Dr. Kirt Onthank, an Associate Professor of Biology at WWU, studies ocean acidification and its effects on the physiology of marine invertebrates, especially cephalopods. 
The goal of this project is to create a script using Python that, provided files containing potential edit sites and corresponding RNA sequences, will locate secondary structures at the site of each edit. This will facilitate Dr. Onthank’s ability to identify actual RNA edits, which will in turn aid his research on the effects of ocean acidification on RNA editing in octopuses.
## Getting Started
### Prerequisites
- [Python](https://docs.python.org/3/using/index.html)
- [Package Installer for Python (pip)](https://pip.pypa.io/en/stable/installation/)
### Clone Repository 
Utilizing either SSH or HTTPS, clone the [repository](https://github.com/wwu-cs/Octopus-RNA-Analysis).
#### SSH
```console
git clone git@github.com:wwu-cs/Octopus-RNA-Analysis.git
```
#### HTTPS
```console
git clone https://github.com/wwu-cs/Octopus-RNA-Analysis.git
```
### Install Dependencies
```console
pip install -r dependencies.txt
```
### Input
#### FASTA
Each item in the FASTA file contains a sequence id followed by the sequence. 
```
>lcl|TRINITY_DN119711_c0_g1_i1:220-300
ATGTCAAGCGACATCCCACGAGAATCTATGCCAGTGCTATATAGGATGGTGTCCATTCTGGTAGTGATTCATGGTGCTTAA
```
#### CSV
Each item in the CSV file must contain (at least) a sequence id (`orf`) and edit site (`pos`).
|orf                                 |pos |
|------------------------------------|----|
|lcl&#124;TRINITY_DN119711_c0_g1_i1:220-300   |32|
## Usage
### Parameters
Use the command below to view required arguments and additional flags. 
```console
python src/main/main.py -h 
```
```
usage: main.py [-h] -f  -c  -of  -st  [-ml] [-ll] [-nl] [-bl] [-nb]

Identify secondary structures in genetic sequence

required arguments:
  -f , --fasta          fasta file containing genetic sequences
  -c , --csv            csv file containing edit positions
  -of , --outputFile    name of output csv file
  -st , --structType    type of structure ('hairpin'|'int_loop'|'bulge')

optional arguments:
  -h, --help            show this help message and exit
  -ml , --minLength     minimum length of reverse complement (Default:5)
  -ll , --loopLength    max num of mismatches in loop (Default:1)
  -nl , --numLoops      max number of loops allowed in loop structure (Default:1)
  -bl , --bulgeLength   max num of mismatches in loop (Default:1)
  -nb , --numBulges     max number of loops allowed in loop structure (Default:1)
```
### Examples
A couple of examples utilizing the sample files in the data folder. 
```console
python src/main/main.py -f data/swissprotORF.fasta -c data/aes_profile.csv -of sample.csv -st 'hairpin' -ml 12
```
|id                                    |position|length|base_string   |base_string_loc|rev_comp      |rev_comp_loc|
|--------------------------------------|--------|------|--------------|---------------|--------------|------------|
|lcl&#124;TRINITY_DN149051_c0_g1_i1:7-1731  |718     |14    |GAATACACATTACT|[707, 720]     |AGTAATGTGTATTC|[356, 369]  |
|lcl&#124;TRINITY_DN139047_c0_g1_i1:115-2832|2457    |12    |CATCAGCAGCAG  |[2452, 2463]   |CTGCTGCTGATG  |[316, 327]  |
|lcl&#124;TRINITY_DN141652_c0_g1_i1:981-1274|230     |14    |AAAAAGAAAAAAAA|[225, 238]     |TTTTTTTTCTTTTT|[277, 290]  |
|lcl&#124;TRINITY_DN102602_c0_g1_i1:128-3073|2835    |12    |TCTCCAAGAGCG  |[2832, 2843]   |CGCTCTTGGAGA  |[1381, 1392]|
```console
python src/main/main.py -f data/swissprotORF.fasta -c data/aes_profile.csv -of sample.csv -st 'int_loop' -ml 8 -ll 2 -nl 1
```
|id                                   |position|length|base_string|base_string_loc|rev_comp   |rev_comp_loc|
|-------------------------------------|--------|------|-----------|---------------|-----------|------------|
|lcl&#124;TRINITY_DN115027_c0_g1_i1:269-580|119     |8     |AATT..TT   |[118, 125]     |AA..AATT   |[59, 66]    |
|lcl&#124;TRINITY_DN131245_c0_g1_i1:16-339 |224     |8     |TAAA.AGA   |[224, 231]     |TCT.TTTA   |[265, 272]  |
|lcl&#124;TRINITY_DN131271_c0_g1_i1:11-253 |195     |11    |TGCACA.CAGG|[185, 195]     |CCTG.TGTGCA|[142, 152]  |
|lcl&#124;TRINITY_DN141161_c0_g1_i1:70-897 |212     |8     |CACTTC.T   |[209, 216]     |A.GAAGTG   |[6, 13]     |
