# Secondary Structure Identification in RNA Sequences
## About The Project

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
>lcl|TRINITY_DN119711_c0_g1_i1:220-300 ORF1_TRINITY_DN119711_c0_g1_i1:219:299
ATGTCAAGCGACATCCCACGAGAATCTATGCCAGTGCTATATAGGATGGTGTCCATTCTGGTAGTGATTCATGGTGCTTAA
#### CSV
Each item in the CSV file must contain (at least) a sequence id (`orf`) and edit site (`pos`).
|orf                                 |pos |
|------------------------------------|----|
|lcl&#124;TRINITY_DN119711_c0_g1_i1:220-300   |32|
## Usage
Use the command below to view required arguments and additional flags. 
```console
python src/main/main.py -h 
```
```console
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

