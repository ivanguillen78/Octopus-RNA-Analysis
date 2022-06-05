# Secondary Structure Identification in RNA
## Background
### Introduction
At Walla Walla University, and many other universities, senior computer science majors are required to work on a capstone project to show what they're learned over the course of their college experience. When it came time to select my project, I was unsure about what I wanted to work on. Professor James Foster, one of my advisors for this course, introduced me to Dr. Kirt Onthank, an associate professor of Biology at WWU. Dr. Onthank's research includes the study of ocean acidification and its impact on marine life, especially cephalopods. For this project, Dr. Onthank was looking for someone to implement a program that could identify secondary structures in octopus RNA sequences around a specific edit site. This would allow him to identify real edits as opposed to false positives, which, in turn, would allow him to observe the effects of ocean acidification on RNA editing in octopuses. 
### What are secondary structures? 
Nucleic acid secondary structures are base-pairing interactions that occur within the same sequence. A base pair is comprised of either A <-> T or G <-> C. For this project, my focus was on three types of structures:  

**Hairpin loops** are a common type of secondary structure that are created when a sequence of RNA folds upon itself and forms base pairs with another section of the same sequence.  

![hairpin loop secondary structure](assets/hairpin.png "Hairpin")

**Internal Loops** are similar, but feature a short sequence of unpaired bases within a larger sequence of paired bases.

![internal loop secondary structure](assets/int_loop.png "Internal Loop")

**Bulges** are also similar, but feature regions on one side of the folded structure that have extra bases with no corresponding bases on the opposite side.

![bulge secondary structure](assets/bulge.png "Bulge")

As seen in the images above, for every structure, there is a "base" string of bases and a "reverse complement" string of bases. The base string is continually expanded around the provided edit site, provided that the reverse complement of that base string exists somewhere else in the sequence. 
### Reverse Complements
To find the reverse complement of a string of bases, we simply:
#### Reverse the String
The string "AGATC" would become "CTAGA". 
#### Switch each Base to its Complement
As stated above, A's become T's and G's become C's (and vice versa). The string "CTAGA" would become "GATCT" Thus, the reverse complement of "AGATC" is "GATCT". 
## Input
There are two files that we are working with. A FASTA file containing sequences and a CSV file containing edit sites corresponding to the sequences in the FASTA file. 
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
## Tools
### Python vs R
Prior to beginning this project, I was unsure about what language to use. I had some experience with Python and R from previous courses and I was leaning towards one of the two. Both languages certainly could have been used for this program. However, I found it difficult to justify R’s lack of speed in comparison to Python, especially when dealing with large datasets. Additionally, while both languages featured imports that would assist in the locating of secondary structures, I found that Python’s libraries were much more intuitive, complete, and accessible. Thus, I decided to implement this project in Python, in addition to using Git as source control and Visual Studio Code as my IDE. 
### Additional Libraries
- [biopython](https://biopython.org)
- [fastaparser](https://pypi.org/project/fastaparser/)
- [alive_progress](https://pypi.org/project/alive-progress/)
## Implementation
Because our goal is to expand around the edit site, my initial thought was to symmetrically expand outwards from that site. The implementation looked something similar to:
```
# sequence (seq) and edit site (E)

Set L to E - 1 and R to E + 1
If the reverse complement of seq[L:E] exists in the sequence:
    Decrement L
If the reverse complement of seq[E:R] exists in the sequence:
    Increment R
Set longest to seq[L:R]
```
Essentially, we would set some pointers that would iterate outwards until they could no longer do so. Although this does work, the implementation is obviously flawed. Because we are search to the left **and** to the right each iteration, there is a possibility that decrementing L or incrementing R on one iteration will hinder future iterations. Thus, we must ensure that we account for all possibilities. 

With some guidance from Dr. Preston Carman, I implemented the following algorithm:
```
# Note: Some steps, such as bounds checking, have been omitted for simplicity

# Step 1: Find Longest Reverse Complement to the Left
Set L = E
While the reverse complement of seq[L:E] exists in either seq[:L] or seq[E:],
    Decrement L by 1
Set longest to seq[L:E]

# Step 2: Find Longest Reverse Complement
Set R = E + 1
While L ≤ E,
    While the reverse complement of seq[L:R] exists in either seq[:L] or seq[R:],
        Set longest to seq[L:R]
        Increment R by 1
    Increment L and R by 1
Return longest
``` 
### Step 1
As seen above, we cannot search to the left and right at the same time. Therefore, we must choose to go in one direction first. In this case, the first step is to find the longest reverse complement to the left. We do this by initializing our left index `L` to `E`. We continually decrement `L` while the reverse complement of `seq[L:E]` exists in the sequence, making sure to avoid overlap by searching only to the left of `L` and to the right of `E`. Once L cannot go any farther to the left, we set `longest` to `seq[L:E]`.
### Step 2
The second step is similar to the first, but now our goal is to find the longest overall reverse complement, which means we must search to the right. While we could simply add on to `longest` from Step 1, we must remember that the addition of a base on either end (the left end in this case) could hinder our search to the right. Therefore, to account for all possibilities, we iterate through every index from `L` to `E`, ensuring that all base strings that include `E` have been checked. Additionally, because our goal is to produce the longest structure, we increment `R` alongside `L` at the end of each iteration to ensure that at each iteration, we are searching for a reverse complement of **at least** the same length as the current longest.
### Visualization
![algorithm visualization](assets/algo_visual.png "Algorithm Visualization")

### Additional Modifications for Internal Loops and Bulges
While the base algorithm is similar for all three types, additional modifications were needed to account for unpaired/extra bases in internal loops and bulges. Regular expression searches were used instead of Python’s built-in string.find() and periods (“.”) were used to simulate gaps between paired bases. A boolean flag was set if a jump over the unpaired/extra bases could be made. The jump would only be made if necessary.

## Usage
### Parameters
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
## Conclusion
By iteratively improving upon the algorithm, this project successfully allows users to efficiently identify three common types of secondary structures. This secondary structure identification will play a key role in allowing Dr. Onthank, and other researchers, to sort out actual edits from false positives. 

While the program itself might not be inherently complex and utilizes a few of Python's built-in libraries to accomplish its goal, I am proud of the work I put into this project. Before starting college, I had zero coding experience. My first experience was on the first day of my first year at WWU. I am thankful for all that I have learned at WWU and I hope to continue learning, improving, and expanding my skill set. For anyone that might be nervous to begin coding, I'm certain that if I can do it, you can as well! 