# Assignment 7: Global and Local Alignment
Assignment Date: Friday, Oct. 22, 2021 <br>
Due Date: Friday, Oct. 29, 2021 @ 1pm ET <br>

## Lecture

[Lecture slides](https://docs.google.com/presentation/d/1IRm-2vsaJLWN2YV0us_UHHwVVDEfrvXu8zW-9zc0Jec/edit?usp=sharing)

Live-coding blank notebook: `wget https://github.com/bxlab/qbb2021/raw/main/week7/needleman-wunsch_livecoding_empty.py`

Live-coding master notebook: `wget https://github.com/bxlab/qbb2021/raw/main/week7/needleman-wunsch_livecoding_master.py`


## Basic Exercises: The Needleman-Wunsch Algorithm

The goal of today's lab is to implement the Needleman-Wunsch algorithm we discussed during class in a Python script. You will then use your implementation to align two DNA sequences, and then to align two protein sequences. Given our discussion regarding chromosome capture in previous lectures, you'll be aligning the CTCF gene between human and mouse genomes. Specifically, you will align both CTCF nucleotide and amino acid sequences from [GENCODE](https://www.gencodegenes.org/) (GENCODE version 38 for human and GENCODE version M27 for mouse).

Write a script to perform global alignment between two sequences using a given scoring matrix and gap penalty. Your script will take four inputs:
1. A FASTA-style file containing two sequences to align
2. A text file containing the scoring matrix you'd like to use for this alignment
3. The penalty for gaps in your alignment
4. The filepath to write your alignment to

Additionally, your script should print out the number of gaps in the first sequence, the number of gaps in the second sequence, and the score of the final alignment.

You'll run your script twice:
1. Align the CTCF DNA transcript sequences from human and mouse using the [HOXD70](https://pubmed.ncbi.nlm.nih.gov/11928468/) scoring matrix and a gap penalty of 300.
2. Align the CTCF amino acid sequences from human and mouse using the [BLOSUM62](https://www.pnas.org/content/89/22/10915) scoring matrix and a gap penalty of 10.

**NOTE**: The DNA sequences are fairly long, and as such the DNA alignment may take a few minutes to run. We recommend testing your code with the protein alignment first, and then running the DNA alignment when you're confident it's working.

### Getting your data

Everything you need for this assignment is located [here](https://github.com/bxlab/qbb2021/raw/main/week7/needleman-wunsch.tar.gz) in a tarball. Download it. To extract it, run the following command: `tar -xzf needleman-wunsch.tar.gz`. You should get 4 files:
1. CTCF_38_M27_AA.faa
2. CTCF_38_M27_DNA.fna
3. BLOSUM62.txt
4. HOXD70.txt


### Submitting your assignment

For this assignment, you should submit four things:
1. Your Needleman-Wunsch Python Script
2. A text file containing your DNA sequence alignment
3. A text file containing your amino acid sequence alignment
4. A markdown or text file detailing the number of gaps and alignment score for both of the above alignments

### Building your script

#### Step 1: Read in your parameters

Use `sys.argv` to read in the parameters you need to run your script from the command line. When reading in the fasta file, it will probably be useful to use the `FASTAReader` function you made during Bootcamp. Either `import` that function, or just copy it into your script.

For the scoring matrix, it makes sense to store it in either a pandas dataframe or a numpy array

#### Step 2: Initializing matrices

You'll need two matrices to carry out the Needleman-Wunsch algorithm: an F-matrix that stores the score of each "optimal" sub-alignment, and a traceback matrix that allows you to determine the optimal global alignment (as a path through this matrix). Initialize two empty matrices for these purposes.

*Hint*: With sequence 1 of length *m* and sequence 2 of length *n*, both matrices should be of size (*m+1*)×(*n+1*), to account for potential leading gaps in either sequence.

#### Step 3: Populating the matrices

Follow the steps of the needleman-wunsch algorithm discussed in class to populate the two matrices.

When generating the traceback matrix: if at any point there is a tie between aligning, a gap in sequence 1, or a gap in sequence 2, resolve the tie in the order (aligning -> gap in sequence 1 -> gap in sequence 2).

#### Step 4: Find the optimal alignment

Use the traceback matrix to find the optimal alignment between the two sequences. Start at the bottom right corner and follow a path backwards through the traceback matrix until you reach the top left of the matrix. You should end up with two strings of the same length, one for each sequence. Gaps in the sequences should be denoted with a hyphen (`-`).

#### Step 5: Write the alignment to the output

Write the alignment to the output file specified in the command line, and print out the additional information requested.
