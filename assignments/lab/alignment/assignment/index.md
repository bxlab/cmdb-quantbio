# Sequence Alignment

## Assignment Overview

The goal of today's lab is to implement the Needleman-Wunsch algorithm we discussed during class in a Python script. You will then use your implementation to align two DNA sequences, and then to align two protein sequences. For the assignment, you'll be aligning the CTCF gene between human and mouse genomes. Specifically, you will align both CTCF nucleotide and amino acid sequences from [GENCODE](https://www.gencodegenes.org/) (GENCODE version 38 for human and GENCODE version M27 for mouse).<br><br>

## Data

All the data you need for this assignment is in a zipped archive on github. You can download the data into your current directory with the following command:

```bash
wget https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/lab/alignment/extra_data/needleman-wunsch.tar.gz
```

After downloading the zipped folder, you'll need to extract it with `tar -zxvf <filename.tar.gz>`. You should get 5 files:
1. CTCF_38_M27_AA.faa
2. CTCF_38_M27_DNA.fna
3. BLOSUM62.txt
4. HOXD70.txt
5. assignment_framework.py
6. fasta.py

You'll be working with FASTA files in today's lab. This is the standard file format used to store nucleotide and amino acid sequences, and there's a good chance you'll encounter this format again in your Ph.D. (read more [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp)). Unfortunately, reading these files into Python in a meaningful way is not a trivial task, so we've put together some code to make this easier found in `fasta.py`. In your python script for this assignment, you can import the `readFASTA`  function using `from fasta import readFASTA`. If you're interested in how this function works, you're certainly encouraged to take a look at it in the `fasta.py` file. Note that for the import to work this script must be in the same directory as the script you are writing.

You also have a file `assignment_framework.py` which serves as a scaffold for this assignment. It includes some parts of the code already written as well as a roadmap for the steps to implementing the Needleman-Wunsch algorithm.
<br><br>

## Exercises

There are five parts to the exercise in this assignment:

1. Read in the parameters to your script
2. Initialize the scoring and traceback matrices
3. Fill in both matrices
4. Find the optimal sequence alignment
5. Write the resulting alignment to a file and report the requested statistics 
<br><br>

## Submission

Before you begin, create a `week10` folder in your `QBB2025-answers` repo. You will be expected to turn in 4 files for this assignment.

1. Python script that performs alignment based on user input
2. DNA alignment file
3. AA alignment file
4. `README.md` with values from step 1.5
<br><br>

### Exercise 1: Needleman-Wunsch Algorithm

Write a script to perform global alignment between two sequences using a given scoring matrix and gap penalty. Your script will take four inputs:
1. A FASTA-style file containing two sequences to align
2. A text file containing the scoring matrix you'd like to use for this alignment
3. The penalty for gaps in your alignment (so if users wanted to penalize gaps by subtracting 10 from the alignment score, they would input **-10**)
4. The filepath to write your alignment to

Additionally, your script should print out 5 things:

1. The number of gaps in the first sequence
2. The number of gaps in the second sequence
3. The percent sequence identity of the first sequence (the percentage of the bases that have an exact match)
4. The percent sequence identity of the second sequence
5. The score of the final alignment

You'll run your script twice:
1. Align the CTCF DNA transcript sequences from human and mouse using the [HOXD70](https://pubmed.ncbi.nlm.nih.gov/11928468/) scoring matrix and a gap penalty of **-300**.
2. Align the CTCF amino acid sequences from human and mouse using the [BLOSUM62](https://www.pnas.org/content/89/22/10915) scoring matrix and a gap penalty of **-10**.

**NOTE**: The DNA sequences are fairly long, and as such the DNA alignment may take a minute or more to run. We recommend testing your code with the protein alignment first (or even just a couple of small test sequences), and then running the DNA alignment when you're confident it's working.
<br><br>

#### **Step 1.1: Read in your parameters**

First `import sys` and then use `sys.argv` to read in the parameters you need to run your script from the command line. When reading in the fasta file, you can use `readFASTA` to store the sequence ids and sequences as follows:

```python
input_sequences = readFASTA(open(<fasta_file>))

seq1_id, sequence1 = input_sequences[0]
seq2_id, sequence2 = input_sequences[1]
```

For the scoring matrix, code has been provided to read these values into a dictionary keyed by a tuple of the pair of bases or amino acids being compared (e.g. `sigma[("A", "G")]`).
<br><br>

#### **Step 1.2: Initializing matrices**

You'll need two matrices to carry out the Needleman-Wunsch algorithm: an F-matrix that stores the score of each "optimal" sub-alignment (this is the one we created in class), as well as a traceback matrix that allows you to determine the optimal global alignment (as a path through this matrix). Initialize two empty matrices for these purposes.

*HINT*: With sequence 1 of length *m* and sequence 2 of length *n*, both matrices should be of size (*m+1*)×(*n+1*), to account for potential leading gaps in either sequence.
<br><br>

#### **Step 1.3: Populating the matrices**

Follow the steps of the Needleman-Wunsch algorithm discussed in class to populate the two matrices.

When generating the traceback matrix: if at any point there is a tie between aligning, a gap in sequence 1, or a gap in sequence 2, resolve the tie in this order: aligning -> gap in sequence 1 -> gap in sequence 2.
<br><br>

#### **Step 1.4: Find the optimal alignment**

Use the traceback matrix to find the optimal alignment between the two sequences. Start at the bottom right corner and follow a path backwards through the traceback matrix until you reach the top left of the matrix, building the alignment as you go. You should end up with two strings of the same length, one for each sequence. Gaps in the sequences should be denoted with a hyphen (`-`). For example, if your input sequences were `TACGATTA` and `ATTAACTTA` your final alignment might look something like:

```
Sequence 1 alignment: '--TACGA-TTA'
Sequence 2 alignment: 'ATTA--ACTTA'
```

*HINT*: A `while` loop will probably be helpful for this part.
<br><br>

#### **Step 1.5: Write the alignment to the output**

Write the alignment to the output file specified in the command line.

**ALSO**, make sure your script prints out the additional requested information:
1. The number of gaps in each sequence
2. The percent sequence identity for each sequence
3. The score of the  alignment 

For *both* alignments (DNA and AA), record these values in your `README.md`.
<br><br>

### Exercise 2: Smith-Waterman (OPTIONAL)

If you finish early and you want to try something a little bit harder, try to implement the Smith-Waterman algorithm for local alignment in a new script (or as a separate function in the same script).
<br><br>

## Submission

1. Python script that performs alignment based on user input (**8 points total**)
	* Step 1.1: Read in inputs (**1.5 points**)
	* Step 1.2: Initialize matrices (**1 point**)
	* Step 1.3: Populate matrices (**2 points**)
	* Step 1.4: Perform traceback (**2 points**) 
	* Step 1.5: Write alignment to file
	* Step 1.5: Print out number of gaps and percent sequence identity in each sequence and alignment score (**1.5 points**)
2. DNA alignment file (**0.5 points**)
3. AA alignment file (**0.5 points**)
4. `README.md` with alignment score, sequence identities, and number of gaps in each sequence *FOR EACH ALIGNMENT* (**1 point**)

**Total Points: 10**

<br><br>
