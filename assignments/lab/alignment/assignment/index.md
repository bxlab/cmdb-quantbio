# Assignment 2: Sequence Alignment
Assignment Date: Friday, Sept. 16, 2022 <br>
Due Date: Friday, Sept. 23, 2022 @ 1:00pm ET <br>

## Lecture and Live-Coding

**Slides** are available here: [Lecture slides](https://docs.google.com/presentation/d/1IRm-2vsaJLWN2YV0us_UHHwVVDEfrvXu8zW-9zc0Jec/edit?usp=sharing)

**Blank Live-coding script** is available here: `~/cmdb-quantbio/assignments/lab/alignment/slides_asynchronous_or_livecoding_resources/needleman-wunsch_livecoding_empty.py`

**Master Live-coding script** is available here: `~/cmdb-quantbio/assignments/lab/alignment/slides_asynchronous_or_livecoding_resources/needleman-wunsch_livecoding_master.py`<br><br>


## Assignment Overview

The goal of today's lab is to implement the Needleman-Wunsch algorithm we discussed during class in a Python script. You will then use your implementation to align two DNA sequences, and then to align two protein sequences. For the assignment, you'll be aligning the CTCF gene between human and mouse genomes. Specifically, you will align both CTCF nucleotide and amino acid sequences from [GENCODE](https://www.gencodegenes.org/) (GENCODE version 38 for human and GENCODE version M27 for mouse).<br><br>

## Data

All the data you need for this assignment is in a zipped folder here: `~/cmdb-quantbio/assignments/lab/alignment/extra_data/needleman-wunsch.tar.gz`. Copy this file to the `answers` directory you made for this assignment.

After copying the zipped folder, you'll need to extract it with `tar -zxvf <filename.tar.gz>`. You should get 4 files:
1. CTCF_38_M27_AA.faa
2. CTCF_38_M27_DNA.fna
3. BLOSUM62.txt
4. HOXD70.txt

You'll be working FASTA files in today's lab. This is the standard file format used to store nucleotide and amino acid sequences, and there's a good chance you'll encounter this format again in your Ph.D. (read more [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp)). Unfortunately, reading these files into Python in a meaningful way is not a trivial task, so we've put together some code to make this easier. Copy the `~/cmdb-quantbio/resources/code/fasta.py` file into your current assignment directory. In your python script for this assignment, you can import the `readFASTA`  function using `from fasta import readFASTA`. If you're interested in how this function works, you're certainly encouraged to take a look at it in the `fasta.py` file.<br><br>

## Assignment

Write a script to perform global alignment between two sequences using a given scoring matrix and gap penalty. Your script will take four inputs:
1. A FASTA-style file containing two sequences to align
2. A text file containing the scoring matrix you'd like to use for this alignment
3. The penalty for gaps in your alignment
4. The filepath to write your alignment to

Additionally, your script should print out the number of gaps in the first sequence, the number of gaps in the second sequence, and the score of the final alignment.

You'll run your script twice:
1. Align the CTCF DNA transcript sequences from human and mouse using the [HOXD70](https://pubmed.ncbi.nlm.nih.gov/11928468/) scoring matrix and a gap penalty of **300**.
2. Align the CTCF amino acid sequences from human and mouse using the [BLOSUM62](https://www.pnas.org/content/89/22/10915) scoring matrix and a gap penalty of **10**.

**NOTE**: The DNA sequences are fairly long, and as such the DNA alignment may take a few minutes to run. We recommend testing your code with the protein alignment first (or even just a couple of small test sequences), and then running the DNA alignment when you're confident it's working.<br><br>

#### Step 1: Read in your parameters

Use `sys.argv` to read in the parameters you need to run your script from the command line. When reading in the fasta file, you can use `readFASTA` to store the sequence ids and sequences as follows:

```
input_sequences = readFASTA(open(<fasta_file>))

seq1_id, sequence1 = input_sequences[0]
seq2_id, sequence2 = input_sequences[1]
```

For the scoring matrix, it probably makes sense to store it in a numpy array.<br><br>

#### Step 2: Initializing matrices

You'll need two matrices to carry out the Needleman-Wunsch algorithm: an F-matrix that stores the score of each "optimal" sub-alignment (this is the one we created in class), as well as a traceback matrix that allows you to determine the optimal global alignment (as a path through this matrix). Initialize two empty matrices for these purposes.

*HINT*: With sequence 1 of length *m* and sequence 2 of length *n*, both matrices should be of size (*m+1*)×(*n+1*), to account for potential leading gaps in either sequence.<br><br>

#### Step 3: Populating the matrices

Follow the steps of the needleman-wunsch algorithm discussed in class to populate the two matrices.

When generating the traceback matrix: if at any point there is a tie between aligning, a gap in sequence 1, or a gap in sequence 2, resolve the tie in the order (aligning -> gap in sequence 1 -> gap in sequence 2).<br><br>

#### Step 4: Find the optimal alignment

Use the traceback matrix to find the optimal alignment between the two sequences. Start at the bottom right corner and follow a path backwards through the traceback matrix until you reach the top left of the matrix, building the alignment as you go. You should end up with two strings of the same length, one for each sequence. Gaps in the sequences should be denoted with a hyphen (`-`). For example, if your input sequences were `TACGATTA` and `ATTAACTTA` your final alignment might look something like:

```
Sequence 1 alignment: '--TACGA-TTA'
Sequence 2 alignment: 'ATTA--ACTTA'
```

*HINT*: A `while` loop will probably be helpful for this part.<br><br>

#### Step 5: Write the alignment to the output

Write the alignment to the output file specified in the command line, and print out the additional information requested (the number of gaps in each sequence and the score of the  alignment).<br><br>


## Submission

For this assignment, you should submit four things:
1. Your Needleman-Wunsch Python Script
2. A text file containing your DNA sequence alignment
3. A text file containing your amino acid sequence alignment
4. A markdown or text file detailing the number of gaps in each sequence (including leading and trailing gaps) and alignment score for both of the above alignments<br><br>

## Just for fun

If you finish early and you want to try something a little bit harder, try to implement the Smith-Waterman algorithm for local alignment in a new script.<br><br>