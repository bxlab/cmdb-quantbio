# Sequence Alignment

## Assignment Overview

The goal of today's lab is to implement the Needleman-Wunsch algorithm we discussed during class in a Python script. You will then use your implementation to align two DNA sequences, and then to align two protein sequences. For the assignment, you'll be aligning the CTCF gene between human and mouse genomes. Specifically, you will align both CTCF nucleotide and amino acid sequences from [GENCODE](https://www.gencodegenes.org/) (GENCODE version 38 for human and GENCODE version M27 for mouse).<br><br>

## Data

All the data you need for this assignment is in a zipped folder on Dropbox. You can download the data into your current directory with the following command:

```
wget --content-disposition "https://www.dropbox.com/scl/fi/cq6mj7h34pnzkkabgz4ks/needleman-wunsch.tar.gz?rlkey=00wi3ypsfcfnml1psfyl7p9vu&dl=0"
```

After downloading the zipped folder, you'll need to extract it with `tar -zxvf <filename.tar.gz>`. You should get 4 files:
1. CTCF_38_M27_AA.faa
2. CTCF_38_M27_DNA.fna
3. BLOSUM62.txt
4. HOXD70.txt

You'll be working FASTA files in today's lab. This is the standard file format used to store nucleotide and amino acid sequences, and there's a good chance you'll encounter this format again in your Ph.D. (read more [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp)). Unfortunately, reading these files into Python in a meaningful way is not a trivial task, so we've put together some code to make this easier. Copy the `~/cmdb-quantbio/resources/code/fasta.py` file into your current assignment directory. In your python script for this assignment, you can import the `readFASTA`  function using `from fasta import readFASTA`. If you're interested in how this function works, you're certainly encouraged to take a look at it in the `fasta.py` file.<br><br>

## Exercises

### Exercise 1: Needleman-Wunsch Algorithm

Write a script to perform global alignment between two sequences using a given scoring matrix and gap penalty. Your script will take four inputs:
1. A FASTA-style file containing two sequences to align
2. A text file containing the scoring matrix you'd like to use for this alignment
3. The penalty for gaps in your alignment (so if users wanted to penalize gaps by subtracting 10 from the alignment score, they would input **-10**)
4. The filepath to write your alignment to

Additionally, your script should print out the number of gaps in the first sequence, the number of gaps in the second sequence, and the score of the final alignment.

You'll run your script twice:
1. Align the CTCF DNA transcript sequences from human and mouse using the [HOXD70](https://pubmed.ncbi.nlm.nih.gov/11928468/) scoring matrix and a gap penalty of **-300**.
2. Align the CTCF amino acid sequences from human and mouse using the [BLOSUM62](https://www.pnas.org/content/89/22/10915) scoring matrix and a gap penalty of **-10**.

**NOTE**: The DNA sequences are fairly long, and as such the DNA alignment may take a few minutes to run. We recommend testing your code with the protein alignment first (or even just a couple of small test sequences), and then running the DNA alignment when you're confident it's working.<br><br>

#### **Step 1.1: Read in your parameters**

Use `sys.argv` to read in the parameters you need to run your script from the command line. When reading in the fasta file, you can use `readFASTA` to store the sequence ids and sequences as follows:

```
input_sequences = readFASTA(open(<fasta_file>))

seq1_id, sequence1 = input_sequences[0]
seq2_id, sequence2 = input_sequences[1]
```

For the scoring matrix, it probably makes sense to store it in a pandas dataframe.

*HINT*: Uh oh! Looks like the scoring matrices don't use a consistent field separator... Maybe `pd.read_csv()` has an argument that lets you separate columns based on an arbitrary amount of whitespace...<br><br>

#### **Step 1.2: Initializing matrices**

You'll need two matrices to carry out the Needleman-Wunsch algorithm: an F-matrix that stores the score of each "optimal" sub-alignment (this is the one we created in class), as well as a traceback matrix that allows you to determine the optimal global alignment (as a path through this matrix). Initialize two empty matrices for these purposes.

*HINT*: With sequence 1 of length *m* and sequence 2 of length *n*, both matrices should be of size (*m+1*)×(*n+1*), to account for potential leading gaps in either sequence.<br><br>

#### **Step 1.3: Populating the matrices**

Follow the steps of the needleman-wunsch algorithm discussed in class to populate the two matrices.

When generating the traceback matrix: if at any point there is a tie between aligning, a gap in sequence 1, or a gap in sequence 2, resolve the tie in the order (aligning -> gap in sequence 1 -> gap in sequence 2).<br><br>

#### **Step 1.4: Find the optimal alignment**

Use the traceback matrix to find the optimal alignment between the two sequences. Start at the bottom right corner and follow a path backwards through the traceback matrix until you reach the top left of the matrix, building the alignment as you go. You should end up with two strings of the same length, one for each sequence. Gaps in the sequences should be denoted with a hyphen (`-`). For example, if your input sequences were `TACGATTA` and `ATTAACTTA` your final alignment might look something like:

```
Sequence 1 alignment: '--TACGA-TTA'
Sequence 2 alignment: 'ATTA--ACTTA'
```

*HINT*: A `while` loop will probably be helpful for this part.<br><br>

#### **Step 1.5: Write the alignment to the output**

Write the alignment to the output file specified in the command line.

**ALSO**, make sure your script prints out the additional requested information:
1. the number of gaps in each sequence and
2. the score of the  alignment 

For *both* alignments (DNA and AA), record these values in your `README.md`.<br><br>

### Exercise 2: Smith-Waterman (OPTIONAL)

If you finish early and you want to try something a little bit harder, try to implement the Smith-Waterman algorithm for local alignment in a new script (or as a separate function in the same script).

## Submission

1. Python script that performs alignment based on user input (**8 points total**)
	* Step 1.1: Read in inputs (**2 points**)
	* Step 1.2: Initialize matrices (**1 point**)
	* Step 1.3: Populate matrices (**2 points**)
	* Step 1.4: Perform traceback (**2 points**) 
	* Step 1.5: Write alignment to file (**0.5 point**)
	* Step 1.5: Print out alignment score and number of gaps in each sequence (**0.5 points**)
2. DNA alignment file (**0.5 point**)
3. AA alignment file (**0.5 point**)
4. `README.md` with alignment score and number of gaps in each sequence *FOR EACH ALIGNMENT* (**1 point**)

**Total Points: 10**

<br><br>
