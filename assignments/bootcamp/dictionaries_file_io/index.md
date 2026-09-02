# Assigment 2 - Python dictionaries and file I/O

### Overview
------------

Most FASTA files you will encounter will contain multiple sequences. Further, for genes you will often get 3' and 5' untranslated sequences (UTRs) flanking the coding sequence. In this exercise you will be writing a function for parsing a multiple-sequence FASTA file and translating each sequence, locating the coding sequence in between the UTRs. Finally, you will look at the prevalence of amino acids across the different sequences.

Biological Learning Objectives

- Use a codon table and DNA sequences to identify and translate coding sequences
- Count and compare amino acid usage for the population of sequences 

Computational Learning Objectives

- Create and modify a dicitionary
- Use keys to retrieve values stored in a dictionary
- Step through each key and value in a dictionary inside a `for` loop
- Use conditionals to alter the behaviour or a `for` or `while` loop, terminating early if necessary
- Open and read a text file
- Write to a text file

### Instructions
----------------

- Save a copy of this notebook as `~/qbXX-answers/day2-afternoon/python_dictionaries.ipynb`.
- Copy the files `sequences.fa` and `codons.tsv` into `~/qbXX-answers/day2-afternoon/` (but do not submit them with your answer)
- Fill in answers in the available code/markdown cells below.
- Remember to comment your code to help yourself and us know what each part is intended to do.

### What to turn in
-------------------

- This filled-in notebook
- Your codon usage file

### Scoring
-----------

- 2.0 pts - Multiple-sequence fasta load function
- 0.5 pts - Codon table load function
- 4.0 pts - Translation function
    - 1.5 pts - Loop stepping through sequence with correct reading frame
    - 1.0 pts - Skipping 5' UTR
    - 0.5 pts - Converting coding sequence to amino acids
    - 1.0 pts - Skipping 3' UTR
- 2.0 pts - Counting amino acid usage
- 0.5 pts - Converting counts to percentages
- 1.0 pts - Writing usage results file

10 pts total

-------------------

1. Start by setting the correct working directory (this is important if you don't want to use full paths for your file names)

```python
%cd ~/qxx-answers/day2-afternoon/
```

2. Building on the code for reading in a single FASTA sequence, adapt it to read in multiple sequences from a single FASTA file, storing each sequence in a dictionary using the sequence name as the key. Wrap this code in a function such that it takes a file name as the only function argument and returns the dictionary of sequences.

- To check if a line represents the start of a new sequence, consider using the string method `.startswith()`
- Don't forget to close your filestream

3. Wrap the code for reading in the codon table into a function, taking a file name in as the argument and returning the dictionary of codon/amin acid pairs.

4. Write a function for translating the CDS sequences into amino acid sequences. This function will need to take in two arguments, the codon table and the DNA sequence. The DNA sequences contain untranslated sequences (UTRs) at the start and end so you will need to step through the sequences to find the first methionine (M). Likewise, you will need to stop when you encounter the first step codon (*) rather than translating through the end of the sequence.

- Using a `while` loop may be useful for this task, but it is not required as `for` loops can also work
- You will need to consider three different parts reading the sequence:
    1. Have you reached the start of the coding sequence
    2. Do you need to record the current codon's amino acid
    3. Have you reached the end of the coding sequence

5. Finally, put it all together, loading in the FASTA sequences and codon table, and translating the into amino acids. Once you have the amino acid sequences, count the number of times each amino acid is used. Finally, convert these counts into percentages and write them to a tab-separated file with the first column being the amino acid letter and the second column being the percent usage.

- To get the percentages, it will helpful to keep a running total of the number of amino acids as you find the counts
- The dictionary method `.setdefault` may be useful for intializing you count dictionary for each new amino acid
- Don't forget to close your filestream
