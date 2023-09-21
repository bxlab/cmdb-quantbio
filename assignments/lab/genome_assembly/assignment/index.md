# Genome Assembly

## Assignment Overview

***FILL IN***

## Exercises 

### Exercise 1: Coverage Simulator 

#### **Step 1.1**
How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?

#### **Step 1.2**
Using Python, simulate sequencing 3x coverage of a 1Mbp genome with 100bp reads and plot the histogram of coverage. Note you do not need to actually output the sequences of the reads, you can just uniform randomly sample positions in the genome and record the coverage. You do not need to consider the strand of each read. The start position of each read should have a uniform random probabilty at each possible starting position (1 through 999,901). You can record the coverage in an array of 1M positions. Overlay the histogram with a Poisson distribution with lambda=3. Also overlay the distribution with a Normal distribution with a mean of 3 and a standard deviation of 1.73 (which is the square root of 3). Here is the pseudocode for the simulator:

```
num_reads = calculate_number_of_reads(genomesize, readlength, coverage)
​
## use an array to keep track of the coverage at each position in the genome
genome_coverage = initialize_array_with_zero(genomesize)
​
for (i = 0; i < num_reads; i++)
{
  startpos = uniform_random(1,genomelength-readlength)
  endpos = startpos + readlength - 1
  for (x = startpos; x <= endpos; x++)
  {
    genomecoverage[x] = genomecoverage[x] + 1
  }
}
​
maxcoverage = max(genomecoverage)
​
## use an array count how many positions have 0x coverage, have 1x coverage, have 2x coverage, ...
histogram = initialize_array_with_zero(maxcoverage)
​
for (x = 0; x < genomelength; x++)
{
  cov = genomecoverage[x]
  histogram[cov] = histogram[cov] + 1
}
​
## now plot the histogram
...
```

#### **Step 1.3**
Using the histogram from Q1.2, how much of the genome has not been sequenced (has 0x coverage)? How well does this match Poisson expectations? How well does the normal distribution fit the data?

#### **Step 1.4**
Now repeat the analysis with 10x coverage: 1. simulate the appropriate number of reads, 2. make a histogram, 3. overlay a Poisson distribution with lambda=10, 4. overlay with a Normal distribution with mean=10, standard deviation=3.16. 5. compute the number of bases with 0x coverage, and 6. evaluate how well it matches the Poisson expectation and Normal expectations.

#### **Step 1.5**
Now repeat the analysis with 30x coverage: 1. simulate the appropriate number of reads, 2. make a histogram, 3. overlay a Poisson distribution with lambda=30, 4. overlay with a Normal distribution with mean=30, standard deviation=5.47 5. compute the number of bases with 0x coverage, and 6. evaluate how well it matches the Poisson expectation and Normal expectations.


### Exercise 2. de Bruijn Graph construction

#### **Step 2.1**

Write a script to draw the de Bruijn graph for the following reads using k=3 (assume all reads are from the forward strand, no sequencing errors, complete coverage of the genome). You may find [graphviz](https://graphviz.org/) to be helpful (see below).

```
ATTCA
ATTGA
CATTG
CTTAT
GATTG
TATTT
TCATT
TCTTA
TGATT
TTATT
TTCAT
TTCTT
TTGAT
```
Here is some pseudocode for de Bruijn graph construction: 

```
graph = {}

for each read
  for (i = 0; i < read.len() - k - 1; i++)
     kmer1 = read.substr(i, k)
     kmer2 = read.substr(i+1, k)
     graph[kmer1][kmer2] = 1

for each kmer1 in graph
   for each kmer2 in graph[kmer1]
      print "kmer1 -> kmer2"
```

#### **Step 2.2**

Assume that the maximum number of occurrences of any 3-mer in the actual genome is 4 using the k-mers from Q2.1. Write one possible genome sequence

#### **Step 2.3**

In a few sentences, what would it take to accurately reconstruct the sequence of the genome?

### Exercise 3: Why Genomics?

#### **Step 3.1**

Use ChatGPT (or Bard or your favorite LLM) to write an essay on why you are interested in genomics. Make sure to ask for references. Make sure to include both your prompt(s) and the output from the LLM
​
#### **Step 3.2**

Comment on the output from the LLM - does it make logical sense, does it include any phrases you would not have written, do the citated papers exist and support the claims from the LLM?

## Bonus Exercises

### Exercise 4: K-mer Uniqueness 

Download the human chomosome 22 from here: [https://schatz-lab.org/appliedgenomics2023/assignments/assignment1/chr22.fa.gz](https://schatz-lab.org/appliedgenomics2023/assignments/assignment1/chr22.fa.gz)

A kmer is a substring of length k. For example, the string GATTACA, has these 3-mers: GAT, ATT, TTA, TAC, ACA
​
A string of length G has G - k + 1 kmers. For long strings, G - k + 1 is nearly the same as G e.g. for human using 19mers, 3,000,000,000 vs 2,999,999,986
​
While a string of length G has G-k+1 kmers, there may be many fewer *distinct* kmers. For example, in the string "GCATCATCATCATCATCATCAT..." the kmers are: GCA, CAT, ATC, TCA, CAT, ATC, TCA, CAT, ATC, TCA, CAT, ... As such there are only 4 disinct kmers (GCA, CAT, ATC, TCA). Of these GCA occurs once and the others occur many times.

##### **Step 4.1** 

How many As, Cs, Gs, Ts and Ns are found in the entire chromosome? If needed convert lowercase letters to uppercase, and any other character can be converted to N.

#### **Step 4.2**

In the language of your choice, tally the frequency of 19-mers in the chromosome, and output the kmer frequency spectrum upto 1000 e.g. how many kmers occur 1 time, how many occur 2 times, how many occur 3 times, etc. For this, convert lowercase letters to uppercase, and any character that is not ACG or T can be converted to A (especially N characters). We recommend you use a dictionary (or hash table) to tally the frequencies using this pseudocode. In your writeup, show the kmer frequency spectrum for 1 to 20, e.g. how many kmers occur 1 time, how many occur 2, ..., how many occur 20 times. This can be done with the unix command `head`:

```
## initialize kmer length
k=19
​
## read genome, convert to upper canse and convert non-DNA to 'A'
genome_string = read_from_file("chr22.fa")
​
## dictionary that maps a kmer (like GAT) to a frequency (like 3)
kmer_frequency = initialize_dictionary()
​
## now scan the genome, extract kmers, and tally up their frequencies
for(i = 0; i < length(genome_string) - k + 1; i++)
{
  kmer = substring(genome_string, i, k)  
  kmer_frequency[kmer] = kmer_frequency[kmer] + 1
}
​
## now tally the frequencies in a dictionary that maps kmer frequency to count
## also determine the maximum kmer frequency
​
tally = initialize_dictionary()
all_kmers = kmer_frequency.keys()
max_frequency = 0
​
for (i = 0; i < length(all_kmers); i++)
{
  kmer = all_kmers[i]
  freq = kmer_frequency[kmer]
  tally[freq] = tally[freq] + 1
  
  if (freq > max_frequency)
  {
    max_frequency = freq
  }
}
​
## now print in sorted order
for (i = 1; i <= max_frequency; i++)
{
  if (tally.has_key(i))
  {
    freq = tally[i]
    print_to_file(i, freq)
  }
}
```

#### **Step 4.3** 

Using the output from 4.2, plot the kmer frequency spectrum: x-axis is the kmer frequency, and the y-axis is the number of kmers that occur x times. Make sure to plot both the x and y-axis in log space.

#### **Step 4.4** 

**a)** What percent of the genome is unique, e.g. what percent of the kmers occur 1 time. \
**b)** What percent of the genome is repetitive (occurs more than 1 time). c) What percent occurs more than 1000 times? 

Note: For this analysis, you should separately consider all of the kmers in the genome, e.g. the denominator will be G-k+1. When computing the unique percentage, use the number of unique kmers as the numerator. When computing repetitive percentages, make sure to separately count each instance of a repetitve kmer. For example the string "GCATCATCAT" has kmers: GCA, CAT, ATC, TCA, CAT, ATC, TCA, CAT. Of these 1/8 (12.5%) are unique and 7/8 (87.5%) are repetitive


























































