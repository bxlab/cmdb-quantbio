# Genome Assembly

## Assignment Overview

***FILL IN***

## Exercises 

### Exercise 1: Coverage simulator 

#### **Step 1.1**
In your `README.md` for this assignment, answer the following question (show your work):
1. How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?

#### **Step 1.2**
Using Python, simulate sequencing 3x coverage of a 1Mbp genome with 100bp reads. Note that you do not need to actually simulate the sequences of the reads, you can just randomly sample positions in the genome and record the coverage. You do not need to consider the strand of each read.

The start position of each read should have a uniform random probabilty at each possible starting position (0 through 999,900). You can record the coverage in an array of 1M positions.

Now, plot the histogram of coverage across the genome. Overlay the histogram with a Poisson distribution with **lambda = 3**. Also overlay the distribution with a Normal distribution with a **mean of 3 and a std. dev. of 1.73** (which is the square root of 3).
* **HINT**: For the poisson and normal distributions, you’ll need to find the probability of getting a certain coverage. For the poisson distribution, this is called the probability mass function (PMF) of the distribution. Feel free to code this yourself using the appropriate equation, or you can take a look at the `scipy.stats.poisson.pmf()` function (more [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.poisson.html)). Unlike the Poisson distribution, the normal distribution is continuous, and so we will instead want to use the probability density function (PDF). Note that for both of these, this will give you the *probability* of observing each coverage. What do we need to do to transform these probabilities into a frequency count comparable to those in our histogram?

Upload this plot as `ex1_3x_cov.png` in your submission directory. **All plots should be clearly labelled and easily interpretable** (i.e. axis labels, legend describing the three things plotted, etc.).

<details><summary><b><font color="#18BC9C">CLICK HERE FOR PSEUDOCODE</font></b></summary>
  <pre>
    <code>
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
    </code>
  </pre>
</details>

#### **Step 1.3**
Using your results from Step 1.2, answer the following questions in your `README.md`:
1. In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
2. How well does this match Poisson expectations? How well does the normal distribution fit the data?

#### **Step 1.4**
Now, repeat the analysis with 10x coverage:
1. Simulate using the appropriate number of reads
2. Make a histogram. Overlay a Poisson distribution with **lambda = 10**. Overlay a Normal distribution with **mean = 10 and std. dev. = 3.16**. Upload this plot as `ex1_10x_cov.png` in your submission directory.
3. In your `README.md`, answer the following questions:
   1. In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
   2. How well does this match Poisson expectations? How well does the normal distribution fit the data?

#### **Step 1.5**
Now, repeat the analysis with 30x coverage:
1. Simulate using the appropriate number of reads
2. Make a histogram. Overlay a Poisson distribution with **lambda = 30**. Overlay a Normal distribution with **mean = 30 and std. dev. = 5.47**. Upload this plot as `ex1_30x_cov.png` in your submission directory.
4. In your `README.md`, answer the following questions:
   1. In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
   2. How well does this match Poisson expectations? How well does the normal distribution fit the data?

### Exercise 2: De Bruijn graph construction

#### **Step 2.1**

Next, you're going to generate your own de Bruijn graph using a provided set of reads. Copy the list of reads below into your code:

```
reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']
```

Write code to find all of the edges in the de Bruijn graph corresponding to the provided reads using **k = 3** (assume all reads are from the forward strand, no sequencing errors, complete coverage of the genome). Each edge should be of the format `ATT -> TTC`. Write all edges to a file, with each edge as its own line in the file.

<details><summary><b><font color="#18BC9C">CLICK HERE FOR PSEUDOCODE</font></b></summary>
  <pre>
    <code>
graph = {}

for each read
  for (i = 0; i < read.len() - k - 1; i++)
     kmer1 = read.substr(i, k)
     kmer2 = read.substr(i+1, k)
     graph[kmer1][kmer2] = 1

for each kmer1 in graph
   for each kmer2 in graph[kmer1]
      print "kmer1 -> kmer2"
    </code>
  </pre>
</details>

#### **Step 2.2**

Now that we know all of the edges, we can go about actually visualizing the de Bruijn graph. For this task, we'll be using the [graphviz](https://graphviz.org/) command line tool.

Create a `conda` environment for this assignment and install `graphviz`:

```
conda create -n graphviz -c conda-forge graphviz
conda activate graphviz
```

#### **Step 2.3**

Graphviz's command line tool is called `dot`. Read more about how to use `dot` and the file format it's expecting [here](https://graphviz.org/doc/info/command.html).
* **NOTE**: We're going to want to produce a *directed* graph.

Based on what you've read, modify your code from Step 2.1 to output the edges in a format `dot` can use.

#### **Step 2.4**

Now, use `dot` to produce a directed graph. Upload this graph as `ex2_digraph.png` in your submission directory. You do NOT need to upload the text file of edges you used to make the graph.

#### **Step 2.5**

Assume that the maximum number of occurrences of any 3-mer in the actual genome is 4. Using your graph from Step 2.4,  write one possible genome sequence that would produce these reads. Record your answer in your `README.md`.

#### **Step 2.6**

In a few sentences, what would it take to accurately reconstruct the sequence of the genome? Record your answer in your `README.md`.

### Exercise 3: Why genomics?

#### **Step 3.1**

Use ChatGPT (or Bard or your favorite LLM) to write an essay on why you are interested in genomics. Make sure to ask for references. Record both your prompt(s) and the output from the LLMin your `README.md`.
​
#### **Step 3.2**

In your `README.md`, comment on the output from the LLM: Does it make logical sense? Does it include any phrases you would not have written? Do the cited papers exist and support the claims from the LLM?


### Exercise 4: K-mer uniqueness (OPTIONAL)

Download the human chomosome 22 DNA sequence using the following command:

```
wget https://schatz-lab.org/appliedgenomics2023/assignments/assignment1/chr22.fa.gz
```
#### **Background**

A kmer is a substring of length k. For example, the string GATTACA, has these 3-mers: GAT, ATT, TTA, TAC, ACA.

A string of length G has G - k + 1 kmers. For long strings, G - k + 1 is nearly the same as G e.g. for human using 19mers, 3,000,000,000 vs 2,999,999,986.

While a string of length G has G-k+1 kmers, there may be many fewer *distinct* kmers. For example, in the string "GCATCATCATCATCATCATCAT..." the kmers are: GCA, CAT, ATC, TCA, CAT, ATC, TCA, CAT, ATC, TCA, CAT, .... As such, there are only 4 disinct kmers (GCA, CAT, ATC, TCA). Of these, GCA occurs once and the others occur many times.

#### **Step 4.1** 

How many As, Cs, Gs, Ts and Ns are found in the entire chromosome? If needed, convert lowercase letters to uppercase. Any other character can be converted to N. Record your answer in your `README.md`.

#### **Step 4.2**

Using Python, tally the frequency of all of the different 19-mers in the chromosome, and calculate the kmer frequency spectrum up to 1000 (e.g. how many kmers occur 1 time, how many occur 2 times, how many occur 3 times, etc.). Store this in a list. For this, convert lowercase letters to uppercase, and convert any character that is not A, C, G or T to A (especially N characters). We recommend you use a dictionary to tally the frequencies using this pseudocode.

In your `README.md`, show the kmer frequency spectrum for 1 to 20, e.g. how many kmers occur 1 time, how many occur 2, ..., how many occur 20 times.

<details><summary><b><font color="#18BC9C">CLICK HERE FOR PSEUDOCODE</font></b></summary>
  <pre>
    <code>
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
    </code>
  </pre>
</details>

#### **Step 4.3** 

Using your list from Step 4.2, plot the kmer frequency spectrum: x-axis is the kmer frequency, and the y-axis is the number of kmers that occur at each of those frequencies times. Make sure to plot both the x and y-axis in log space. Upload this graph as `ex4_kmer_spec.png` in your submission directory.

#### **Step 4.4** 

Using your list from Step 4.2, answer the following questions. Record your answers in your `README.md`.
1. What percent of the genome is unique (e.g. what percent of the kmers occur 1 time)?
2. What percent of the genome is repetitive (occurs more than 1 time)?
3. What percent occurs more than 1000 times? 

* **Note**: For this analysis, you should separately consider all of the kmers in the genome, e.g. the denominator will be G-k+1. When computing the unique percentage, use the number of unique kmers as the numerator. When computing repetitive percentages, make sure to separately count each instance of a repetitve kmer. For example the string "GCATCATCAT" has kmers: GCA, CAT, ATC, TCA, CAT, ATC, TCA, CAT. Of these 1/8 (12.5%) are unique and 7/8 (87.5%) are repetitive

<br><br>
