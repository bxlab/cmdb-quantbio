# Assignment 11: The Wright-Fisher model
Assignment Date: Friday, Dec. 2, 2022 <br>
Due Date: Friday, Dec. 9, 2022 @ 1pm ET <br>

## The Wright-Fisher model

Recall, the Wright-Fisher model of DNA sequence evolution:

- A genetic locus with two possible alleles
- Non-overlapping generations
- Fixed population of size \(*N*\) (thus \(2*N*\) chromosomes (observed non-unique alleles) per generation)
- Random mating

Allele frequencies range from 0 to 1 within a population. The frequency of an allele is defined as the fraction of all chromosomes in the population that carry that allele.

Let's say that we have two alleles, A and B. If the A allele's frequency is 0, then the B allele's frequency is 1. If the A allele's frequency is 0.8, then the B allele's frequency is 0.2. An allele is said to be fixed in a population when it has reached a frequency of 1. It is fixed because it is the only allele that will be sampled in the next generation (unless mutation creates a new allele).

Because of these characteristics (two alleles, etc.), the pool of alleles at each generation can be described by a stochastic process with a binomial transition probability. We can determine the ratio of alleles in the next generation by drawing randomly (with replacement) from the previous generation's alleles.

What parameters do we need for randomly drawing from a binomial process?

Numpy allows us to draw samples from a binomial distribution: [`numpy.random.binomial`](https://numpy.org/doc/stable/reference/random/generated/numpy.random.binomial.html) and it has the parameters `n` and `p`

* `n` is the number of trials (or the number of alleles in the generation).
* `p` is the probability of success or we can say the probability of selecting an A allele

But how do we know the values for the parameters we should use for this binomial process?

* We have already described the number of alleles in the generation

* The probability of selecting an A allele: If \(*i*\) is the number of A alleles at generation \(*n* - 1\), then the probability of selecting an A allele from the population is \(*p*\_A = *i* / 2*N*\). (Note: The probability of selecting a B allele from the population is then \(*p*\_B = \(\(2*N*\) - *i*\) / 2*N*\).)

The probability of observing *j* A alleles at a generation \(*n*\), assuming no selection, and given *i* A alleles at generation \(*n-1*\) is described by the following statement:

![\Large binomial](https://latex.codecogs.com/svg.latex?P%28%20X_n%20%7C%20X_%7Bn-1%7D%20%3D%20i%20%29%20%3D%20%5Cbinom%7B%202N%20%7D%7B%20j%20%7D%20p_i%5Ej%20%28%201%20-%20%7Bp_i%7D%20%29%5E%7B2N-j%7D%20%5Csim%20Binomial%28%20n%3D2N%2C%20p%3D%5Cfrac%7Bi%7D%7B2N%7D%29)

(Note *i* may or may not equal *j*)

By sampling with `numpy.random.binomial`, you are applying the previous equation without having to explicitly code the math.

We can introduce selection (strength determined by the variable *s*) by changing the probability of sampling the A allele from a generation with *i* A alleles. For example:

![\Large p_A = \frac{ i ( 1 + s ) }{ 2N - i + i ( 1 + s )}](https://latex.codecogs.com/svg.latex?%5Cdpi%7B150%7D%20%5Clarge%20p_A%20%3D%20%5Cfrac%7B%20i%20%28%201%20&plus;%20s%20%29%20%7D%7B%202N%20-%20i%20&plus;%20i%20%28%201%20&plus;%20s%20%29%7D)

So to introduce selection into a Wright-Fisher simulation, you will have to change the `p` parameter using the \(p\_A\) equation provided above.

## Assignment

1. Implement a Wright-Fisher simulation of allele frequencies for an arbitrary starting allele frequency and population size. The simulation should run until one of the alleles becomes fixed (reaches a frequency of 1). Implement this simulation as a python function that takes two input arguments: the starting allele frequency, and the population size. As an output, the function should return a list that contains the allele frequency at each generation.

	HINT: You definitely want to do this in Python. Remember that you are randomly drawing allele samples from a binomial distribution using `numpy.random.binomial`.

2. Write a function that plots a line showing the allele frequency (y-axis) versus generation number (x-axis) for the entirety of your simulation. Produce such a plot for one simulation. For that simulation, use any starting allele frequency and population size you would like, but notate these values on the plot. Make sure you label your axes.

3. For a starting allele frequency of 0.5, and a population size of 100, produce a histogram with density showing time to fixation over (at least) 1000 independent runs of the simulation.

    HINT: The time to fixation is the number of generations it takes for an allele to fix. Since the function from part 1 is returning a list that contains the allele frequency at each generation, and your function runs until one of the alleles is fixed, how can you find the time to fixation?

4. For a starting allele frequency of 0.5, vary the population size and produce a line plot that shows fixation time (y-axis) vs \(*N*\) (x-axis). A reasonable range of population sizes is 100 to 10 million. Use at least 6 different population sizes. Be sure to label your ticks with the population sizes you used (and maybe consider setting a log scale for the x-axis). Note: simulations with large population sizes will take a while to run.

5. Run the simulation under a range of different starting allele frequencies, but a constant population size. Use 10 different starting allele frequencies ranging from 0 to 1. Produce a plot showing starting allele frequency (x-axis) vs. time to fix (y-axis). Do (at least) 100 simulations for each starting allele frequency and include the variability in your plot. To show variability in the time to fixation, use something like seaborn stripplot ([documentation](https://seaborn.pydata.org/generated/seaborn.stripplot.html), [example](https://towardsdatascience.com/jitter-plots-with-pythons-seaborn-62188bf511b8)) to make a jitter plot, a violin plot ([documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.violinplot.html)), or a boxplot ([documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html)). Make sure to label your ticks with the starting allele frequencies used and notate on the plot what constant population size was used.

## Advanced Exercise

1. Introduce selection to your function from Part 1 (as an additional parameter that can be specified) and plot the allele frequency trajectory for some chosen parameters. On your plot, make sure you note what your selection coefficient for the simulation was. Additionally, plot selection coefficient vs time to fixation for a fixed population size of your choice. On your plot, make sure you note what your population size was.

## Submitting your assignment
Submit your script and plots
