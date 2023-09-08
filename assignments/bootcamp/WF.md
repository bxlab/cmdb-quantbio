# In-Class Exercise: Wright-Fisher Simulation

The Wright-Fisher model is commonly used to investigate the effect of genetic drift: random fluctuation of alleles from generation to generation caused by a finite population size. In this model, we assume:

* Constant population size
* Everyone reproduces once per generation, at the same time
* No selection
* No mutation
* Random mating

We can think of Wright-Fisher as modeling an allele that is evolutionarily neutral.

Since random effects are the only thing changing the frequency of an allele from generation to generation, we can simulate the change frequency from generation by sampling from the binomial distribution. The binomial states that the probability of observing $j$ alleles in the next generation is:

\begin{align}
&\binom{2N}{j}\;p^j_i\;(1-p_i)^{2N-j} \\[1em]
&\binom{n}{k} = \frac{n!}{ k!( n-k)!}
\end{align}

with $N$  as the population size and  $i$ as the current allele frequency. Note that the quantity we care about here is  $2N$: the number of chromosomes.

The easiest way to draw from the binomial distribution is to use the function `np.random.binomial(n, p)`, where `n` is the number of trials and `p` is the probability of success.

## Exercise 1: The Wright-Fisher Model

In the space below, create a function that implements the Wright-Fisher model. Your function should accept two arguments:

1. a starting allele frequency
2. the population size

Your function should run until one of the two alleles reaches fixation (i.e., your allele frequency hits 0 or 1).

In each generation, your function should take the current allele frequency and generate the allele frequency of the next generation by sampling from the binomial distribution. 

Your function should return a list containing the allele frequencies at every generation, including the first and last generations.
