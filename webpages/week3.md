# Resources for the Variant Calling Assignment

Assignment Date: Friday, Sept. 27, 2024

Due Date: Friday, Oct. 4, 2024

## Lecture Slides

Lecture slides are available through Dropbox [here](https://www.dropbox.com/scl/fi/dukux5o61wfko0smxc48y/20250926_variant_calling.pptx?rlkey=fmu6e18dfbxxai4gv6kvt5toi&st=2zyz7r5v&dl=0).


## Live Coding

Code for binomial simulations, tabulation, and plotting is provided below:

```{r}

library(tidyverse)

n_experiments <- 30
n_tosses_per_experiment <- 5
prob_tails <- 0.5

data.frame(ntails = rbinom(n_experiments, n_tosses_per_experiment, prob_tails)) %>%
  count(ntails) %>%
  ggplot(aes(x = ntails, y = n)) +
    geom_bar(stat = "identity")

```

## Homework Assignment

As always, before you do anything else, create a `week3` directory in your `qbb2025-answers` directory for this assignment.

[Homework assignment](https://bxlab.github.io/cmdb-quantbio/assignments/lab/variant_calling/assignment)
