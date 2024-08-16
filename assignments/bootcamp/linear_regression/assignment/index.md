## Assignment Overview

The goal of this exercise is to use linear regression and other statistical methods to investigate the relationship between paternal age, maternal age, and the number of de novo mutations (DNMs) in a proband (offspring). Today's assignment will build familiarity with manipulating tabular datasets containining mixed data types. Specifically, you will import a table of de novo mutations and will manipulate it to calculate the number of maternal and paternal DNMs per individual, which you will then model with linear regression.

## Data

Data are taken from [Halldorsson, B. V., Palsson, G., Stefansson, O. A., Jonsson, H., Hardarson, M. T., Eggertsson, H. P., ... & Gudjonsson, S. A. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science, 363(6425)](https://science.sciencemag.org/content/363/6425/eaau1043.abstract).

Read the abstract from the above paper to understand the context of the datasets you will be using. There are two files we'll be using for this assignment: 1. information about the number and parental origin of each de novo mutation detected in an offspring individual (i.e. "proband"), stored in `aau1043_dnm.csv` 2. ages of the parents of each proband, stored in `aau1043_parental_age.csv`

Before beginning the assignment, you should examine the two files (with `less -S` perhaps) to make sure you understand how they're organized.

## Exercises

### Exercise 1: Wrangle the data

#### **Step 1.1**

You'll start by exploring the data in `aau1043_dnm.csv`. First, load this data into a tibble.

#### **Step 1.2**

Use `group_by()` and `summarize()` to tabulate the number of paternally and maternally inherited DNMs in each proband. Note that the maternal versus paternal origin of the mutations are recorded in the column titled `Phased_combined`.

#### **Step 1.3**

Now, load the data from `aau1043_parental_age.csv` into a new `pandas` dataframe.

#### **Step 1.4**

You now have two dataframes with complementary information. It would be nice to have all of this in one data structure. Use the `left_join()` function to combine your dataframe from step 2 with the dataframe you just created in step 3 based on the shared column `Proband_id`.

### Exercise 2: Fit and interpret linear regression models with Python

Using the merged dataframe from the previous section, you will be exploring the relationships between different features of the data.

#### **Step 2.1**

First, you're interested in exploring if there's a relationship between the number of DNMs and parental age. Use ggplot2 to plot the following. All plots should be clearly labelled and easily interpretable.

1\. the count of maternal de novo mutations vs. maternal age

2\. the count of paternal de novo mutations vs. paternal age

#### **Step 2.2**

Now that you've visualized these relationships, you're curious whether they're statistically significant. Fit a linear regression model to the data using the `lm()` function.

1.  What is the "size" of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 2.1?

There is an average increase of 1.35 paternal DNMs per year of paternal age.

2.  Is this relationship significant? How do you know? In your own words, what does this mean?

The relationship is statistically significant with a p-value of 1.6e-84 which is less than our threshold of 0.05. We interpret this to mean that there is a very small probability of observing a relationship this extreme or more extreme assuming that the null hypothesis is true (i.e., that there is no relationship between paternal age and paternal DNM count).

#### **Step 2.3**

As before, fit a linear regression model, but this time to test for an association between *paternal* age and *paternally* inherited de novo mutations.

1.  What is the "size" of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?

There is an average increase of 0.38 maternal DNMs per year of maternal age.

2.  Is this relationship significant? How do you know? In your own words, what does this mean?

The relationship is statistically significant with a p-value of 6.9e-24 which is less than our threshold of 0.05. We interpret this to mean that there is a very small probability of observing a relationship this extreme or more extreme assuming that the null hypothesis is true (i.e., that there is no relationship between maternal age and maternal DNM count).

#### **Step 2.4**

Using your results from **step 2.3**, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband's time of birth. Record your answer and your work (i.e. how you got to that answer).

#### **Step 2.5**

Next, you're curious whether the number of paternally inherited DNMs match the number of maternally inherited DNMs. Plot the distribution of maternal DNMs per proband (as a histogram). In the same panel (i.e. the same set of axes) plot the distribution of paternal DNMs per proband. Make sure to make the histograms semi-transparent so you can see both distributions.

#### **Step 2.6**

Now that you've visualized this relationship, you want to test whether there is a *significant* difference between the number of maternally vs. paternally inherited DNMs per proband. What would be an appropriate statistical test to test this relationship?

After performing your test, answer the following questions:

1.  What statistical test did you choose? Why?

2.  Was your test result statistically significant? Interpret your result as it relates to the number of paternally and maternally inherited DNMs.

### Exercise 3 (OPTIONAL)

Note that standard linear regression assumes a continuous response variable. When we want to work with response variables that are "counts", such as the number of de novo mutations, we should technically use an approach such as "Poisson regression" that is designed for count data. To fit a Poisson regression model use the `glm()` function with the argument `family = "poisson"`.

#### **Step 3.1**

Re-fit the models above (steps 2 and 3 in Exercise 2) using Poisson regression.

#### **Step 3.2**

The interpretation of parameter estimates from Poisson regression differs from that of ordinary least squares, as the reported coefficients are on the log scale. Predictions can be converted back to the normal response scale through exponentiation: `exp(x)`.

Using the relevant Poisson regression model that you fit, predict the number of paternal de novo mutations for a proband with a father who was 40.2 years old at the proband's time of birth. Record your answer and your work.


### Exercise 4 (OPTIONAL)

Do the analyses we tell you to do is surely fun, but isn't it more fun to do your *own* analyses?

#### **Step 4.1** 

Select a new dataset from those listed at the bottom of [this website](https://github.com/rfordatascience/tidytuesday). The corresponding data can generally be found as a `.csv` file in the `tidytuesday/data/<year>/<date>` subdirectory of the GitHub repository. Record which dataset you picked.

#### **Step 4.2**

Generate figures to explore these data. What patterns do you notice? Record your observations and submit any figures you make.

#### **Step 4.3**

Pose a hypothesis about the data that can be tested with a linear regression model.

Fit your model, evaluate the model fit, and test your hypothesis. Record your hypothesis and results.
