# Linear Regression

## Assignment Overview

<!-- ADD THIS (use DNM acronym) --> 

## Data

Data are taken from [Halldorsson, B. V., Palsson, G., Stefansson, O. A., Jonsson, H., Hardarson, M. T., Eggertsson, H. P., ... & Gudjonsson, S. A. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science, 363(6425)](https://science.sciencemag.org/content/363/6425/eaau1043.abstract).

Read the abstract from the above paper to understand the context of the datasets we will be using. The data you need for this assignment has already been loaded onto your laptop in two files:
1. information about the number and parental origin of each de novo mutation detected in an offspring individual (i.e. "proband"), stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_dnm.csv`
2. ages of the parents of each proband, stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_parental_age.csv`

You can use this data as is, or make copies of it in your submission directory for this assignment. If you do make copies in your submission directory, don't forget to add them to your `.gitignore` file within the submission directory.

Before beginning the assignment, you should examine the two files (with `less -S` perhaps) to make sure you understand how they're organized.

## Exercises

### Exercise 1: Wrangle the data

1. We'll start by exploring the data in `aau1043_dnm.csv`. First, load this data into a `pandas` dataframe. Using this dataframe, we want to count the number of paternally and maternally inherited DNMs in each proband.
2. Using this dataframe, create a dictionary where the keys are the proband IDs and the value associated with each key is a list of length 2, where the first element in the list is the number of maternally inherited DNMs and the second element in the list is the number of paternally inherited DNMs for that proband. You can ignore DNMs without a specified parent of origin.
3. Use the follwing code snippet to convert this dictionary into a new pandas dataframe (this assumes your dictionary from step 2 is called `deNovoCount`):
`deNovoCountDf = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal', 'paternal'])`
Feel free to ask questions about how this code is working, or if you're interested, you can try to figure it out yourself.
4. Now, load the data from `aau1043_parental_age.csv` into a new `pandas` dataframe.
5. You now have two dataframes with complementary information. It would be nice to have all of this in one data strucutre. Use the [`pd.concat()`](https://pandas.pydata.org/docs/reference/api/pandas.concat.html) function to combine your dataframe from step 3 with the dataframe you just created in step 4 to create a new merged dataframe.
  * You will need to specify the `axis` and `join` arguments in `pd.concat()`

### Exercise 2: Fit and interpret linear regression models with Python

Using the merged dataframe from the previous section, you will be exploring the relationships between different features of the data. [`statsmodels`](https://www.statsmodels.org/stable/index.html) is a Python package for conducting statistical tests as well as linear regressions. As such, it is especially appropriate for the types of questions we're interested in here. For this assignment, we'll be using the `formula` api from `statsmodels` to run some regressions between variables in our dataset (read more [here](https://www.statsmodels.org/stable/example_formulas.html)). You can load this tool into Python with `import statsmodels.formula.api as smf`. 

6. First, you're interested in exploring if there's a relationship between the number of DNMs and parental age. Use `matplotlib` to plot the following. All plots should be clearly labelled and easily interpretable.
 * the count of maternal de novo mutations vs. maternal age (upload as `ex2_a.png` in your submission directory)
 * the count of paternal de novo mutations vs. paternal age (upload as `ex2_b.png` in your submission directory)

7. Now that you've visualized these relationships, you're curious whether they're statistically significant. Perform ordinary least squares using the `smf.ols()` function to test for an association between *maternal* age and *maternally* inherited de novo mutations. In your `README.md` for this assignment, answer the following questions:
 * What is the "size" of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?
 * Is this relationship significant? How do you know?

8. As before, perform ordinary least squares using the `smf.ols()` function, but this time to test for an association between *paternal* age and *paternally* inherited de novo mutations. In your `README.md` for this assignment, answer the following questions:
 * What is the "size" of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?
 * Is this relationship significant? How do you know?

9. Using the results of step 8, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband's time of birth. Record your answer and your work (i.e. how you got to that answer) in your `README.md`.

10. Next, you're curious whether the number of paternally inherited DNMs match the number of maternally inherited DNMs. Using `matplotlib`, plot the distribution of maternal DNMs per proband (as a histogram). In the same panel (i.e. the same `axes`) plot the distribution of paternal DNMs per proband. Make sure to make the histograms semi-transparent so you can see both distributions. Upload as `ex2_c.png` in your submission directory.

11. Now that you've visualized this relationship, you want to test whether there is a *significant* difference between the number of maternally vs. paternally inherited DNMs per proband. What would be an appropriate statistical test to test this relationship? Choose a statistical test, and find a Python package that lets you perform this test. If you're not sure where to look, [the `stats` module from `scipy`](https://docs.scipy.org/doc/scipy/reference/stats.html) provides tools to perform several different useful statistical tests. After performing your test, answer the following answers in your `README.md` for this assignment:
 * What statistical test did you choose? Why?
 * Was your test result statistically significant? Interpret your result as it relates to the number of paternally and maternally inherited DNMs.


### Optional Exercise 3

Note that standard linear regression assumes a continuous response variable. When we want to work with response variables that are "counts", such as the number of de novo mutations, we should technically use an approach such as "Poisson regression" that is designed for count data. To fit a Poisson regression model with Python statsmodels, simply use `smf.poisson()` in place of `smf.ols()`.

12. Re-fit the models (steps 7 and 8) above using Poisson regression.

13. The interpretation of parameter estimates from Poisson regression differs from that of OLS. Using the relevant Poisson regression model that you fit, predict the number of paternal de novo mutations for a proband with a father who was 40.2 years old at the proband's time of birth. Record your answer and your work (i.e. how you got to that answer) in your `README.md`.

### Optional Exercise 4

14. Select a new dataset from those listed at the bottom of this website: https://github.com/rfordatascience/tidytuesday. If not obvious, the corresponding data can generally be found as a `.csv` file in the `tidytuesday/data/<year>/<date>` subdirectory of the GitHub repository. Record which dataset you picked in your `README.md`.
  
15. Generate figures to explore these data. What patterns do you notice? Record your observations in your `README.md`.

16. Pose a question about the data that can be tested with a linear regression model.

17. Fit your model, evaluate the model fit, and test your hypothesis. Record your hypothesis and results in your `README.md`
