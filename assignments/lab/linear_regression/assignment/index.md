# Linear Regression

## Assignment Overview

<!-- ADD THIS (use DNM acronym) --> 

## Data

Data are taken from Halldorsson, B. V., Palsson, G., Stefansson, O. A., Jonsson, H., Hardarson, M. T., Eggertsson, H. P., ... & Gudjonsson, S. A. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science, 363(6425). [link](https://science.sciencemag.org/content/363/6425/eaau1043.abstract)

Read the abstract from the above paper to understand the context of the datasets we will be using. The data you need for this assignment has already been loaded onto your laptop in two files:
1. information about the number and parental origin of each de novo mutation detected in an offspring individual (i.e. "proband"), stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_dnm.csv`
2. ages of the parents of each proband, stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_parental_age.csv`

You can use this data as is, or make copies of it in your submission directory for this assignment. If you do make copies in your submission directory, don't forget to add them to your `.gitignore` file within the submission directory.

Before beginning the assignment, you should examine the two files (with `less -S` perhaps) to make sure you understand how they're organized.

## Exercises

### Exercise 1: Wrangle the data

1. We'll start by exploring the data in `aau1043_dnm.csv`. First, load this data into a `pandas` dataframe. Using this dataframe, we want to count the number of paternally and maternally inherited DNMs in each proband.
2. Using this dataframe, create a dictionary where the keys are the proband IDs and the associated value with each key is a list of length 2, where the first element in the list is the number of maternally inherited DNMs and the second element in the list is the number of paternally inherited DNMs for that proband. You can ignore DNMs without a specified parent of origin.
3. Use the follwing code snippet to convert this dictionary into a new pandas dataframe (this assumes your dictionary from step 2 is called `deNovoCount`):
`deNovoCountDf = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal', 'paternal'])`
Feel free to ask questions about how this code is working, or if you're interested, you can try to figure it out yourself.
4. Now, load the data from `aau1043_parental_age.csv` into a new `pandas` dataframe.
5. You now have two dataframes with complementary information. It would be nice to have all of this in one data strucutre. Use the `pd.concat()` function to combine your dataframe from step 3 with the dataframe you just created in step 4 to create a new merged dataframe.
  * You will need to specify the `axis` and `join` arguments in `pd.concat()`

### Exercise 2: Fit and interpret linear regression models with Python

4. Use numpy `genfromtxt` to load the "joined" data from step 3 into a numpy array. Use the `Names =` option to give your fields informative names.

5. Use matplotlib to plot:
* the count of maternal de novo mutations vs. maternal age (upload as `ex2_a.png`)
* the count of paternal de novo mutations vs. paternal age (upload as `ex2_b.png`)

6. Use ordinary least squares `smf.ols()` to test for an association between maternal age and maternally inherited de novo mutations.
* Is this relationship significant?
* What is the size of this relationship?

7. Use ordinary least squares `smf.ols()` to test for an association between paternal age and paternally inherited de novo mutations.
* Is this relationship significant?
* What is the size of this relationship?

8. Plot a histogram of the number of maternal de novo mutations and paternal de novo mutations per proband on a single plot with semi-transparency (and upload as `ex2_c.png`).

9. Test whether the number of maternally inherited de novo mutations per proband is significantly different than the number of paternally inherited de novo mutations per proband.

10. Predict the number of paternal de novo mutations for a proband with a father who was 50.5 years old at the proband's time of birth.

### Optional exercise

Note that standard linear regression assumes a continuous response variable. When we want to work with response variables that are "counts", such as the number of de novo mutations, we should technically use an approach such as "Poisson regression" that is designed for count data. To fit a Poisson regression model with Python statsmodels, simply use `smf.poisson()` in place of `smf.ols()`.

11. Re-fit the models (questions 6, 7, and 9) above using Poisson regression.

12. The interpretation of parameter estimates from Poisson regression differs from that of OLS. Using the relevant Poisson regression model that you fit, predict the number of paternal de novo mutations for a proband with a father who was 40.2 years old at the proband's time of birth.

### Optional exercise

13. Select a new dataset from those listed at the bottom of this website: https://github.com/rfordatascience/tidytuesday

If not obvious, the corresponding data can generally be found as a `.csv` file in the `tidytuesday/data/<year>/<date>` subdirectory of the GitHub repository.
  
14. Generate figures to explore these data. What patterns do you notice?

15. Pose a question about the data that can be tested with a linear regression model.

16. Fit your model, evaluate the model fit, and test your hypothesis with the `.summary()` method.
