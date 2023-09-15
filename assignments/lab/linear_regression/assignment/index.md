# Linear Regression

## Assignment Overview

The goal of today's lab is to use linear regression and other statistical methods to investigate the relationship between paternal age, maternal age, and the number of de novo mutations (DNMs) in a proband (offspring). Today's assignment will build familiarity with manipulating tabular datasets containining mixed data types using the `pandas` library in Python. Specifically, you will import a table of de-novo mutations and will manipulate it to calculate the number of maternal and paternal DNMs per individual. This assignment will also introduce the `statsmodels` package, which you will use to perform statistical analysis of your data. 

## Data

Data are taken from [Halldorsson, B. V., Palsson, G., Stefansson, O. A., Jonsson, H., Hardarson, M. T., Eggertsson, H. P., ... & Gudjonsson, S. A. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science, 363(6425)](https://science.sciencemag.org/content/363/6425/eaau1043.abstract).

Read the abstract from the above paper to understand the context of the datasets you will be using. The data you need for this assignment has already been loaded onto your laptop. There are two files we'll be using for this assignment:
1. information about the number and parental origin of each de novo mutation detected in an offspring individual (i.e. "proband"), stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_dnm.csv`
2. ages of the parents of each proband, stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_parental_age.csv`

You can use this data as is, or make copies of it in your submission directory for this assignment. If you do make copies in your submission directory, don't forget to add them to your `.gitignore` file within the submission directory.

Before beginning the assignment, you should examine the two files (with `less -S` perhaps) to make sure you understand how they're organized.<br><br>

## Exercises

### Exercise 1: Wrangle the data

#### **Step 1.1**
You'll start by exploring the data in `aau1043_dnm.csv`. First, load this data into a `pandas` dataframe.

#### **Step 1.2**
You first want to count the number of paternally and maternally inherited DNMs in each proband. Using this dataframe, create a dictionary where the keys are the proband IDs and the value associated with each key is a list of length 2, where the first element in the list is the number of maternally inherited DNMs and the second element in the list is the number of paternally inherited DNMs for that proband. You can ignore DNMs without a specified parent of origin.

#### **Step 1.3**
Use the following code snippet to convert this dictionary into a new pandas dataframe (this assumes your dictionary from step 1.2 is called `deNovoCount`):

`deNovoCountDF = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal_dnm', 'paternal_dnm'])`

Feel free to ask questions about how this code is working or, if you're interested, you can try to figure it out yourself.

#### **Step 1.4**
Now, load the data from `aau1043_parental_age.csv` into a new `pandas` dataframe.

#### **Step 1.5** 
You now have two dataframes with complementary information. It would be nice to have all of this in one data structure. Use the `pd.concat()` function (more [here](https://pandas.pydata.org/docs/reference/api/pandas.concat.html)) to combine your dataframe from step 3 with the dataframe you just created in step 4 to create a new merged dataframe.

**NOTE**: You will need to specify the `axis` and `join` arguments in `pd.concat()`

### Exercise 2: Fit and interpret linear regression models with Python

Using the merged dataframe from the previous section, you will be exploring the relationships between different features of the data. The `statsmodels` package (more [here](https://www.statsmodels.org/stable/index.html)) is an incredibly useful package for conducting statistical tests and running regressions. As such, it is especially appropriate for the types of questions we're interested in here. For this assignment, we'll be using the `formula` api from `statsmodels` (more [here](https://www.statsmodels.org/stable/example_formulas.html)) to run some regressions between variables in our dataset. You can load this tool into Python with `import statsmodels.formula.api as smf`. 

#### **Step 2.1**
First, you're interested in exploring if there's a relationship between the number of DNMs and parental age. Use `matplotlib` to plot the following. **All plots should be clearly labelled and easily interpretable**.
1. the count of maternal de novo mutations vs. maternal age (upload as `ex2_a.png` in your submission directory)
2. the count of paternal de novo mutations vs. paternal age (upload as `ex2_b.png` in your submission directory)

#### **Step 2.2**
Now that you've visualized these relationships, you're curious whether they're statistically significant. Perform ordinary least squares using the `smf.ols()` function to test for an association between *maternal* age and *maternally* inherited de novo mutations. In your `README.md` for this assignment, answer the following questions:
1. What is the "size" of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?
2. Is this relationship significant? How do you know?

#### **Step 2.3**
As before, perform ordinary least squares using the `smf.ols()` function, but this time to test for an association between *paternal* age and *paternally* inherited de novo mutations. In your `README.md` for this assignment, answer the following questions:
1. What is the "size" of this relationship? In your own words, what does this mean? Does this match what you observed in your plots in step 6?
2. Is this relationship significant? How do you know?

#### **Step 2.4**
Using your results from **step 2.3**, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband's time of birth. Record your answer and your work (i.e. how you got to that answer) in your `README.md`.

#### **Step 2.5**
Next, you're curious whether the number of paternally inherited DNMs match the number of maternally inherited DNMs. Using `matplotlib`, plot the distribution of maternal DNMs per proband (as a histogram). In the same panel (i.e. the same `axes`) plot the distribution of paternal DNMs per proband. Make sure to make the histograms semi-transparent so you can see both distributions. Upload as `ex2_c.png` in your submission directory.

#### **Step 2.6**
Now that you've visualized this relationship, you want to test whether there is a *significant* difference between the number of maternally vs. paternally inherited DNMs per proband. What would be an appropriate statistical test to test this relationship? Choose a statistical test, and find a Python package that lets you perform this test. If you're not sure where to look, the `stats` module from `scipy` (more [here](https://docs.scipy.org/doc/scipy/reference/stats.html)) provides tools to perform several different useful statistical tests. After performing your test, answer the following answers in your `README.md` for this assignment:
1. What statistical test did you choose? Why?
2. Was your test result statistically significant? Interpret your result as it relates to the number of paternally and maternally inherited DNMs.

### Exercise 3 (OPTIONAL)
Note that standard linear regression assumes a continuous response variable. When we want to work with response variables that are "counts", such as the number of de novo mutations, we should technically use an approach such as "Poisson regression" that is designed for count data. To fit a Poisson regression model with Python statsmodels, simply use `smf.poisson()` in place of `smf.ols()`.

#### **Step 3.1**
Re-fit the models above (steps 2 and 3 in Exercise 2) using Poisson regression.

#### **Step 3.2**
The interpretation of parameter estimates from Poisson regression differs from that of OLS. Using the relevant Poisson regression model that you fit, predict the number of paternal de novo mutations for a proband with a father who was 40.2 years old at the proband's time of birth. Record your answer and your work (i.e. how you got to that answer) in your `README.md`.

### Exercise 4 (OPTIONAL)

Do the analyses we tell you to do is surely fun, but isn't it more fun to do your *own* analyses?

#### **Step 4.1**
Select a new dataset from those listed at the bottom of [this website](https://github.com/rfordatascience/tidytuesday). The corresponding data can generally be found as a `.csv` file in the `tidytuesday/data/<year>/<date>` subdirectory of the GitHub repository. Record which dataset you picked in your `README.md`.

#### **Step 4.2**
Generate figures to explore these data. What patterns do you notice? Record your observations in your `README.md`, and submit any figures you make (upload as `ex4_<something>.png`).

#### **Step 4.3**
Pose a hypothesis about the data that can be tested with a linear regression model.

Fit your model, evaluate the model fit, and test your hypothesis. Record your hypothesis and results in your `README.md`<br><br>

## Submission

1. Python script with all code for the assignment (**5.5 points total**)
  * Code to create dataframe with paternal and maternal DNMs per proband (**1 point**)
  * Code to generate merged dataframe with DNMs per proband AND parental ages (**0.5 point**)
  * Code to generate `ex2_a.png` plot in Step 2.1 (**0.5 point**)
  * Code to generate `ex2_b.png` plot in Step 2.1 (**0.5 point**)
  * Code to run maternal OLS in Step 2.2 (**1 point**)
  * Code to run paternal OLS in Step 2.3 (**1 point**)
  * Code to run maternal DNM vs paternal DNM statistical test in Step 2.6 (**1 point**)
2. `README.md` file with answers to questions in the assignment (**3.5 points total**)
  * Answer to questions in Step 2.2 (**1 point**)
  * Answer to questions in Step 2.3 (**1 point**)
  * Answer to question in Step 2.4 (**0.5 point**)
  * Answer to questions in Step 2.6 (**1 point**)
4. `ex2_a.png` clearly labelled and easily interpretable (**0.5 point**)
5. `ex2_b.png` clearly labelled and easily interpretable (**0.5 point**)

**Total Points: 10**

<br><br>
