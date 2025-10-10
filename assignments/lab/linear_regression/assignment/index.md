# Linear Regression

**Deadline:** Monday, October 20\
**Resubmission deadline:** Friday, October 31

## Assignment Overview

The goal of today's lab is to use linear regression and related statistical methods to investigate the relationship between paternal age, maternal age, and the number of de novo mutations (DNMs) in a proband (offspring). Today's assignment will build familiarity with manipulating tabular datasets containing mixed data types using the **tidyverse** in **R**. Specifically, you will import a table of de novo mutations and manipulate it to calculate the number of maternal and paternal DNMs per individual. You will then fit and interpret linear models with **stats::lm** (and optionally **stats::glm** for Poisson regression), and tidy results with **broom**.

## Data

Data are taken from [Halldorsson, B. V. et al. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. *Science*, 363(6425)](https://science.sciencemag.org/content/363/6425/eaau1043.abstract).

Read the abstract from the above paper to understand the context of the datasets you will be using. The data you need for this assignment are available from Dropbox at:

1.  Information about the number and parental origin of each de novo mutation detected in a proband (offspring): [aau1043_dnm.csv](https://www.dropbox.com/scl/fi/6e28a3dow872fi02cp537/aau1043_dnm.csv?rlkey=l3gs7fb6igff4el5ov6wg96ai&dl=0)
2.  Ages of the parents of each proband: [aau1043_parental_age.csv](https://www.dropbox.com/scl/fi/sjrq1x1g30h0j10ktxysi/aau1043_parental_age.csv?rlkey=e9g3m9iq4tfsm9ski7w4vb0bf&dl=0)

You may copy these into your submission directory (and add to your `.gitignore`).

Before beginning the assignment, take a quick look at both files (e.g., with `less -S` in Unix) to confirm their structure.

------------------------------------------------------------------------

## Getting started (R packages)

Load the tidyverse and broom packages.

```r
```

------------------------------------------------------------------------

## Exercises

### Exercise 1: Wrangle the data

#### **Step 1.1 — Load DNMs**

Load `aau1043_dnm.csv` into a tibble.

```r
```

#### **Step 1.2 — Count DNMs by parental origin per proband**

Create a **per-proband** summary with counts of maternally and paternally inherited DNMs. Ignore DNMs without a specified parent of origin.

```r
```

#### **Step 1.3 — Load parental ages**

Load `aau1043_parental_age.csv`.

```r
```

#### **Step 1.4 — Merge counts with ages**

Join the two tibbles by proband ID.

```r
```

------------------------------------------------------------------------

### Exercise 2: Fit and interpret linear regression models with R

Use your merged data frame for the following. All plots should be clearly labeled and easily interpretable.

#### **Step 2.1 — Visualize relationships**

1)  Create a scatter plot of the count of maternal DNMs vs. maternal age → save as `ex2_a.png`\
2)  Create a scatter plot of the count of paternal DNMs vs. paternal age → save as `ex2_b.png`

```r
```

#### **Step 2.2 — OLS: maternal age vs. maternal DNMs**

Fit a simple linear regression model relating maternal age to the number of maternal de novo mutations. 

In `README.md`, answer: 1. What is the "size" (i.e., slope) of this relationship? Interpret the slope in plain language. Does it match your plot? 2. Is the relationship significant? How do you know? Explain the p-value in plain but precise language.


```r
```

#### **Step 2.3 — OLS: paternal age vs. paternal DNMs**

Repeat the step above but for paternal age vs. paternal DNMs.

```r
```

#### **Step 2.4 — Predict for a 50.5-year-old father**

Use the paternal regression model to predict the expected number of paternal DNMs for a father of age 50.5. You are welcome to do this manually or using a built-in function, but show your work in `README.md`.

```r
```

#### **Step 2.5 — Compare distributions of maternal vs. paternal DNMs**

Plot both distributions on the **same axes** as semi-transparent histograms; save as `ex2_c.png`.

```r
```

#### **Step 2.6 — Statistical test: maternal vs. paternal DNMs per proband**

We have **paired** observations per proband (maternal vs. paternal). The paired t-test assumes that the within-pair differences are approximately normally distributed.

Apply a paired t-test in R using `t.test(merged$maternal_dnm, merged$paternal_dnm, paired = TRUE)`.

In `README.md`, answer: 1. What is the "size" of this relationship (i.e., the average difference in counts of maternal and paternal DNMs)? Interpret the difference in plain language. Does it match your plot? 2. Is the relationship significant? How do you know? Explain the p-value in plain but precise language.

Note that the paired t-test is equivalent to using the difference between the maternal and paternal DNM counts per proband as the response variable and fitting a model with only an intercept term (indicated with `1` on the right side of the model formula). Fit this model using `lm()` and compare to the results of the paired t-test. How would you interpret the coefficient estimate for the intercept term? 

------------------------------------------------------------------------

### Exercise 3: Explore a new dataset

#### **Step 3.1 — Pick a TidyTuesday dataset**

Choose a dataset from the bottom of the [TidyTuesday README](https://github.com/rfordatascience/tidytuesday). Record which one you chose in `README.md`.

#### **Step 3.2 — Explore and visualize**

Generate figures and note any interesting patterns in `README.md`; save figures as `ex4_<something>.png`.

#### **Step 3.3 — Pose and test a linear-model hypothesis**

State a hypothesis, fit a linear model, evaluate fit, and report results in `README.md`.

------------------------------------------------------------------------

## Submission
- R script or R markdown file with analysis code.
- `README.md` file with answers to questions.

## Grading Rubric (Total = 10 points)

### **Exercise 1 — Wrangle the Data (2 points)**
- Load and inspect DNM data (**0.5 pt**)
- Create per-proband maternal and paternal DNM counts (**0.5 pt**)
- Load and inspect parental age data (**0.5 pt**)
- Join counts with ages into a merged table (**0.5 pt**)

### **Exercise 2 — Fit and Interpret Linear Models (6 points)**
- Step 2.1: Scatter plots for maternal and paternal DNMs vs. parental age (**1 pt**)
- Step 2.2: Fit and interpret maternal OLS model (**1 pt**)
- Step 2.3: Fit and interpret paternal OLS model (**1 pt**)
- Step 2.4: Predict paternal DNMs for age 50.5 (**0.5 pt**)
- Step 2.5: Plot distributions of maternal vs. paternal DNMs (**1 pt**)
- Step 2.6: Paired t-test (`t.test` or `lm(diff ~ 1)`) and interpret results (**1.5 pts**)

### **Exercise 3 — Explore a New Dataset (2 points)**
- Step 3.1: Choose and document TidyTuesday dataset (**0.5 pt**)
- Step 3.2: Produce exploratory figure(s) (**0.5 pt**)
- Step 3.3: Pose and test a linear-model hypothesis and interpret results (**1 pt**)

**Total Points: 10**

------------------------------------------------------------------------
