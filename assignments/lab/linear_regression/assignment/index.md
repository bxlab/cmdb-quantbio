# Linear Regression

**Deadline:** Friday, October 9\
**Resubmission deadline:** Friday, October 30

## Assignment Overview

The goal of today's lab is to use linear regression and related statistical methods to investigate the relationship between paternal age, maternal age, and the number of de novo mutations (DNMs) in a proband (offspring). Today's assignment will build familiarity with manipulating tabular datasets containing mixed data types using the **tidyverse** in **R**. Specifically, you will import a table of de novo mutations and manipulate it to calculate the number of maternal and paternal DNMs per individual. You will then fit and interpret simple and multiple linear regression models with **stats::lm**, compare nested models with **stats::anova**, and tidy results with **broom**.

This assignment is an R Markdown notebook. Write your code in the empty code chunks and your answers where you see *Your answer:*.

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

```{r}
```

------------------------------------------------------------------------

## Exercises

### Exercise 1: Wrangle the data

#### **Step 1.1 — Load DNMs**

Load `aau1043_dnm.csv` into a tibble.

```{r}
```

#### **Step 1.2 — Count DNMs by parental origin per proband**

Create a **per-proband** summary with counts of maternally and paternally inherited DNMs. Ignore DNMs without a specified parent of origin. Name the count columns `maternal_dnm` and `paternal_dnm`, as later steps use these names.

```{r}
```

#### **Step 1.3 — Load parental ages**

Load `aau1043_parental_age.csv`.

```{r}
```

#### **Step 1.4 — Merge counts with ages**

Join the two tibbles by proband ID. Name the result `merged`, as later steps use this name.

```{r}
```

------------------------------------------------------------------------

### Exercise 2: Fit and interpret linear regression models with R

Use your merged data frame for the following. All plots should be clearly labeled and easily interpretable.

#### **Step 2.1 — Visualize relationships**

**2.1.1** Create a scatter plot of the count of maternal DNMs vs. maternal age → save as `ex2_a.png`

```{r}
```

**2.1.2** Create a scatter plot of the count of paternal DNMs vs. paternal age → save as `ex2_b.png`

```{r}
```

**2.1.3** Create a scatter plot of paternal age vs. maternal age → save as `ex2_c.png`

```{r}
```

#### **Step 2.2 — OLS: maternal age vs. maternal DNMs**

Fit a simple linear regression model relating maternal age to the number of maternal de novo mutations.

```{r}
```

Answer the following questions:

1.  What is the "size" (i.e., slope) of this relationship? Interpret the slope in plain language. Does it match your plot?

    *Your answer:*

2.  Is the relationship significant? How do you know? Explain the p-value in plain but precise language.

    *Your answer:*

#### **Step 2.3 — OLS: paternal age vs. paternal DNMs**

Repeat the step above but for paternal age vs. paternal DNMs.

```{r}
```

Answer the following questions:

1.  What is the "size" (i.e., slope) of this relationship? Interpret the slope in plain language. Does it match your plot?

    *Your answer:*

2.  Is the relationship significant? How do you know? Explain the p-value in plain but precise language.

    *Your answer:*

#### **Step 2.4 — Predict for a 50.5-year-old father**

Use the paternal regression model to predict the expected number of paternal DNMs for a father of age 50.5. You are welcome to do this manually or using a built-in function, but show your work in the code chunk.

```{r}
```

*Your answer:*

#### **Step 2.5 — Multiple regression: correlated predictors**

Maternal DNMs arose in the mother's germline, so there is no obvious reason that the father's age should affect them.

**2.5.1** Fit a simple linear regression model relating **paternal** age to the number of **maternal** DNMs.

```{r}
```

1.  Is paternal age a significant predictor of maternal DNMs on its own? Using your plot from 2.1.3, explain why this might happen.

    *Your answer:*

**2.5.2** Fit a multiple linear regression model with both maternal age and paternal age as predictors of maternal DNMs, just as we added year as a covariate in the penguin example.

```{r}
```

2.  What happened to the paternal age coefficient and its p-value?

    *Your answer:*

3.  Interpret the maternal age coefficient in plain language. How does its meaning differ from the slope you estimated in Step 2.2? (Hint: what is being held constant?)

    *Your answer:*

**2.5.3** Use `anova()` to compare your two-predictor model from 2.5.2 to your model from Step 2.2 (maternal age only).

> **Hint:** In the live coding, we used `anova()` to test whether island improved a model that already contained sex and year. The same approach works for any pair of models where the smaller model is the larger model with some terms removed (these are called "nested" models). Here, the smaller model is the one that leaves out paternal age.

```{r}
```

4.  Does adding paternal age improve the model? Compare the p-value from `anova()` to the p-value for paternal age in the `summary()` of the two-predictor model. Why are they the same here, when for island we needed `anova()` to get a single answer? (Hint: how many coefficients did island add to the model, and how many does paternal age add?)

    *Your answer:*

#### **Step 2.6 — Compare distributions of maternal vs. paternal DNMs**

Plot both distributions on the **same axes** as semi-transparent histograms; save as `ex2_d.png`.

```{r}
```

#### **Step 2.7 — Statistical test: maternal vs. paternal DNMs per proband**

We have **paired** observations per proband (maternal vs. paternal). The paired t-test assumes that the within-pair differences are approximately normally distributed.

**2.7.1** Apply a paired t-test in R using `t.test(merged$maternal_dnm, merged$paternal_dnm, paired = TRUE)`.

```{r}
```

1.  What is the "size" of this relationship (i.e., the average difference in counts of maternal and paternal DNMs)? Interpret the difference in plain language. Does it match your plot?

    *Your answer:*

2.  Is the relationship significant? How do you know? Explain the p-value in plain but precise language.

    *Your answer:*

**2.7.2** The paired t-test is equivalent to using the difference between the maternal and paternal DNM counts per proband as the response variable and fitting a model with only an intercept term (indicated with `1` on the right side of the model formula). Fit this model using `lm()`.

```{r}
```

3.  How do the results compare to the paired t-test? How would you interpret the coefficient estimate for the intercept term?

    *Your answer:*

------------------------------------------------------------------------

### Exercise 3: Explore a new dataset

#### **Step 3.1 — Pick a TidyTuesday dataset**

Choose a dataset from the bottom of the [TidyTuesday README](https://github.com/rfordatascience/tidytuesday) and load it.

```{r}
```

1.  Which dataset did you choose?

    *Your answer:*

#### **Step 3.2 — Explore and visualize**

Generate figures; save them as `ex3_<something>.png`.

```{r}
```

1.  What interesting patterns do you notice?

    *Your answer:*

#### **Step 3.3 — Pose and test a linear-model hypothesis**

1.  State a hypothesis that you can test with a linear model.

    *Your answer:*

Fit a linear model with at least two predictors to test your hypothesis, and evaluate its fit.

```{r}
```

2.  Why did you include each predictor?

    *Your answer:*

3.  Report and interpret your results. How well does the model fit?

    *Your answer:*

------------------------------------------------------------------------

## Submission

- Your completed R Markdown notebook (`.Rmd`), with code in each chunk and your answers to all questions.
- The figures you saved (`ex2_a.png` through `ex2_d.png`, and your `ex3_` figures).

## Grading Rubric (Total = 10 points)

### **Exercise 1 — Wrangle the Data (2 points)**

- Load and inspect DNM data (**0.5 pt**)
- Create per-proband maternal and paternal DNM counts (**0.5 pt**)
- Load and inspect parental age data (**0.5 pt**)
- Join counts with ages into a merged table (**0.5 pt**)

### **Exercise 2 — Fit and Interpret Linear Models (6 points)**

- Step 2.1: Scatter plots for DNMs vs. parental age and paternal vs. maternal age (**1 pt**)
- Step 2.2: Fit and interpret maternal OLS model (**1 pt**)
- Step 2.3: Fit and interpret paternal OLS model (**0.5 pt**)
- Step 2.4: Predict paternal DNMs for age 50.5 (**0.5 pt**)
- Step 2.5: Multiple regression with both parental ages, and nested model comparison with `anova()` (**1.5 pts**)
- Step 2.6: Plot distributions of maternal vs. paternal DNMs (**0.5 pt**)
- Step 2.7: Paired t-test (`t.test` and `lm(diff ~ 1)`) and interpret results (**1 pt**)

### **Exercise 3 — Explore a New Dataset (2 points)**

- Step 3.1: Choose and document TidyTuesday dataset (**0.5 pt**)
- Step 3.2: Produce exploratory figure(s) (**0.5 pt**)
- Step 3.3: Pose and test a multiple linear regression hypothesis and interpret results (**1 pt**)

**Total Points: 10**

------------------------------------------------------------------------
