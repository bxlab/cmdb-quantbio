# Day 5 Lunch: Linear Regression

## Human de novo mutations

Data are taken from Halldorsson, B. V., Palsson, G., Stefansson, O. A., Jonsson, H., Hardarson, M. T., Eggertsson, H. P., ... & Gudjonsson, S. A. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science, 363(6425). [link](https://science.sciencemag.org/content/363/6425/eaau1043.abstract)

### Exercise 1: Wrangle the data with Unix

1. Read the abstract from the above paper to understand the context of the datasets we will be using. The relevant data are stored in two files:
* information about the number and parental origin of each de novo mutation detected in an offspring individual (i.e. "proband"), stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_dnm.csv`
* ages of the parents of each proband, stored in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/statistical_modeling/extra_data/aau1043_parental_age.csv`

2. Starting with the data in `aau1043_dnm.csv`, use Unix commands to count the number of de novo mutations per proband. The `Phase_combined` column records the inferred parent of origin of each de novo mutation. Break the counts of de novo mutations down into maternally inherited versus paternally inherited de novo mutations (ignore mutations of unknown parental origin). Store these counts in a **new file** with 3 fields containing: proband ID, number of paternally inherited de novo mutations for that proband, and number of maternally inherited de novo mutations for that proband. **Hint:** Parts of this problem will likely require `sort`, `uniq`, `join`, and `cut`.

3. Use the Unix `join` command to combine the above dataset that you created with the dataset containing ages of the mother and father of each proband (`aau1043_parental_age.csv`) at the probands time of birth.

### Exercise 2: Fit and interpret linear regression models with Python

4. Use numpy `genfromtxt` to load the "joined" data from step 3 into a numpy array. Use the `Names = ` option to give your fields informative names.

5. Use matplotlib to plot:
* the count of maternal de novo mutations vs. maternal age (upload as `ex2_a.png`)
* the count of paternal de novo mutations vs. paternal age (upload as `ex2_b.png`)

6. Use ordinary least squares `smf.ols()` to test for an association between maternal age and maternally inherited de novo mutations.
* Is this relationship significant?
* What is the size of this relationship?

7. Use ordinary least squares `smf.ols()` to test for an association between paternal age and paternally inherited de novo mutations.
- Is this relationship significant?
- What is the size of this relationship?

8. Plot a histogram of the number of maternal de novo mutations and paternal de novo mutations per proband on a single plot with semi-transparency.

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

