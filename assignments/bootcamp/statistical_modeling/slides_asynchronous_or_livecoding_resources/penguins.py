#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy import stats

# read in penguins file
df = np.genfromtxt("penguins.csv", delimiter = ",",
                   dtype = None, encoding = None, names = True)    
                   
# identify rows where the species is Adelie
rows = np.where(df["species"] == "Adelie")
# make dataframe of only Adelie rows
df_adelie = df[rows]

# get lists of flipper length for male and female penguins
df_adelie_m = df_adelie[np.where(df_adelie["sex"] == "male")]
df_adelie_f = df_adelie[np.where(df_adelie["sex"] == "female")]

# plot distribution of flipper length for males and females
fig, ax = plt.subplots()
ax.hist(df_adelie_m["flipper_length_mm"], alpha = 0.5, label = "male")
ax.hist(df_adelie_f["flipper_length_mm"], alpha = 0.5, label = "female")
ax.set_xlabel("Flipper Length (mm)")
ax.set_ylabel("Number of penguins")
ax.legend()
plt.show()


# Is the distribution of flipper lengths different for male/female penguins?
# test using t-test
print(stats.ttest_ind(df_adelie_m["flipper_length_mm"],
                      df_adelie_f["flipper_length_mm"])))
# test using regression model, controlling for other variables
# (here we use the full penguins dataset so we can see if species has an effect)
model = smf.ols(formula = "flipper_length_mm ~ 1 + species + sex + year + island",
                data = df).fit()
results = model.fit()
print(results.summary())


# The previous model tested whether flipper length for each species was different
# than flipper length for Adelie penguins. How can we ask if there are flipper
# length differences between species, more generally?
full_model = smf.ols(formula = "flipper_length_mm ~ 1 + species + sex",
                     data = df).fit()
reduced_model = smf.ols(formula = "flipper_length_mm ~ 1 + sex",
                        data = df).fit()
print(sm.stats.anova_lm(reduced_model, full_model, typ = 1))


# What is predicted flipper length of a female penguin, island of Dream, year 2008?
# look at full model to get parameters for effects of island, etc.
print(full_model.summary())
# calculate by hand:
print(-5275.0999 + 0 + (2.7197 * 2008) + 1.4111)

# calcute using `predict` function
new_data = df[0] # select first row of dataframe
new_data.fill(0) # replace values in row with zeros
# add in prediction values of interest
new_data['island'] = 'Dream'
new_data['sex'] = 'female'
new_data['year'] = 2008
# predict
print(full_model.predict(new_data))