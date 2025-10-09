library(tidyverse)
library(palmerpenguins)
library(broom)

head(penguins)

glimpse(penguins)

table(penguins$species)
table(penguins$species, penguins$island)
table(penguins$species, penguins$island, penguins$sex)

#What questions might we want to ask?
#  Some examples...
#  Do male and female Adelie penguins have different bill lengths?
#  Do Adelie penguins on different islands have different flipper lengths?
#  Have bill lengths of Chinstrap penguins changed over the years?

adelie_penguins <- penguins %>%
  filter(species == "Adelie")

ggplot(data = adelie_penguins %>% filter(!is.na(sex)), aes(x = flipper_length_mm, fill = sex)) +
  geom_histogram(position = "identity", alpha = 0.5)

male_adelie_penguins <- adelie_penguins %>%
  filter(sex == "male")

female_adelie_penguins <- adelie_penguins %>%
  filter(sex == "female")

# Student's t-test
t.test(male_adelie_penguins$flipper_length_mm, female_adelie_penguins$flipper_length_mm, var.equal = TRUE)

model_results <- lm(data = adelie_penguins, 
   formula = flipper_length_mm ~ sex) %>%
  summary() %>%
  tidy()

model_results$p.value[2] # extract the p-value

# what if this is confounded by year?
# include year as a covariate in a multiple linear regression

lm(data = adelie_penguins, 
   formula = flipper_length_mm ~ sex + year) %>%
  summary() %>%
  tidy()

# what if island is also potentially important?

lm(data = adelie_penguins, 
   formula = flipper_length_mm ~ sex + year + island) %>%
  summary() %>%
  tidy()

# what is being treated as the "reference level" of island?
table(adelie_penguins$island)

# how can we test if not the difference between each island and Biscoe, 
# but whether there are any differences across islands at all?
# in other words, is island a significant predictor?

full_model <- lm(data = adelie_penguins, 
                 formula = flipper_length_mm ~ sex + year + island)

reduced_model <- lm(data = adelie_penguins, 
                    formula = flipper_length_mm ~ sex + year)

anova(full_model, reduced_model)

## Using your model to make predictions
# What is the expected flipper length of a female penguin on the island of Dream in year 2008?

summary(full_model)
-5275.0999 + (2.7197 * 2008) + 1.4111

new_obs <- tibble(sex = "female", island = "Dream", year = 2008)

predict(full_model, newdata = new_obs)


ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm, color = species)) +
  geom_point() +
  stat_smooth(method = "lm")




