# =============================================================================
# BLOCK 1 (LIVE CODING) -- Data exploration in R
#
# Highlight a line or block and press Cmd/Ctrl + Return to run it.
#
# A note on chatbots: you used HopGPT yesterday, and you should keep using it
# today. Two rules for the rest of the course:
#   1. You may not submit code you cannot explain.
#   2. Anything it tells you about *your data* gets checked against your data.
# When questions come up as we go, we'll take them to HopGPT together.
# =============================================================================


# R has many built-in operators (similar to Python)
1 + 1
10 ^ 3      # exponents
5 %/% 2     # integer (floor) division

# comparison
1 == 1
1 != 2      # not equal to
1 >= -9     # greater than or equal to

# assignment: R commonly uses <- (arrow), though = also works in many contexts
x <- 1/40
x
x <- x + 1  # notice the variable gets updated in the Environment tab

# R is vectorized: many operations work elementwise on vectors
# (in Python you would need a loop here)
1:5                 # create a vector of the numbers 1 through 5
2^(1:5)
x <- 1:5
2^x

# vectors are a basic data structure; all elements share the same type
# create with c()
c(10, 4, 2, 45.5)
c(10, 4, 2, TRUE)     # coerces to numeric (TRUE -> 1)
c(10, 4, 2, "hello")  # coerces to character

# indexing a vector
my_vector <- c(10, 4, 2, 45.5)
my_vector[2] # note that R is 1-based, where Python is 0-based
my_vector[length(my_vector)]
my_vector[1:2]
my_vector[c(TRUE, FALSE, TRUE, FALSE)]

# common atomic types in R: double, logical, integer, character, complex
typeof(2)
typeof(2.1)
typeof("hello, world")
# typeof(True)   # error; R uses TRUE/FALSE (all caps)
typeof(TRUE)
typeof("TRUE")
typeof(2L)       # L denotes integer
typeof(1 + 4i)   # complex

# lists can store heterogeneous types
my_list <- list(10, "hello", c(1, 2, 3), TRUE)
print(my_list)

# data frames: columns can have different types (each column itself is a vector)
data.frame(
  fruits = c("raspberry", "strawberry", "mango"),
  is_berry = c(TRUE, FALSE, FALSE),
  inventory = c(40, 5, 3)
)

fruits_df <- data.frame(
  fruits = c("raspberry", "strawberry", "mango"),
  is_berry = c(TRUE, FALSE, FALSE),
  inventory = c(40, 5, 3)
)

# index elements of a data.frame as df[row, col]
fruits_df[1, 2]
fruits_df[1, 1:2]
fruits_df[, 1]

# extract an entire column as a vector with df$column_name
fruits_df$fruits

# tidyverse is a popular set of packages for data manipulation and more
# it includes an upgraded data.frame called a tibble (e.g., better printing)
library(tidyverse)

fruits_tibble <- as_tibble(fruits_df)
fruits_tibble[1, 2]
fruits_tibble[1, 1:2]
fruits_tibble[, 1]

## 10-20 min: load tabular data into a tibble with read_delim()

library(tidyverse)  # tidyverse is actually a collection of other packages (dplyr, ggplot2, readr, etc.)

# read_delim() will accept a file path or even a URL as its first argument
penguins <- read_delim(
  "https://gist.githubusercontent.com/slopp/ce3b90b9168f2f921784de84fa445651/raw/4ecf3041f0ed4913e7c230758733948bc561f434/penguins.csv",
  delim = ","
)

penguins            # tibbles print nicely
dim(penguins)       # rows, columns
head(penguins, 6)   # first 6 rows
colnames(penguins)  # column names
glimpse(penguins)   # compact summary

# Look at the column types every time you load a file. If a column that should
# be numeric loads as character (e.g., due to an NA), use
# `read_delim(..., na = c("", "NA"))`

## 10 min: select() and filter() from dplyr

# %>% is the pipe character
# R also now has a built-in pipe operator, |>, but we will use %>% as it is the original from the tidyverse

# select a few columns
penguins %>%
  select(species, island, bill_length_mm)

# exclude columns by name
penguins %>%
  select(-year, -sex)

# filter rows by equality
penguins %>%
  filter(species == "Adelie")

# filter with multiple conditions (and)
penguins %>%
  filter(species == "Adelie" & island == "Biscoe")

# filter with multiple conditions (or)
penguins %>%
  filter(species == "Adelie" | island == "Biscoe")

# filter using comparisons
penguins %>%
  filter(bill_length_mm > 45)

# filter using %in%
penguins %>%
  filter(species %in% c("Adelie", "Gentoo"))

## a short but important detour: filtering and missing data

# some penguins have a missing (NA) value for sex
penguins %>%
  count(sex)

# now watch what happens when we ask for the penguins that are NOT male
nrow(penguins)
penguins %>% filter(sex == "male") %>% nrow()
penguins %>% filter(sex != "male") %>% nrow()

# these two do not add up to the total. `NA != "male"` is not TRUE and it is
# not FALSE -- it is NA -- and filter() keeps only the rows that are TRUE.
# The penguins with missing sex silently disappeared. No error, no warning.

# if you want them, say so:
penguins %>%
  filter(is.na(sex) | sex != "male") %>%
  nrow()

# Habit to build: check your row count before and after any filter or join.
# Quietly losing rows is the most common way to get a wrong answer that
# looks completely reasonable.

## 10 min: mutate() and arrange()

# create a new column based on old columns with mutate
penguins %>%
  mutate(bill_ratio = bill_length_mm / bill_depth_mm)

# multiple new columns at once
penguins %>%
  mutate(
    bill_ratio = bill_length_mm / bill_depth_mm,
    flipper_cm = flipper_length_mm / 10
  )

# arrange (sort) ascending
penguins %>%
  arrange(bill_length_mm)

# Arrange descending
penguins %>%
  arrange(desc(bill_length_mm))

penguins %>%
  arrange(-bill_length_mm)

## 10 min: group_by() and summarise()

# Mean bill length by species
penguins %>%
  group_by(species) %>%
  summarize(mean_bill = mean(bill_length_mm, na.rm = TRUE))

# Multiple summaries by species and island
penguins %>%
  group_by(species, island) %>%
  summarise(
    mean_flipper = mean(flipper_length_mm, na.rm = TRUE),
    sd_flipper = sd(flipper_length_mm, na.rm = TRUE),
    n = n()
  )
# note that the tibble output by the above code is still grouped by species
# to remove all grouping from the output (e.g., if you want to use it in a subsequent step),
# you can use the .groups = "drop" argument


# example: proportion of samples for which body mass measurements are missing by species
# (this works because the mean of a TRUE/FALSE vector is a proportion)
penguins %>%
  group_by(species) %>%
  summarise(
    prop_missing_mass = mean(is.na(body_mass_g))
  )

## 10 min: left_join(), pivot_longer(), pivot_wider()

# create some species metadata
nests <- tibble(
  species = c("Adelie", "Gentoo", "Chinstrap"),
  nest_type = c("pebbles", "stones", "grass")
)

# Left join by the key column 'species'
penguins_joined <- penguins %>%
  left_join(nests, by = "species")

# again: check that the join did what you expected
nrow(penguins)
nrow(penguins_joined)

# pivot_longer(): go from wide to long format
penguins_long <- penguins %>%
  pivot_longer(
    cols = c(flipper_length_mm, body_mass_g),
    names_to = "trait",
    values_to = "value"
  )

penguins_long %>% head()

# pivot_wider(): go from long to wide format
penguins_wide <- penguins_long %>%
  pivot_wider(
    names_from = trait,
    values_from = value
  )

penguins_wide %>% head()


## 20 min: writing your own functions

# You have written functions in Python. The idea is identical here; only the
# syntax differs. Write one whenever you find yourself copying a block of code
# and editing a name or a number in the copy -- that edit is exactly where
# mistakes hide, because the result is a plausible number rather than an error.

# In Python:                        In R:
#
#   def celsius_to_f(temp):           celsius_to_f <- function(temp) {
#       return temp * 9/5 + 32            return(temp * 9/5 + 32)
#                                       }
#
# A function is a value, so we assign it with <- like anything else.

celsius_to_f <- function(temp) {
  return(temp * 9/5 + 32)
}

celsius_to_f(0)
celsius_to_f(100)

# R will also return the last expression automatically, so you will see
# functions written without return(). We will always write return() explicitly:
# it is clearer, and it matches what you already know from Python.

# The vectorization from the top of the script comes along for free:
celsius_to_f(c(-40, 0, 37, 100))


## arguments and default values

# a function that counts how many values are above some cutoff
count_above <- function(values, cutoff) {
  n <- sum(values > cutoff, na.rm = TRUE)
  return(n)
}

count_above(penguins$body_mass_g, 4500)
count_above(penguins$flipper_length_mm, 200)

# give an argument a default by assigning to it in the function definition
count_above <- function(values, cutoff = 4000) {
  n <- sum(values > cutoff, na.rm = TRUE)
  return(n)
}

count_above(penguins$body_mass_g)          # uses the default
count_above(penguins$body_mass_g, 4500)    # overrides it

# arguments match by position first, then by name -- name them when it helps
count_above(values = penguins$body_mass_g, cutoff = 4500)


## always test a new function on an input whose answer you already know

count_above(c(1, 2, 3, 4, 5), cutoff = 3)   # should be 2

# and check what it does at the boundary and with missing values
count_above(c(1, 2, 3), cutoff = 3)         # is > what you wanted, or >= ?
count_above(c(1, 2, NA, 4), cutoff = 1)

# This habit is worth more than it looks. It is also exactly how you check
# code that a chatbot wrote for you.


## functions can return more than a single number: return a tibble

trait_summary <- function(values) {
  result <- tibble(
    n = sum(!is.na(values)),
    mean = mean(values, na.rm = TRUE),
    sd = sd(values, na.rm = TRUE),
    median = median(values, na.rm = TRUE)
  )
  return(result)
}

trait_summary(penguins$body_mass_g)
trait_summary(penguins$flipper_length_mm)
trait_summary(penguins$bill_length_mm)

# One function, three summaries, no copied code to keep in sync.


## functions can call your other functions

body_mass_report <- function(values, cutoff = 4000) {
  result <- trait_summary(values)
  result$n_above_cutoff <- count_above(values, cutoff)
  return(result)
}

body_mass_report(penguins$body_mass_g)
body_mass_report(penguins$body_mass_g, cutoff = 5000)


# -----------------------------------------------------------------------------
# Now: the assignment, using median gene expression across human tissues
# from GTEx.
# -----------------------------------------------------------------------------
