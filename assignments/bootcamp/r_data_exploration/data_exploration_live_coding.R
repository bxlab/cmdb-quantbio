# highlight a line or block and press **Cmd/Ctrl + Return** to run it

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
1:5                 # create a vector of the numbers 1 through 5
2^(1:5)
x <- 1:5
2^x

# vectors are a basic data structure; all elements share the same type
# create with c()
c(10, 4, 2, 45.5)
c(10, 4, 2, TRUE)     # coerces to numeric (TRUE -> 1)
c(10, 4, 2, "hello")  # coerces to character

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
  fruits = c("raspberry", "banana", "mango"),
  is_berry = c(TRUE, FALSE, FALSE),
  inventory = c(40, 5, 3)
)

fruits_df <- data.frame(
  fruits = c("raspberry", "banana", "mango"),
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

## 10–20 min: load tabular data into a tibble with read_delim()

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

# Note: if a column that should be numeric loads as character (e.g., due to an NA), use `read_delim(..., na = c("", "NA"))`

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