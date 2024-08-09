# QBB2024 - Day 1 - Lunch Exercises

## Overview

Biological Learning Objectives
- Explore the [Adult GTEx project](https://gtexportal.org) sample metadata
- Describe the results of your [data analysis](https://r4ds.had.co.nz/transform.html)

Computational Learning Objectives
- Explore text files using R and [dplyr](https://dplyr.tidyverse.org)
- Document your work on [GitHub](https://swcarpentry.github.io/git-novice/07-github.html)

## Instructions

Document your answers in a single file stored in `~/qbb2024-answers/day1-lunch`.

Please `git push` after each exercise and **do not wait** until the end of the session e.g.

```
git add explore-samples.Rmd
git commit -m "Add answer for exercise 1"
git push
# Confirm at github.com
``` 

## Exercises

0. Browse data dictionary

    - Navigate using a web browser to https://gtexportal.org/home/downloads/adult-gtex/overview
    - Select data type: `Metadata`
    - Download and take a quick look at `GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx`

1. Prepare your working environment

    - Reset R by opening the Session menu and selecting "Terminate R..."
    - Create either a new R Script (e.g. `explore-samples.R`) or a new R Notebook (e.g. `explore-samples.Rmd`) using the provided template
    - Load the `tidyverse` package

2. Wrangle the sample metadata

    - Load the `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` file and assign to variable `df`
    - Create a `SUBJECT` column using the following code

        ```
        df <- df %>%
          mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 )
        ```

    - Confirm that the first column in `df` is `SUBJECT`

3. Which two `SUBJECT`s have the most samples?  The least?

    - Use `group_by()`, `summarize()`, and `arrange()`
    - Describe your results as either a comment (`#` in R Script) or a bullet point (`-` in R Notebook)

4. Which two `SMTSD`s (tissue types) have the most samples?  The least?  Why?

    - Solve and document as above

5. For subject `GTEX-NPJ8`

    - Filter for samples from this subject and save as a new object (e.g. `df_npj8`)
    - Which tissue has the most samples?
    - For that tissue, what is different between the samples?  Scroll to the 15th through 20th columns ...

6. Explore `SMATSSCR` (autolysis score)

    - Filter out `NA` values in this column to avoid `mean()` returning `NA`
        ```
        df %>%
            filter( !is.na(SMATSSCR) )
        ```
    - How many samples have a mean score of 0?
    - What other observations can you make about the distribution of mean scores?
    - What are possible ways to present this information in a report?

## Just for fun

A. Identify another handful of columns using the data dictionary, explore the data, and describe your results

