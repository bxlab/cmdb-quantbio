# Assigment 1 - Python basics

### Overview
------------

Genomic annotations are often represented in the browser extensible data (BED) file format. This consists of genomic location data and may include labels, scores, strandedness, and other features for each entry. In this exercise, you will be looking at data from a basic BED file containing gene boundaries and expression levels to calculate some basic descriptive statistics.

Biological Learning Objectives

- Understand how interval data is represented in the BED file format
- Use descriptive statistics to characterize a set of gene annotations

Computational Learning Objectives

- Create, modify, and access data from lists
- Use built-in string functions to manipulate text
- Convert between different data types
- Step through list items using a `for` loop
- Condition execution of code based on an `if` statement

### Instructions
----------------

- Save a copy of this notebook as `~/qbXX-answers/python_basics.ipynb`.
- Fill in answers in the available code cells below.
- Remember to comment your code to help yourself and us know what each part is intended to do.
- You may find the functions `sorted()`, `sum()`, and `len()`, which all accept a list as their argument, useful for completing this exercise

```python
# Use this as your input data

text_data = "chrI\t8377406\t8390027\tNM_059873.7\t6\t-\nchrI\t8377598\t8392758\tNM_182066.7\t502\t-\nchrI\t8377600\t8392768\tNM_001129046.3\t7\t-\nchrMt\t1041473\t1049600\tNM_058410.5\t476\t+\nchrI\t3144409\t3147793\tNM_058707.5\t531\t-\nchrI\t4193240\t4203303\tNM_058924.6\t673\t+\nchrI\t6284972\t6294057\tNM_059412.6\t532\t+\nchrI\t6289432\t6294068\tNM_001374900.1\t615\t+\nchrI\t6290315\t6293988\tNM_001374901.1\t709\t+\nchrMt\t7339231\t7345684\tNM_001136303.4\t7\t+\nchrI\t9431248\t9441017\tNM_060105.7\t481\t+\nchrI\t9435464\t9441029\tNM_060106.7\t3\t+\nchrMt\t11526917\t11557854\tNM_060545.6\t498\t-\nchrI\t11526922\t11552160\tNM_001383731.1\t536\t-\nchrI\t11527027\t11557792\tNM_001306313.4\t6\t-\nchrI\t11527076\t11541127\tNM_001306312.3\t7\t-\nchrI\t11527076\t11546567\tNM_001306311.3\t694\t-\nchrI\t11527076\t11552106\tNM_001306310.3\t2\t-"

# This is what the data would look like in a file
# chrI    8377406 8390027 NM_059873.7     6       -
# chrI    8377598 8392758 NM_182066.7     502     -
# chrI    8377600 8392768 NM_001129046.3  7       -
# chrMt   1041473 1049600 NM_058410.5     476     +
# chrI    3144409 3147793 NM_058707.5     531     -
# chrI    4193240 4203303 NM_058924.6     673     +
# chrI    6284972 6294057 NM_059412.6     532     +
# chrI    6289432 6294068 NM_001374900.1  615     +
# chrI    6290315 6293988 NM_001374901.1  709     +
# chrMt   7339231 7345684 NM_001136303.4  7       +
# chrI    9431248 9441017 NM_060105.7     481     +
# chrI    9435464 9441029 NM_060106.7     3       +
# chrMt   11526917        11557854        NM_060545.6     498     -
# chrI    11526922        11552160        NM_001383731.1  536     -
# chrI    11527027        11557792        NM_001306313.4  6       -
# chrI    11527076        11541127        NM_001306312.3  7       -
# chrI    11527076        11546567        NM_001306311.3  694     -
# chrI    11527076        11552106        NM_001306310.3  2       -

```

### Excercises
--------------

1. Parse the data from the data string, converting values to appropriate data types as necessary
    
    - Split the string by the newline character (`\n`) to create a list containing one line per list-entry
    - Use a `for` loop to step through the data, line by line
    - Split each line by the tab characer (`\t`) to get values for ecah column
    - Convert position values to integers and expression values to floats

```python
# Split the data string into a list of individual lines

# Step through each line in the list

    # Split the line into a list of individual fields

    # Assign each field to a variable (chrom, start, stop, etc.), converting the data type if necessary

```

2. Create and populate lists with the parsed data using one to record expression values, one to record gene lengths, and one to record strandedness, skipping any genes from the mitochondrial chromosome

    - Use your above code as a starting point for this, adding lines inside the `for` loop for recording relevant values
    - It will be much easier to analyze the strandedness data if you convert it into ones and zeros (you can use an `if` statement for this)

```python
# Split the data string into a list of individual lines

# Create a list for recording expression data, one for recording gene lengths, and one for recording strandedness

# Step through each line in the list

    # Split the line into a list of individual fields

    # Assign each field to a variable (chrom, start, stop, etc.), converting the data type if necessary

    # Check if gene is from "chrMt" and if so, skip

    # Add the expression value, gene size, and strandedness to their corresponding lists

```

3. Calculate and report the 1st, 2nd, and 3rd quartiles of gene expression

    - Your expression values will need to be in order to determine quartiles
    - The first quartile is the 0.25 * (n + 1)th item, the second quartile is the 0.5 * (n + 1)th item, etc.
    - Because python starts counting at zero, the ith item is at index i - 1
    - Remember that the get a value from the list at specific index, you can only use intergers, not floats
    - Use a `print` statement to report your answer

```python
# Sort your expression values from lowest to highest

# Determine which position in the list corresponds to the first, second, and third quartiles

# Find the value from those positions in your expression list

# Report your results

```

4. Calculate and report the minimum and maximum gene size

```python
# Sort your gene sizes values from lowest to highest

# Find the largest and smallest gene sizes

# Report your results

```

5. Calculate and report the percent of genes that are on the positive strand

```python
# Calculate the percentage of positive strand genes

# Report your results

```
