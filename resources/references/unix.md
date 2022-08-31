---
title: Unix Cheat Sheet
layout: default
---

<!-- Make columns even -->
<style>
  td:first-child { width: 20% }
</style>

## Unix Cheat Sheet

#### FINDING HELP

|:---:|:--------------------------------------------|
| `man` | format and display the on-line manual pages |
|     | `man [section] name ...`                      |
|     | (`q`)uit (`/`) search                           |
{:.table.table-striped}

#### SPECIAL KEY (COMBINATIONS)

|:---------:|:---------------------------------|
| `<TAB>`   | Auto-complete                    |
|  `<UP>`   | Scroll through previous commands |
|`<CTRL>-C` | Cancel                           |
|`<CTRL>-A` | Jump to beginning of line        |
{:.table.table-striped}

#### NAVIGATION and ORGANIZATION

|:-----:|:----------------------------------------------------------------------|
| `ls`  | list directory contents                                               |
|       | `ls [-ABCFGHLOPRSTUW@abcdefghiklmnopqrstuwx1] [file ...]`             |
| `cd`  | Change the current directory to DIR                                   |
|       | `cd [-L\|-P] [dir]`                                                   |
| `cp`  | copy files                                                            |
|       | `cp [-R [-H \| -L \| -P]] [-fi \| -n] [-apvX] source_file target_file`|
| `mv`  | move files                                                            |
|       | `mv [-f \| -i \| -n] [-v] source target`                              |
| `rm`  | remove directory entries                                              |
|       | `rm [-dfiPRrvW] file ...`                                             |
| `mkdir`| make directories                                                      |
|       | `mkdir [-pv] [-m mode] directory_name ...`                            |
|  `.`  | Current directory                                                     |
| `..`  | Parent directory                                                      |
|  `~`  | Home directory                                                        |
|  `/`  | Root directory                                                        |
|  `*`  | Wildcard                                                              |
{:.table.table-striped}

#### EXPLORE TEXT FILES

|:----:|:---------------------------------------------------------------------------|
|`head`| display first lines of a file                                              |
|      | `head [-n count \| -c bytes] [file ...]`                                   |
|`less`| program used to view  contents of text files one screen at a time          |
|      | `less [-...N...S...]`                                                      |
|      | (`q`)uit  (`/`)search                                                      |
| `wc` | word, line, character, and byte count                                      |
|      | `wc [-clmw] [file ...]`                                                    |
|`grep`| file pattern searcher                                                      |
|      | `grep [-...c...vw...] [-A num] [-B num] [-f file] [pattern] [file ...]`    |
|`cut` | cut out selected portions of each line of a file                           |
|      | `cut -f list [-d delim] [-s] [file ...]`                                   |
|`sort`| sort lines of text files                                                   |
|      | `sort [OPTION]... [FILE]...`                                               |
|`uniq`| report or filter out repeated lines in a file                              |
|      | `uniq [-c \| -d \| -u] [-i] [-f num] [-s chars] [input_file [output_file]]`|
{:.table.table-striped}

#### INPUT OUTPUT REDIRECTION

|:----:|:-------------------------------------------------------|
| <code>&#124;</code>  | send output of one command to another command as input |
| `>`  | save output to a new file                              |
| `>>` | append output to an existing file                      |
| `<`  | use file as input                                      |
{:.table.table-striped}

#### bash scripts

#### awk

`awk` can be used to `cat` a file, or return all of its contents, but it can also be used to subset and only print parts of a file. Specific lines or specific columns based on some condition.

`awk` expressions use single quotation marks and curly braces.

**Print the whole file:**

`awk '{print}' input_file`

**Print specific columns:**

This example will print just the 1st, 3rd, and 5th columns of every line of the input file.

`awk '{print $1,$3,$5}' input_file`

**Checking a condition:**

If you want to print a line only if a column contains a specific value, then you would use an if conditional statement within the awk expression.

This example will only print a line from the input_file if the value in the first column is "chr21".

`awk '{if ($1 == "chr21") {print}}' input_file`

**Passing a variable:**

If you want to use a variable that is predefined within a script within your awk command, you have to set the identity of that variable with the `-v` flag. Variables set within a script can't be referenced within the single parentheses of the `awk` expression. In the example below, a line will be printed if the first column has "apple" in it.

```
FRUIT="apple"
awk -v fruit=$FRUIT '{if ($1 == fruit) {print}}' input_file
```

**Specifying the delimiter or field separator of output:**

This example sets the output field separator (OFS) to make the printed output tab-delimited ("\t")

```
awk 'BEGIN{OFS="\t"} {$1=$1; print}' input_file
```

This example sets the output field separator (OFS) to make the printed output (only the first 3 columns of the input files) separated by a comma (",")

```
awk 'BEGIN{OFS=","} {print $1,$2,$3}' input_file
```
