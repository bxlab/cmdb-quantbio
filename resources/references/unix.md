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
