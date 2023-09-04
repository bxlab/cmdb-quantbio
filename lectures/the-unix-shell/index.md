# Software Carpentry: The Unix Shell

- https://swcarpentry.github.io/shell-novice

|  Time  | Activity |
|-------:|:---------|
|  8:45 | Arrive |
|  9:00 | [welcome-to-qbb2023](https://docs.google.com/presentation/d/1TJpwKrwHDkiC_0HTydT_3UtmkJ4zhokD_YRrL-mxsEU) |
|  9:30 | Summary and Setup |
|  9:45 | Introducing the Shell |
|       | Navigating Files and Directories |
| 10:30 | Break |
| 10:45 | Working With Files and Directories |
|       | Additional Content |
| 12:00 | Lunch |

## Summary and Setup

1. Search for Safari.app [using Spotlight](https://support.apple.com/guide/mac-help/search-with-spotlight-mchlp1008/mac)
1. Navigate to Software Carpentry's [The Unix Shell](https://swcarpentry.github.io/shell-novice) lesson
1. Download shell-lesson-data.zip
1. Use Finder.app to move the file to your Desktop
1. Launch Terminal.app

## Introducing the Shell

- Summary -- Prompt ($), ls, command not found

## Navigating Files and Directories

- Summary -- pwd, file system, ls -F, \<UP\>\<DOWN\>, man, cd, hidden files, paths, syntax, \<TAB\>
- Exercises
    - [Exploring more ls options](https://swcarpentry.github.io/shell-novice/02-filedir.html#exploring-more-ls-options)
    - [Listing in reverse chronological order](https://swcarpentry.github.io/shell-novice/02-filedir.html#listing-in-reverse-chronological-order)
    - [Absolute vs relative paths](https://swcarpentry.github.io/shell-novice/02-filedir.html#absolute-vs-relative-paths)
    - [Relative path resolution](https://swcarpentry.github.io/shell-novice/02-filedir.html#relative-path-resolution)
    - [ls reading comprehension](https://swcarpentry.github.io/shell-novice/02-filedir.html#ls-reading-comprehension)

## Working With Files and Directories

- Summary -- mkdir, ls -R, filenames, nano, mv, cp, rm, wildcards
- Exercises
    - [Creating files a different way](https://swcarpentry.github.io/shell-novice/03-create.html#creating-files-a-different-way)
    - [Moving files to a new folder](https://swcarpentry.github.io/shell-novice/03-create.html#moving-files-to-a-new-folder)
    - [Renaming files](https://swcarpentry.github.io/shell-novice/03-create.html#renaming-files)
    - [Moving and copying](https://swcarpentry.github.io/shell-novice/03-create.html#moving-and-copying)
    - [Using rm safely](https://swcarpentry.github.io/shell-novice/03-create.html#using-rm-safely)
    - [Copy with multiple filenames](https://swcarpentry.github.io/shell-novice/03-create.html#copy-with-multiple-filenames)
    - [List filenames matching a pattern](https://swcarpentry.github.io/shell-novice/03-create.html#list-filenames-matching-a-pattern)
    - [More on wildcards](https://swcarpentry.github.io/shell-novice/03-create.html#more-on-wildcards)
    - [Organizing directories and files](https://swcarpentry.github.io/shell-novice/03-create.html#organizing-directories-and-files)
    - [Reproduce a folder structure](https://swcarpentry.github.io/shell-novice/03-create.html#reproduce-a-folder-structure)

## Pipes and Filters

- Summary -- wc,  >, cat, sort, head, |
- Exercises
    - [What does sort -n do?](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#what-does-sort--n-do)
    - [What does >> mean?](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#what-does-mean)
    - [Appending data](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#appending-data)
    - [Piping commands together](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#piping-commands-together)
    - [Pipe reading comprehension](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#pipe-reading-comprehension)
    - [Pipe construction](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#pipe-construction)
    - [Which pipe?](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#which-pipe)
    - [Removing unneeded files](https://swcarpentry.github.io/shell-novice/04-pipefilter.html#removing-unneeded-files)

## Loops

- Summary -- for filename in *.dat; do echo $filename; done
- Exercises
    - [Write your own loop](https://swcarpentry.github.io/shell-novice/05-loop.html#write-your-own-loop)
    - [Variables in loops](https://swcarpentry.github.io/shell-novice/05-loop.html#variables-in-loops)
    - [Limiting sets of files](https://swcarpentry.github.io/shell-novice/05-loop.html#limiting-sets-of-files)
    - [Saving to a file in a loop - Part one](https://swcarpentry.github.io/shell-novice/05-loop.html#saving-to-a-file-in-a-loop---part-one)
    - [Saving to a file in a loop - Part two](https://swcarpentry.github.io/shell-novice/05-loop.html#saving-to-a-file-in-a-loop---part-two)
    - [Doing a dry run](https://swcarpentry.github.io/shell-novice/05-loop.html#doing-a-dry-run)
    - [Nested loops](https://swcarpentry.github.io/shell-novice/05-loop.html#nested-loops)

## Shell Scripts

- Summary -- nano middle.sh, bash middle.sh, $1, #
- Exercises
    - [List unique species](https://swcarpentry.github.io/shell-novice/06-script.html#list-unique-species)
    - [Why record commands in the history before running them?](https://swcarpentry.github.io/shell-novice/06-script.html#why-record-commands-in-the-history-before-running-them)
    - [Variables in shell scripts](https://swcarpentry.github.io/shell-novice/06-script.html#variables-in-shell-scripts)
    - [Find the longest file with a given extension](https://swcarpentry.github.io/shell-novice/06-script.html#find-the-longest-file-with-a-given-extension)
    - [Script reading comprehension](https://swcarpentry.github.io/shell-novice/06-script.html#script-reading-comprehension)
    - [Debugging scripts](https://swcarpentry.github.io/shell-novice/06-script.html#debugging-scripts)

## Finding Things

- Summary -- grep, find
- Exercises
    - [Using grep](https://swcarpentry.github.io/shell-novice/07-find.html#using-grep)
    - [Tracking a species](https://swcarpentry.github.io/shell-novice/07-find.html#tracking-a-species)
    - [Little women](https://swcarpentry.github.io/shell-novice/07-find.html#little-women)
    - [Matching and substracting](https://swcarpentry.github.io/shell-novice/07-find.html#matching-and-subtracting)
    - [find pipeline reading comprehension](https://swcarpentry.github.io/shell-novice/07-find.html#find-pipeline-reading-comprehension)

