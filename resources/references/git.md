---
title: Git Cheat Sheet
layout: default
---

<!-- Make columns even -->
<style>
  td:first-child { width: 20% }
</style>

## Git Cheat Sheet [by ZeroTurnaround](http://zeroturnaround.com/wp-content/uploads/2016/05/Git-Cheat-Sheet-by-RebelLabs.png)

##### FINDING HELP

|  `man git-add`   | join git command with a dash (`-`) to format and display the on-line manual pages |
| `man git-commit` |                                                                                   |
|       etc        |                                                                                   |
{:.table.table-striped}

##### CORE COMMANDS

| `git status` | Show the working tree status                     |
|              | `git status [<options>...] [--] [<pathspec>...]` |
|  `git add`   | Add file contents to the index                   |
|              | `git add [options] [--] [<pathspec>...]`         |
| `git commit` | Record changes to the repository                 |
|              | `git commit [-m <msg>]`                          |
|  `git push`  | Update remote refs along with associated objects |
|              | `git push [<repository> [<refspec>...]]`         |
{:.table.table-striped}

##### MORE COMMANDS

| `git log`  | Show commit logs                                                              |
|            | `git log [<options>] [<revision range>] [[--] <path>...]`                     |
| `git diff` | Show changes between commits, commit and working tree, etc                    |
|            | `git diff [options] [<commit>] [--] [<path>...]`                              |
|  `git rm`  | Remove files from the working tree and from the index                         |
|            | `git rm [-f] [-n] [-r] [--cached] [--ignore-unmatch] [--quiet] [--] <file>...`|
|  `git mv`  | Move or rename a file, a directory, or a symlink                              |
|            | `git mv [-v] [-f] [-n] [-k] <source> <destination>`                           |
{:.table.table-striped}
