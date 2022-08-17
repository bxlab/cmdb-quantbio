---
title: Git Recipes
layout: default
---


# Git Recipes

### Standard Usage

```
$ git status
$ git add file
$ git commit -m "Short comment on changes made"
$ git push
```

### Checkout older version of file

```
$ git log
$ git checkout <hash> path/to/file
$ git commit -m "Check out older version of file"
$ git push
```

### Un-commit large file (most recent commit)

```
$ git rm --cached path/to/largefile
$ git commit --amend -CHEAD
$ git push
```

### Un-commit large file (from a few commits ago)

```
$ mkdir backup
$ mv path/to/largefile backup
$ git filter-branch --tree-filter 'rm -f path/to/largefile' HEAD
$ git push

$ mv backup/largefile new/path/for/largefile
$ rm -r backup
```

### Rename/move files in Git repository

```
git mv path/to/file new/path/to/file
git commit -m "Move/rename file comment"
git push
```
