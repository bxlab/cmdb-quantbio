# Software Carpentry: Version Control with Git

- https://swcarpentry.github.io/git-novice

|  Time  | Activity |
|-------:|:---------|
| 13:00 | Automated Version Control |
|       | Setting Up Git |
|       | Creating a Repository |
|       | Tracking Changes |
|       | Exploring History |
|       | Ignoring Things |
| 14:10 | Break |
| 14:25 | Remotes in GitHub |
|       | Conflicts |
|       | Additional Content |
| 15:00 | [day2-homework](https://github.com/bxlab/cmdb-quantbio/blob/main/assignments/bootcamp/version-control-with-git/index.md) |

## Automated Version Control

- Images
    - ![](https://swcarpentry.github.io/git-novice/fig/phd101212s.png "notFinal.doc by Jorge Cham")
    - ![](https://swcarpentry.github.io/git-novice/fig/play-changes.svg "track changes")
    - ![](https://swcarpentry.github.io/git-novice/fig/versions.svg "branch")
    - ![](https://swcarpentry.github.io/git-novice/fig/merge.svg "merge")
- Exercises
    - [Paper Writing](https://swcarpentry.github.io/git-novice/instructor/01-basics.html#paper-writing)

## Setting Up Git

- Summary -- git config --global, vim, git help

## Creating a Repository

- Summary -- git init, .git, git checkout, git status
- Images
    - ![](https://swcarpentry.github.io/git-novice/fig/motivatingexample.png "Werewolf vs dracula")
- Exercises
    - [Places to create git repositories](https://swcarpentry.github.io/git-novice/instructor/03-create.html#places-to-create-git-repositories)
    - [Correcting git init mistakes](https://swcarpentry.github.io/git-novice/instructor/03-create.html#correcting-git-init-mistakes)

## Tracking Changes

- Summary -- nano mars.txt, git status, git add, git commit -m, git log, git diff
- Images
    - ![](https://swcarpentry.github.io/git-novice/fig/git-staging-area.svg "staging area")
    - ![](https://swcarpentry.github.io/git-novice/fig/git-committing.svg "staging area, two files")
- Exercises
    - [Choosing a commit message](https://swcarpentry.github.io/git-novice/instructor/04-changes.html#choosing-a-commit-message)
    - [Commiting changes to git](https://swcarpentry.github.io/git-novice/instructor/04-changes.html#committing-changes-to-git)
    - [Committing multiple files](https://swcarpentry.github.io/git-novice/instructor/04-changes.html#committing-multiple-files)
    - [bio repository](https://swcarpentry.github.io/git-novice/instructor/04-changes.html#bio-repository)

## Exploring History

- Summary -- git diff HEAD~1 mars.txt, git show, git diff \<commit\> mars.txt, git checkout
- Images
    - ![](https://swcarpentry.github.io/git-novice/fig/git-checkout.svg "git checkout")
    - ![](https://swcarpentry.github.io/git-novice/fig/git_staging.svg "how git works a cartoon")
- Exercises
    - [Recovering older versions of a file](https://swcarpentry.github.io/git-novice/instructor/05-history.html#recovering-older-versions-of-a-file)
    - [Reverting a commit](https://swcarpentry.github.io/git-novice/instructor/05-history.html#reverting-a-commit)
    - [Understanding workflow and history](https://swcarpentry.github.io/git-novice/instructor/05-history.html#understanding-workflow-and-history)
    - [Getting rid of staged changes](https://swcarpentry.github.io/git-novice/instructor/05-history.html#getting-rid-of-staged-changes)

## Ignoring Things

- Summary -- nano .gitignore
- Exercises
    - [Ignoring nested files](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#ignoring-nested-files)
    - [Including specific files](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#including-specific-files)
    - [Ignoring nested files: variation](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#ignoring-nested-files-variation)
    - [Ignoring all data files in a directory](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#ignoring-all-data-files-in-a-directory)
    - [Ignoring all data files in the repository](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#ignoring-all-data-files-in-the-repository)
    - [The order of rules](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#the-order-of-rules)
    - [Log files](https://swcarpentry.github.io/git-novice/instructor/06-ignore.html#log-files)

## Remotes in GitHub

- Summary -- https://github.com, git remote add origin, git push origin main, git pull origin main
- Images
    - ![](https://swcarpentry.github.io/git-novice/fig/git-freshly-made-github-repo.svg "~/vlad/planets/.git and github.com/vlad/planets.git")
    - ![](https://swcarpentry.github.io/git-novice/fig/github-repo-after-first-push.svg "git push origin main")
- Exercises
    - [Github gui](https://swcarpentry.github.io/git-novice/instructor/07-github.html#github-gui)
    - [Github timestamp](https://swcarpentry.github.io/git-novice/instructor/07-github.html#github-timestamp)
    - [Push vs commit](https://swcarpentry.github.io/git-novice/instructor/07-github.html#push-vs.-commit)
    - [Github license and readme files](https://swcarpentry.github.io/git-novice/instructor/07-github.html#github-license-and-readme-files)

## Conflicts

- Summary -- new line via https://github.com, new line locally and git push, git pull, nano mars.txt
- Images
    - ![](https://swcarpentry.github.io/git-novice/fig/conflict.svg "conflict")
- Exercises
    - [Conflicts on non-textual files](https://swcarpentry.github.io/git-novice/instructor/09-conflict.html#conflicts-on-non-textual-files)
    - [A typical work session](https://swcarpentry.github.io/git-novice/instructor/09-conflict.html#a-typical-work-session)
