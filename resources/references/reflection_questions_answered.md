# Answers to your FAQs

Each day in Bootcamp, we asked you to fill out a reflection and write any questions you had at the end of the day. We tried to address each of these throughout the course of the week - and hopefully many of these will sound familiar to you. We also thought we'd create a permanent resource for you to refer back to throughout Quant Bio. Please remember that these TAs and instructors will remain resources for you for, well, forever - so slack or email them any time. The [course website](https://andrew-bortvin.github.io/cmdb-bootcampNotes/index.html) also has extensive resources and examples for the materials we covered in class. 


## Python Syntax

**Indexing...**

Indexing is the way you'd subset your data. It lets you grab some element (e.g., a list, a string, etc.) by its position. For example, if you have `my_list = ['a', 'b', 'c']`, and you wanted to print b, you'd do `my_list[1]`. You can also index to get the last item, e.g., `my_list[-1]` would show 'c'.

**How can I use for loops instead of numpy to open and read files?**

In class so far we've implemented this as a 3-step process: 1) `open(filename)`. 2) `readlines()`. 3) Use a for loop to walk through each line of the file and do what you want (e.g., strip the newline character, split on whitespace or some delimiter, etc.). `.readlines()` is used to put the file into a massive list, but you can also loop over the file and do other operations, e.g., 
```
fs = open(filename)
for line in fs:
   print(line.rstrip())
fs.close() 
```
More on this in the [textbook](https://andrew-bortvin.github.io/cmdb-bootcampNotes/reading-in-data.html#parsing-a-file).

**What are some common applications of while loops?**

Tons of examples and you'll get a better insight the more you do it. One is if you want to get the first 1000 bases in a genome sequence. You'd walk through your sequence and record each nucleotide, and you'd keep count of how many nucleotides you'd seen so far: `num_nucs = 0`. `nucleotides = []`. `while num_nucs < 1000: for i in range(len(genome_seq)): num_nucs.append(genome_seq[i] \n num_nucs += 1`. In this case you'll append each nucleotide to your list, as long as you haven't yet done it 1000 times. 

**What is the use of a tuple, other than being able to serve as a key for a dictionary?**

A tuple can be useful for storing multiple items in a single variable. For example, if you are plotting the graphical coordinates of multiple objects and you want both the x and y position for each, you could store it as a list of tuples e.g., `my_objects = [(2, 3), (4, 1)]` would tell you that object 1 is at x = 2 and y = 3; and object 2 is at x = 4, y = 1. 

**How could I call out and manipulate a set of values from a dictionary like I would with a list?**

Since a list is ordered, you can index by position, e.g., `my_list[2]` will show you the third (because Python indexes from 0) object in your list. A dictionary is not ordered, so you can't look at objects quite the same way. Instead, you can look up a given value using its key. E.g., if you had a dictionary with amino acid/codon info, you could do `codon_table['ATG']` and it would print Methionine. To change the AA from Methionine to Start, since it's a start codon, you would do `codon_table['ATG'] = "Start"`. And if you somehow discovered a new amino acid named acidone and it's UDA (new nucleotide too!) you would add it to your dictionary via `codon_table['UDA'] = "acidone"`. 

**In what contexts would I use a dictionary in the real world instead of a list?**

A dictionary is useful for paired information: you have a key and a corresponding value: Amino acid name and its codon sequence. First name and last name. Individual and their height. The alternative is to have this structure be unpaired, which is worse; a dictionary keeps the info together and lets you use both pieces of information. 

**What does += mean, and how is it used?**

It is a shorthand: `x += 1` is the same as writing `x = x + 1`. You are assigning to the variable x a new value, which is one plus the existing value for x. It's commonly used in a for loop as an iterator or counter (e.g., to keep track of how many times you've done something). 

**Are there any times that you should not use a dictionary?**

A dictionary is most useful for when you have paired information. You'll build an instinct for which data structure is the best to use the more you code. 


## Functions 

**What are the best ways to share functions among Python files? Say you generate a function and save the script in a specific directory. If you're working in a separate directory, how do you call that function? Is it still just import, or do you need to specify where the function script is stored more precisely?.**

In Python, you can import a function from another file into your current script (as long as the other function is in a file with a .py suffix). Maybe you made a useful function yesterday and you want to use it again today. You'll import the name of the file the function was in, as well as the name of the function itself. Then you can use it in your current script freely. See [tutorial](https://problemsolvingwithpython.com/07-Functions-and-Modules/07.05-Calling-Functions-from-Other-Files/#:~:text=To%20use%20the%20functions%20written,of%20the%20filename%20during%20import.).

**Can I define a function within another function that I am defining? If not, why?**

You [can](https://realpython.com/inner-functions-what-are-they-good-for/#:~:text=Inner%20functions%2C%20also%20known%20as,closure%20factories%20and%20decorator%20functions), but for this class it's probably cleaner to keep your functions separate. If you really want to define a new function inside another function, just do an extra few checks (e.g., inspect your outputs via `print()`, `len()` etc.) to make sure both worked as you'd intended. 


## Numpy 

**I am wondering what other functionality numpy has.**

So much functionality! They have built-in functions for different mathematical steps applied to numpy arrays, and ways of reshaping data. More info is in the [manual](https://numpy.org/doc/stable/user/).

**How we can play with numpy when we have super huge dataset?**

Numpy is really handy for doing math across huge datasets and reading these data into a structure (array) with the rows and columns arranged how you want. You can skim the [manual](https://numpy.org/doc/stable/user/) and check out tutorials online. In the fall Quant Bio class you'll be introduced to other uses as well. 


## General Quantitative and Computational Biology 

**Is there a simple way to make more aesthetically pleasing plots this way? Even Excel and Powerpoint graphs look nicer.**

You'll do a ton of this in the fall quant bio course - making aesthetic plots is a huge part of science and communicating your work. Much nicer than Excel, trust me. You'll get a lot more practice with matplotlib and other packages such as seaborns.

**I am wondering why we need to use Sublime and the terminal to execute Python script. What is the benefit of this method over using a program that can both be a text manager for Python and an executor?**

Sublime is a text editor - it's basically microsoft word without any of the annoying formatting things. It's just a place you type your script. You could use dozens of other text editors (some common ones are vim, nano, textmate); each has their own unique features and keyboard shortcuts. We recommended Sublime because if you use a .py suffix on your file, it'll automatically apply Python styling (colored text, comments as # , etc.). You have to execute your scripts from the command line. Terminal is the app on your macbook that is your command line (there are others of these too but for this course you should stick with terminal). You could run Python code directly from the command line (remember if you type `python` it'll open a python window), but then you'll be unable to easily modify and save your code without rerunning it. Scripts are useful because they contain the code you've done and let you quickly modify it and then immediately re-run (in a single command, right? ./myscript.py) from the command line. If you're newer to coding we recommend keeping them separate (i.e., a text editor + Terminal). If you've done this before, and have a preferred program that does both (text editor + command line) feel free to use it. If you're new but want to check out a combined text editor/command line program, Slack a TA. 

**How can I best optimise the things I learned today?**

This is a great question and something we're all still working on! As you get more comfortable coding, you'll learn keyboard shortcuts, best practices, and you'll develop an intuition for what code structures to use when. 


## GitHub and git 

**Are GitHub personal access tokens (PATs) device-specific?**

They are not device-specific. If you log into your git account on another laptop and connect a git repository, you can login using the same token you used on your bootcamp laptop. 

**How do you pull a file from git, once it’s on there? Specifically, if my colleague pushed something, how do I pull it to my computer? What does it look like if I make changes and push it back?**

For the purposes of this course, you probably won't interact live with files that are on GitHub or files that your collaborators are modifying. In real life (maybe in your research), you would modify a file on something called a `git branch`, and then you'd push that file to GitHub. GitHub would note that it's on a certain branch (i.e., not the main branch). If a collaborator wanted to edit it, they could either git pull your new branch, modify the file, and push it back; you'd then git pull you branch again (to collect those updates on to your local) and then deal with the file. Or, a collaborator could make their own branch, modify any files there, and then push their branch back into main; for you to get the updates from main to your local, you'll do a git pull. This is a lot of info and while it's intuitive, it's also new. Feel free to slack any of the TAs any time (even after this course) for more of a git tutorial. We also have a section on Git in the [course website](https://andrew-bortvin.github.io/cmdb-bootcampNotes/git.html). 




