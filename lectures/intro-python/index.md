# Software Carpentry: Python Novice

- [https://swcarpentry.github.io/python-novice-inflammation/instructor/01-intro.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/01-intro.html)

## Motivation (1:00-1:10; PowerPoint)

- State learning objectives
- What is a program?
- Why Python?
- [Setup](https://swcarpentry.github.io/python-novice-inflammation/#obtain-lesson-materials)

## Executing code (1:10-1:40; live coding)

1. Interactively via the interpreter (1:10-1:25)
- demonstrate use as a calculator
- demonstrate a print statement "Hello, World!"
- explain that print is a built-in function and show them the help page
- demonstrate statements that would produce an error
- demonstrate what happens if you press enter before completing a statment

Give students time to experiment on their own and emphasize the power of experimenting!

2. As scripts via the command line (1:25-1:40)
- demonstrate conversion of "Hello, World!" to a script written in Sublime
- what happens if I don't have a print statement?
- demonstrate a script that would produce an error
- have students copy and run more complex script
- note basic aspects of Python syntax (indentation, line breaks, comments)

Allow students to practice and experiment

## Basic data types (1:40-2:10; live coding)

1. Numbers (int, float, complex) (1:40-1:50)
- describe how to specify and convert among types
- demonstrate what happens when you perform operations on mixed types

2. Booleans (and comparison operators) (1:50-1:55)
- show that booleans are specified without quotes
- demonstrate how comparisons produce Booleans
- demonstrate use of various comparison operators
- emphasize difference between `=` for assignment and `==` for comparison

3. Strings (don't introduce this as a sequence/interable yet) (1:55-2:00)
- show that strings are specified with quotes
- show some basic operations that are possible with strings (e.g., `+`)
- show what happens when you try to add a string and a number
- show conversion of a string to a number
- show some string methods, such as `.toupper()` as an intro to the concept of methods

Allow students to practice specifying, converting, and checking types (2:00-2:10)

## Variables (2:10 - 2:25; live coding)

1. Creating variables
- use SWCarpentry example for creating variables (`weight_kg = 65`) https://swcarpentry.github.io/python-novice-inflammation/instructor/01-intro.html#variables
- demonstrate that you can check the type of a variable
- demonstrate that you can perform operations on variables
- show that simply performing the operation does not change the variable
- demonstrate that you can copy the value of one variable to a new variable, but these are subsequently unrelated

2. Changing variables
- use SWCarpentry example of how to create new variables via operations on existing variables (`weight_lb = 2.2 * weight_kg`)
- show that does not alter the original variable
- show how you can update the value of a variable

## Lists (2:25 - 2:55; live coding)

Work through SWCarpentry materials: https://swcarpentry.github.io/python-novice-inflammation/instructor/04-lists.html

1. Creating lists (2:25-2:30)
- from scratch
- with `my_string.split()`
- check the length of the list with `len()`

2. Indexing lists (2:30-2:40)
- basic integer indexing (zero-based)
- negative integers as indices
- slicing

3. Modifying lists (2:40-2:50)
- changing specific elements based on indexing
- list methods such as `.append()`, `.pop()`, and `.reverse()` - use this to reinforce the concept of methods
- combining lists

4. Nested lists (2:50-2:55)
- lists can even contain lists as elements
- how to index a nested list

## `for` loops (2:55 - 3:45; live coding)

1. Motivate the use of loops for performing a repeated action - copy and paste print statements (2:55 - 3:00)

2. Demonstrate `for` loop syntax (3:00-3:10)
- start with basic print statement, interating through a list
- note indentation...what happens if the print statement is not indented (i.e., outside the loop)?
- show that the list can be a variable
- show that the name you use for the index variable can be anything (and use this to discuss what makes a good variable name)
- Demonstrate `range()`

3. Demonstrate the use of loops to iterate through a file line-by-line (3:10-3:30)

- Start simple (just printing lines), then build up...
- print the lines...notice the extra whitespace...why?
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

for i in range(2000): 
    print(lines[i])

f.close()
```

- use `.strip()` to remove the whitespace
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

for i in range(2000):
  line = lines[i].strip()
  print(line)

f.close()
```

- lines are strings...can split them on whitespace
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  print(line_list[1])

f.close()
```

- the elements of the list are still strings...
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  print(type(line_list[1]))

f.close()
```

- convert them to numbers
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  start_coord = int(line_list[1])
  print(start_coord)

f.close()
```

- store in a list
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

start_coord_list = []

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  start_coord = int(line_list[1])
  start_coord_list.append(start_coord)

f.close()

print(start_coord_list)
```

- operate on the list
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

start_coord_list = []

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  start_coord = int(line_list[1])
  start_coord_list.append(start_coord)

f.close()

start_coord_list.reverse()

print(start_coord_list)
```

- compute the mean
```
f = open("mouseBed.bed", "r")

lines = f.readlines()

start_coord_list = []
sum = 0

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  start_coord = int(line_list[1])
  start_coord_list.append(start_coord)
  sum = sum + start_coord

f.close()

print(sum / len(start_coord_list))
```

- import functions from an outside library
```
import numpy

f = open("mouseBed.bed", "r")

lines = f.readlines()

start_coord_list = []

for i in range(2000):
  line_string = lines[i].strip()
  line_list = line_string.split()
  start_coord = int(line_list[1])
  start_coord_list.append(start_coord)

f.close()

start_coord_mean = numpy.mean(start_coord_list)
print(start_coord_mean)
```

Give students time to practice writing and using for loops (SWCarpentry exercises; 3:30-3:45)

## Homework (4:00-5:00)

Exercise prepared by Andrew and Matthew
