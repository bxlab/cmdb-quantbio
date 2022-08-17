---
title: Python Cheat Sheet
layout: default
---

<!-- Make columns even -->
<style>
  td:first-child { width: 20% }
</style>

## Python Cheat Sheet

#### STRINGS
##### A sequence of unicode (Python3) or ASCII (Python2) characters

|:------------------------------------:|:------------------------------------------------------------------------------------------------------------------------|
| `myStr = "Hello World"`              | Strings are specified with double or single quotes                                                                      |
| `newStr = str(variable)`             | Convert the given variable to a string                                                                                  |
| `newLst = myString.split(sep)`       | Split a string into a list with elements divided by the given separator (the default separator is space)                |
| `myStr.startswith(substring)`        | Check whether a string begins with a given substring. Check end of string with endswith(). Returns either True or False |
| `stripStr = myStr.lstrip(substring)` | Strip the given substring from the start (left) of myStr. Use rstrip() to strip from the end of the string              |
| `upStr = myStr.upper()`              | Convert all characters in myStr to uppercase. Use .lower() for converting to lowercase                                  |
{:.table.table-striped}

#### FLOATS & INTEGERS
##### Floating point numbers and integers are two separate numeric types in Python. Integers are whole numbers while floats contain a decimal.

|:-----------------------:|:----------------------------------------------------------------------------|
| `myInt = 6`             | Integer values have no quotations and do not contain a decimal point        |
| `myFloat = 6.0`         | Float values have no quotations and must contain a decimal point            |
| `newFloat = float(int)` | Converts the given int to a float. Also converts numerical strings e.g. "6" |
| `newInt = int(float)`   | Converts the given float to an integer. Also converts numerical strings     |
{:.table.table-striped}

##### Numerical Operations <mark>(Warning: Performing arithmetic with ints will result in outputs rounded to the nearest integer)</mark>

|:--------:|:----------------------|
| `x + y`  | Sum of x and y        |
| `x - y`  | Difference of x and y |
| `x * y`  | Product of x and y    |
| `x / y`  | Quotient of x and y   |
| `x % y`  | Remainder of x and y  |
| `x ** y` | x to the power of y   |
| `abs(x)` | Absolute value of x   |
{:.table.table-striped}

#### LISTS
##### A one dimensional ordered array of Python objects that can be modified (mutable)

|:-----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `myList = ['a', 'b', 'c']`          | Create a new list called 'myList'. Use [] to create an empty list                                                                                                                    |
| `myList.append('new entry')`        | Append a new entry to the end of a list                                                                                                                                              |
| `list1.extend(list2)`               | Concatenate lists                                                                                                                                                                    |
| `newLst = list(variable)`           | Convert variable into a list                                                                                                                                                         |
| `myList[index]`                     | Retrieve an item from myList at the given index. You can change the item at the given index by setting it equal to a new value. An index of -1 will return the last item in the list |
| `myList.insert(index, item)`        | Insert an item into myList at the given index                                                                                                                                        |
| `del myList[index]`                 | Delete an item from myList at the given index                                                                                                                                        |
| `slicedList = myList[start:stop]`   | Create a new list from a slice of myList. Either start or stop can be ommitted (e.g. myList[:3] will give you a slice of myList up to index 3)                                       |
| `myList.sort(key=..., reverse=...)` | Sort myList in a specific order                                                                                                                                                      |
| `newStr = ''.join(myList)`          | Combine all the elements of a list into a single string. Each element will be separated by whatever value is within the quotation marks                                              |
{:.table.table-striped}

#### SETS
##### An unordered collection of unique elements that cannot be modified (immutable)

|:---------------------------:|:----------------------------------------------------------------------------------------------------------| 
| `mySet = set(['a','b','c']` | Create a new set named 'mySet'                                                                            |
| `unionSet = set1 | set2`    | Union operation. Returns a set containing all elements of set1 and set2                                   |
| `interSet = set1 & set2`    | Intersection operation. Returns a set containing only the elements present in both set1 and set2          |
| `diffSet = set1 - set2`     | Difference operation. Returns a set containing only the elements that are different between set1 and set2 |
{:.table.table-striped}

#### DICTIONARIES
##### A data structure consisting of a set of keys that each map to some value. Keys must be unique and immutable but values can be practically any data structure and repeat themselves

|:---------------------------------------------:|:------------------------------------------------------------------------------------------------------------------|
| `myDict = {'key1':'value1', 'key2':'value2'}` | Create a dictionary named myDict. {} will create an empty dictionary                                              |
| `myDict[key]`                                 | Access the value associated with the given key in myDict                                                          |
| `myDict['key'] = 'newValue'`                  | Create a new key:value pair in myDict or update an old key with a new value                                       |
| `myDict.keys()`                               | Returns all the keys within myDict. Use .values() to return all values and .items() to return all key:value pairs |
{:.table.table-striped}

#### PYTHON CONDITIONS AND IF STATEMENTS
##### Python conditions evaluate to either True or False. If statements allow for control over the execution of a block of code once a certain condition is met

|:-----------------:|:---------------------------------------------------------------------------------------------------------------------------------------|
| `x == y`          | Returns True if x and y are equal to each other                                                                                        |
| `x is y`          | Returns True if x and y are identical to each other (i.e. x and y are pointing to the same object)                                     |
| `x != y`          | Returns True if x and y are not equal to each other                                                                                    |
| `x < y`           | Returns True if x is less than y. Check if x is greater than y with x > y                                                              |
| `x <= y`          | Returns True if x is less than or equal to y. Check if x is greater than or equal to y with x >= y                                     |
| `if statements`   | Check whether a condition is True and if so, execute the associated block of code                                                      |
|                   | `if x == y:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print('x equals y')`                                                                          |
| `elif statements` | After an initial if statement evaluates as False another condition can be checked with an elif (else if) statement                     |
|                   | `if x==y:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print('x equals y')`<br>`elif x > y:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print('x is greater than y')` |
| `else statements` | If no previous if or elif statements evaluate to True, another block of text can be executed under control of an else statment         |
|                   | `if x == y:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print('x equals y')`<br>`else:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print('x does not equal y')`      |
{:.table.table-striped}

#### IMPORTING MODULES
##### Modules are files containing Python code. After installation Python functions, classes, and variables can be imported from these files and used in new projects. There are many useful Python Modules out there such as numpy and matplotlib and you can even write your own modules! Commands for importing modules are written at the start of a Python file

|:---:|:---|
| `import moduleName` | Grants access to the entirety of a modules contents. Functions, variables, and classes from the module must be preceeded by the module name (i.e. moduleName.functionName() |
| `import moduleName as alias` | Grants access to the entirety of a modules contents under a user defined alias. Functions, variables, and classes from the module must be preceeded by the alias (i.e. alias.functionName()) |
| `from moduleName import functionName` | Grants access to a specific part of a module. A function, class, or variable imported in this manner does not have to be preceded by the moduleName (i.e. functionName()) |
| `from moduleName import *` | Grants access to the entirety of a modules contents. Functions, classes, and variables imported in this manner do not have to be preceded by the moduleName. Warning: best practice is to be specific with import statements and the wildcard import should be avoided |
{:.table.table-striped}

#### LOOPS
##### A way of repeating a block of code over an iterable object, a range of numbers, or until a certain condition is met

|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
| `for i in range(0, 4):`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print(i, i*100))`<br><br>`Out:`<br>&nbsp;`(0, 0)`<br>&nbsp;`(1, 100)`<br>&nbsp;`(2, 200)`<br>&nbsp;`(3, 300)`                 | For loop through range of numbers                                                   |
| `myList = ['a', 'b', 'c']`<br>`for entry in myList:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print(entry.upper())`<br>`Out:`<br>&nbsp;`A`<br>&nbsp;`B`<br>&nbsp;`C`                           | For loop through entries in a list. You can use a for loop with any iterable object |
| `count = 0`<br>`while count < 4:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print(count)`<br>&nbsp;&nbsp;&nbsp;&nbsp;`count += 1`<br>`Out:`<br>&nbsp;`0`<br>&nbsp;`1`<br>&nbsp;`2`<br>&nbsp;`3` | While loop that stops when condition is no longer met                               |
{:.table.table-striped}

#### FUNCTIONS
##### Functions are blocks of code that accept arguments as input and return some value. After defining a function it can be called repeatedly with new arguments. This makes them great for tasks that are performed multiple times. They're also handy for organizing a script by breaking it up into smaller pieces that each complete a defined task

```
def simpleFxn():
	"""
	Example of a simple function with no arguments and a return statement. 
	If no return statement is specified the function will return None. 
	Text within triple quotations is called a docstring and is used for 
	documenting code (like hashtags)
	"""
	return 'Hello World'

# After defining a function it must be called to run
print(simpleFxn())

Out:
	'Hello World'
```
```
def subtract(x, y):
	"""
	Example of a subtraction function that takes two numbers as arguments and 
	finds the difference between them. This function returns None
	"""
	print('%s minus %s equals...' % (x, y))
	print(x - y)

newVal = subtract(5, 3)
# Since there is no return statement in subtract() this will print None
print(newVal)

Out:
	5 minus 3 equals...
	2
	None
```

##### Some useful functions built-in to Python

|:------------------------:|:----------------------------------------------------------------------------------------------------------------------------------|
| `objType = type(object)` | Returns the class type of the given object                                                                                        |
| `listLen = len(myList)`  | Returns the length (number of items) of the object. This object can be a sequence such as a string or a collection such as a list |
| `listSum = sum(myList)`  | Returns the sum of all numeric entries in a given list                                                                            |
{:.table.table-striped}

#### FILE I/O
##### Python has some built in tools for reading in and writing output to files

|:------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `f = open(filepath, mode)`                                                                | Creates a file object in Python. You must create this object before reading from the file. Mode describes the way the file will be used and is read only ('r') by default             |
| `for line in f:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`print(line)`                                 | Loops through each line in a file object. Looping is often the most efficient way to read files in Python                                                                             |
| `fileString = f.read()`                                                                   | Read the file object and return the file's contents as a string                                                                                                                       |
| `lineString = f.readline()`                                                               | Reads a single line of the file and returns that line as a string. Also progresses you forward one line in the file so the next time readline() is called the next line will returned |
| `f.write(string)`                                                                         | Writes the contents of the string to the open file and returns the number of characters written. This will only work if the file was opened in a write capable mode                   | 
| `f.close()`                                                                               | After you're done with a file you should **always** close it to prevent it from taking up system resources                                                                            |
| `with open(filepath) as f:`<br>&nbsp;&nbsp;&nbsp;&nbsp;`read_data = f.read()`             | Opening a file using the with keyword is good practice and will automatically close the file once the indented code block has completed                                               |
{:.table.table-striped}

#### PANDAS
##### pandas is a popular module in Python for working with dataframes (tables with rows and columns). Below is just the tip of the iceberg in terms of what you can do with pandas. For more information see the [pandas documentation](https://pandas.pydata.org/pandas-docs/stable/)

|:---------------------:|:-------------------------------------------------------------------------------------------------------------------------------------|
| `import pandas as pd` | To use pandas you have to import it at the beginning of your Python file. It's common practice to import pandas with the alias of pd |
{:.table.table-striped}

##### Creating dataframes

|:-----------------------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `df = pd.read_csv(filepath, sep=',')`                 | Create a dataframe from a csv file. Since csv files are comma separated the default sep is ',' but you can change this to any delimiter you want (e.g. use sep='\t' for tab separated files). There are other functions available in pandas for reading in other filetypes that go by the naming scheme pd.read_filetype() |
| `pd.DataFrame(array, headers=['my', 'headers'])`      | Create a pandas dataframe from a numpy array or list of lists                                                                                                                                                                                                                                                              |
| `df.columns = ['my', 'desired', 'column', 'headers']` | Set or change the dataframe column headers to the specified labels                                                                                                                                                                                                                                                         |
{:.table.table-striped}

##### Combining and adding to dataframes

|:---------------------------------------------:|:------------------------------------------------------------|
| `df['newColumn'] = ['new', 'column', 'data']` | Add a list as a new column to the dataframe                 |
| `df.concat([df1, df2], axis=1)`               | Combine two dataframes by rows (axis=0) or columns (axis=1) |
{:.table.table-striped}

##### Subsetting

|:-----------------------:|:--------------------------------------------------------------------------------------|
| `df[['my', 'headers']]` | Returns a df containing only the specified columns                                    |
| `df[df[column] != 0]`   | Filter a dataframe by row based on whether a condition is met in the specified column |
| `df.head(n)`            | Takes the first 'n' rows of the df. You can use df.tail(n) to take the last 'n' rows  | 
{:.table.table-striped}

##### Data selection

|:----------------------------------:|:------------------------------------------------------------------|
| `df.iloc[rowNumber, columnNumber]` | An integer position based way of selecting entries in a dataframe |
| `df.loc[rowName]`                  | A label based way of selecting rows from a dataframe              |
{:.table.table-striped}

##### Operations and inspecting data

|:------------------------:|:----------------------------------------------------------------------------------------|
| `df.sort_values(column)` | Sort a dataframe by values in the specified column                                      |
| `df.mean()`              | Returns the mean for all columns in the dataframe                                       |
| `df.max()`               | Returns the max value for all columns in the dataframe. Use min() for the minimum value |
| `df.std()`               | Returns the standard deviation of each column in the dataframe                          |
| `df.shape`               | Gives the number of rows and columns in the dataframe                                   |
{:.table.table-striped}

##### Saving dataframes to a file

|:------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `df.to_csv(outputFilename, sep=',')` | Writes the contents of a dataframe to a file. Again, the default delimiter here is a comma for csv but you can specify whatever delimiter you wish with the sep argument |
{:.table.table-striped}