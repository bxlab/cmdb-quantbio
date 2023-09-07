1.) Without using any external libraries (such as numpy) write a function that takes a list (of any length) of integers as input and returns the mean (i.e., average). Write a script called `mean.py` where you create the list of integers, compute the mean using your function, and print it. 

2.) Use the function from part 1 to write a Python script called `mean_from_file.py` that computes and prints the mean of a series of integers from a data file (e.g., "my_integers.txt") where the data contain a set of integers, one per line. [Optional: Write your code so that the user can specify the name of that file from the command line (hint: use `sys.argv`).]

Extend the Day-1 homework without using any external libraries (such as numpy) unless otherwise noted to analyze data from `inflammation-01.csv`:

3.) Write a function that takes a patient row index (integer) as input and returns the mean inflammation level across the 40 days (float) for that given patient. Embed this function in a script called `mean_inflammation.py` that defines a patient row index as a variable, executes the function, and prints the output.

4.) Write a function that takes two patient row indices (integers) as input and returns a list of the difference between their inflammation levels on each of the 40 days (floats). Embed this function in a script called `difference_inflammation.py` that defines patient row indices as variables, executes the function, and prints the output.

Optional:

5.) Write a function that takes a patient row index (integer) as input and returns a dictionary with keys "mean", "min", and "max" and corresponding values representing the mean, minimum, and maximum of those integers (floats). Embed this function in a script called `mean_min_max_inflammation.py` that defines the patient row index as a variable, executes the function, and prints the output.

6.) Building on exercise 4, write a second function that uses matplotlib to produce a line plot of the difference between patient inflammation values (day on the x axis, difference in inflammation levels on the y axis). Have the plot appear on the user's screen. Embed this function in a script called `plot_difference.py` that defines the patient row indices as variables, executes the function, and generates the plot.

7.) Revise the function from exercise 6 to take an argument specifying the output file name for the plot. Rather than having the plot appear on screen, save the plot the the specified file path. Embed this function in a script called `plot_to_file_difference.py` that defines the patient row indices as variables, executes the function, and generates the plot.
