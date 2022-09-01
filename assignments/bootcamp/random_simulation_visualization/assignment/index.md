# QBB2022 – Day 4 Homework: Visualizing Simulation Results and Power

## Instructions

* Create a `~/qbb2022-answers/day4-homework` directory.
* Create a working copy of the script `~/cmdb-quantbio/assignments/bootcamp/random_simulation_visualization/slides_asynchronous_or_livecoding_resources/binomial_power_interactive_lecture.py` within your new directory.
* Complete the single homework exercise and submit your edited code, the figure(s), and your `README.md` file to your GitHub repository.

## Homework Exercise -- Plot a heatmap, visualizing power for unfair coin tosses

In the interactive lecture, we computed the power for two unfair coin toss experiments: One experiment used a 0.6 probability of a coin toss to return heads and 500 coin tosses. Another experiment used a 0.95 probability of a coin toss to return heads, but only 10 coin tosses.

You will be editing and adding code to the simulation we wrote within the interactive lecture. Most of your code additions will be within the `run_experiment()` function. You may choose to write a function for plotting, but this is not required and the plotting code may be included within the `run_experiment()` function if you'd prefer.

### A -- Expand simulation conditions
For this exercise, within the `run_exerpiment()` function, expand your code to use a nested for loop to iterate through a range of probabilities controlling how unfair the coin is and to iterate through a range of the number of coin tosses within a simulation.
  * Please use `tosses = numpy.array([10, 50, 100, 250, 500, 1000])` as the range of the number of coin tosses.
  * Please use `probs = numpy.around(numpy.arange(0.55, 1.05, 0.05), decimals=2)[::-1]` as the range of probabilities of the unfair coin to return heads.
    * Take a moment to use manual pages and `print()` statements to break this line down, recording observations and conclusions in your `README.md` file .
    * First, focus on the inside function `numpy.arange(0.55, 1.05, 0.05)`. Based on [the manual page](https://numpy.org/doc/stable/reference/generated/numpy.arange.html), what are 0.55, 1.05, and 0.05 being used for? Use `print()` to look at the resulting array. Comment in your `README.md` on what you observe.
    * Second focus on the outside function `numpy.around( ,decimals=2)`. Based on [the manual page](https://numpy.org/doc/stable/reference/generated/numpy.around.html), what should it be doing? Use `print(numpy.around(numpy.arange(0.55, 1.05, 0.05), decimals=2))` to see how the `numpy.arange(0.55, 1.05, 0.05)` result changes when it is passed to `numpy.around()`. Comment in your `README.md` file on what you observe and why we might want to do this.
    * Finally, focus on the `[::-1]` part of the code. Use a `print()` statement to see how the array has changed. Describe what you see in your `README.md` file.

### B -- Compute and store power
Again, within the `run_experiment()` function, expand your code to not only compute, but also store the power of each experiment.
  * Create a numpy 2-dimensional array (or matrix) where you will store the power of experiment [(consider using `numpy.zeros` to initialize this array)](https://numpy.org/doc/stable/reference/generated/numpy.zeros.html)
  * Compute the power of each experiment (within your nested for loop)
  * Store the power within the 2-dimensional numpy array

### C -- Plot
Use Matplotlib and Seaborn to visualize the power for the unfair coin toss experiments within a heatmap.
  * Create a new `fig,ax` with `plt.subplots()`
  * Use `seaborn.heatmap()` to plot a heatmap of the power of the coin toss experiments.
    * [Use this manual page for `seaborn.heatmap()`](https://seaborn.pydata.org/generated/seaborn.heatmap.html)
    * Use the arguments `ax`,`vmin`, `vmax`, `cmap`, `xticklabels`, and `yticklabels`
      * `ax` is used to specify a matplotlib figure axis that seaborn should add the heatmap to.
      * Power is a probability. Therefore, there are absolute ranges we expect the values could take for any experiment. Use these limits to set `vmin` and `vmax`
      * We can set the colormap that the heatmap function will use to plot with the `cmap` argument. Look at these resource to read about picking a colormap, especially a color-blind friendly color palette.
        * [Seaborn color palette documentation](https://seaborn.pydata.org/tutorial/color_palettes.html)
        * [Stack overflow question about the topic](https://stackoverflow.com/questions/35735470/safe-seaborn-theme)
      * Use your array that stores the range of number of coin tosses (`tosses`) and the array that stores the range of probabilities of the unfair coin to return heads (`probs`) to set the `xticklabels` and `yticklabels` arguments such that the tick labels reflect the experiment conditions.
  * Use the normal matplotlib functions to set the x-axis and y-axis labels, the title, and save the plot to a file.
  * Within your `README.md` file, comment on any trends that you observe within your heatmap for the relation of power with the probability or number of toss parameters.
  * Make a plot that displays power, computed without Multiple Hypothesis Testing Correction and a plot that displays power with Multiple Hypothesis Testing Correction. It is up to you whether these are two separate plots or whether you have one plot with two columns/panels.

### D -- Compare to a real study
In your `README.md` file, comment on how this simulation relates to a real simulation performed recently by the McCoy lab  (specifically grad student Sara Carioscia).

  We briefly discussed the biology behind [this preprint](https://www.biorxiv.org/content/10.1101/2021.11.19.469261v2) within the interactive lecture.

  * Please Read the Abstract and Look at Figure S13 within the Supplementary Materials.
    * S13 is found on Page 17 of 18.
    * This figure is further described on page 2 of the Supplementary Materials as well as within the 'Strict adherence to Mendelian expectations across sperm genomes' Results section within the preprint.
  * If you still have questions about the motivation or data, consider reading the Introduction as well.
  * Within your `README.md`, briefly summarize the biological phenomenon of interest that this study is focusing on.
  * Within your `README.md`, compare the simulation experiment that you performed with the simulation experiment performed for Figure S13.
    * Generally, what parallels do you observe and what differences do you observe?
    * More specifically, which of the parameters within your experiment (probability of heads `prob_heads`, number of tosses `n_toss`, number of iterations `n_iters`) corresponds to the transmission rate axis?
    * Which of the parameters within your experiment corresponds to the Number of sperm axis?
    * Specifically, describe why both simulations use a binomial test.

### E -- Submit

Submit your code, the figure(s), and your `README.md` file to your GitHub repository.

## Advanced Exercise 1 -- Write a simulation for rolling a 6-sided die

Using the `simulate_coin_toss()` function from the interactive lecture as an example, write a function that simulates rolling a 6-sided die.

## Advanced Exercise 2 -- Write a simulation for rolling an n-sided die

Generalize your rolling a die simulation such that the function accepts as an argument how many sides the die is.

## Advanced Exercise 3 -- Use a new range of probabilities within a new unfair coin toss experiment

Add a new experiment to your coin toss experiment script. Consider a case where tails is a more likely result and the probability of heads ranges from 0 to 0.45. Re-visualize your heatmap. How does it compare to your original?
