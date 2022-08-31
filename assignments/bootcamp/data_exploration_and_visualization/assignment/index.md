# QBB2022 -- Day 3 Homework: Principal Component Analysis

## Prepare `day3-homework` directory

a. Create a `/Users/cmdb/qbb2022-answers/day3-homework` directory.

b. Use TextMate.app to create a `README.md` file in the new `day3-homework` directory.

c. Add the following line to `README.md` and save the file. Record the code you used and answers to questions within `README.md`.

```
# QBB2022 - Day 3 - Homework Exercises Submission
```
d. Submit your answers to your `qbb2022-answers` repository as you work on the assignment and complete an exercise; do not wait until all of the exercises are complete. Use informative commit messages. For example:

```
git add README.md ex2_a.png
git commit -m "edited figure and answers for day 3 hw exercise 2"
git push
```
e. Verify that your files were successfully uploaded to your online repository on GitHub.  

## Exercise 1

* Using PLINK, perform principal components analysis on the genotype data stored in `/Users/cmdb/data/vcf_files/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz`. Look at the PLINK documentation / Google examples of how to perform PCA with PLINK.
* Record the code you used to perform PCA above in your README.md file.

## Exercise 2

* Using numpy, read in the principal component coordinates for each sample. (In linear algebra terms, these are the first and second eigenvectors of the covariance matrix, and are therefore stored in the `plink.eigenvec` output file).
* Create a scatter plot using the first vs. second principal component coordinates for your samples.  Your figure should have 2548 total data points.
* Repeat the step above, but plotting the first vs. third principal component.
* Label the axes of each figure as PC1, PC2, or PC3, as appropriate.
* Upload the above figures to your GitHub repo as `ex2_a.png` and `ex2_b.png`.
* Do you notice any structure among the points? What do you think this structure represents?

## Exercise 3

* Use the unix `join` command to intersect the `plink.eigenvec` file with the metadata for those same samples stored in `/Users/cmdb/data/metadata_and_txt_files/integrated_call_samples.panel`. Output this joined table to a new file. (Read the documentation for `join`! **Hint:** for `join` to work, both files need to be sorted on the field you will use to join and both files also have to have the same delimiter character.)
* Read the above file into a numpy array. 
* Use these data to color the plot of PC1 vs. PC2 according to population, superpopulation, and sex (3 separate plots).
* Add a legend to each plot to explain the colors.
* Add an informative title to each plot.
* Label the axes of each figure as PC1 and PC2, as appropriate.
* Upload the above figures to your GitHub repo as `ex3_a.png`, `ex3_b.png`, and `ex3_c.png`.

## Optional exercise

* Try plotting the first three principal components together on a 3D scatter plot.
* See https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
* NOTE: Never put a figure like this in a (2D) paper! 3D plots are only potentially useful if you can animate or rotate them.

## Optional exercise

* Using any of the the metadata files in `/Users/cmdb/data/metadata_and_txt_files/`, develop your own visualizations to convey a point of your choosing.
* For example, using data in `integrated_call_samples.panel`, you could create a "stacked bar plot" to compare the sex ratio within different superpopulations (see https://matplotlib.org/stable/gallery/lines_bars_and_markers/bar_stacked.html).
