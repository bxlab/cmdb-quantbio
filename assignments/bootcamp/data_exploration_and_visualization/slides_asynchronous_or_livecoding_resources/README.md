Link to any outside resources in this README

The first script used in the live interactive lecture was `scatter_plot.py`, it made some line plots, annotated the plot, and used a multi-panel plot by the end. Iterations of the code shown below

Iterations of the `scatter_plot.py` code

1st iteration -- plotting a line
```
#!/usr/bin/env python

import matplotlib.pyplot as plt

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

fig, ax = plt.subplots() # create a figure and axes

ax.plot(x, y)
plt.show()
```

2nd iteration to plot multiple series or lines
```
#!/usr/bin/env python

import matplotlib.pyplot as plt

fig, ax = plt.subplots() # create a figure and axes

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

x2 = [2, 4, 6]
y2 = [8, 64, 216]

ax.plot(x, y)
ax.plot(x2, y2)
plt.show()
```

3rd iteration to annotate the plot with a legend

```
#!/usr/bin/env python

import matplotlib.pyplot as plt

fig, ax = plt.subplots() # create a figure and axes

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

x2 = [2, 4, 6]
y2 = [8, 64, 216]

ax.plot(x, y, label = "x^2")
ax.plot(x2, y2, label = "x^3")
ax.legend()

plt.show()
```

4th iteration to use multiple panels, each with a single line rather than multiple lines in the same panel
```
#!/usr/bin/env python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(nrows = 2) # create a figure and axes
type(ax)
len(ax)

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

x2 = [2, 4, 6]
y2 = [8, 64, 216]

ax[0].plot(x, y, label = "x^2")
ax[1].plot(x2, y2, label = "x^3")
ax[0].legend()

plt.show()
```

Final iteration that saves the plot to a file
```
#!/usr/bin/env python

import matplotlib.pyplot as plt
fig, ax = plt.subplots(nrows = 2) # create a figure and axes
type(ax)
len(ax)

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

x2 = [2, 4, 6]
y2 = [8, 64, 216]

ax[0].plot(x, y, label = "x^2")
ax[1].plot(x2, y2, label = "x^3")
ax[0].legend()
fig.savefig("lineplot.png")
plt.close(fig)

```

The second script used in the live interactive lecture was `histogram.py` and initially it simulated data from a normal distribution and plotted it. This is the earlier commented out code. Later, it used data made with bcftools within a created bcftools conda environment

On terminal beforehand
```
conda create -n bcftools -y python=3 bcftools -c conda-forge
conda activate bcftools
cp ~/data/vcf_files/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz .
bcftools query -f '%AF\n' ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz > chr21_af.txt
conda deactivate
python histogram.py
```  

Then, we used `barchart.py` to look at the number of genetic functional elements on each chromosome

On terminal beforehand
```
cp ~/data/gtf_files/gencode.v41.annotation.gtf.gz .
gzcat gencode.v41.annotation.gtf.gz | grep -v "^##" | cut -f1 | sort | uniq -c > genes_per_chrom.txt
```

Finally, the python script `hardy_weinberg.py` was used. The goal is to for each variant, extract the alternative allele frequency (AC / AN), as well as the frequency of homozygotes for the reference allele (0|0), homozygotes for the alternative allele (1|1), and heterozygotes (0|1 or 1|0). then make a scatter plot of alternative allele frequency versus frequency of each genotype class

On terminal beforehand
```
cp ~/cmdb-quantbio/assignments/bootcamp/annotating_and_writing_variants/slides_asynchronous_or_livecoding_resources/vcfParser.py .
cp ~/data/vcf_files/random_snippet.vcf .
```

The script was built such that we

* imported needed libraries/modules
* initialized some storage lists
* loaded a parsed vcf file whose name would be passed to the script from the command line.
* set up a for statement to loop through the parsed vcf list (skipping the header by starting the `range()` function at 1), and grab a single line or list from the vcf file, storing it in `snp`
* keep editing the for loop to extract the allele frequency from that specific `snp` line, (the 7th element of the SNP list is a dictionary, and we can use the key “AF” to get the allele frequency value for that SNP. We’re going to store or append it in the `af_list` storage list we initialized earlier
* Keep editing the for loop to count the specific genotypes we see across the 2000 some samples for that snp. We’ll use a nested for loop, initializing counters before the second/nested for loop statement. Then in the second for loop, we’ll look at the 10th element through the last element, adding to the counters based on the genotype we observe. The genotype is `snp[j]`
* Now, we divide the counts by the number of samples to get the frequencies, and then we append those to the lists. This is done outside of the second for loop, but within the first loop still. We’ve iterated through all of the samples, but we want to add these counts for every SNP in the vcf list
* Finally, outside of both for loops, we’re going to plot!
  * plot the observed data as three scatter series (one for each genotype, homozygous ref, heterozygous, homozygous alt)
  * plot the predicted/theoretical data from the Hardy Weinberg equation expectations as 3 lines plots
  * Note that mathematically in python:
    * `x**2` will take a number (x) and raise it to the 2nd power
    * `2 * x * (1-x)` will for a number (x), subtract x from 1, then multiply that difference with x and 2 for the final product.
    * The single `*` is used for multiplication here, and the parentheses are used to make sure that x is subtracted from 1 before any multiplication happens
