## Visualizing metadata from the GTEx Project

A description of GTEx is copied from the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex) below.

The Adult Genotype Tissue Expression (GTEx) Project is a comprehensive public resource to study human gene expression and regulation, and its relationship to genetic variation across multiple diverse tissues and individuals.

The project collected samples from up to 54 non-diseased tissue sites across nearly 1,000 deceased individuals. All individuals were densely genotyped to assess genetic variation within their genomes by Whole Genome Sequencing (WGS). Gene expression of each tissue was assessed by RNA sequencing (bulk RNA-seq). Expression quantitative trait loci (eQTLs) were identified as genetic variants that were significantly correlated with changes in the expression of nearby genes. The project provides a comprehensive identification of tissue-shared and tissue-specific human eQTLs, as well as a valuable basis for the mechanistic interpretation of the many non-coding genetic variants that have been associated with common human diseases, such as heart disease, cancer, diabetes, asthma, and stroke.

## Assignment

Building off of your work this morning, we will visualize various aspects of the sample metadata from the GTEx project.

Please submit your answers as a single `.R` script with comments that separate out answers to each question below. (Remember that comment lines start with `#` and are ignored by R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

Q1.  Load the `tidyverse` package, and use the function `read_delim()` to read in sample-level metadata that was obtained from the GTEx Portal (`GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`). In addition, open the data dictionary `GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx` in Excel, which provides a description of each column in the `.txt` file.

Q2.  View the first rows of the tibble by simply entering the variable name in which you stored it. Notice that some of the columns were cut off due to the limits of the display. Use the `glimpse()` function to examine the data types and first entries of all of the columns.

Q3.  Use the `filter()` function to subset the dataset to only the RNA-seq data by selecting rows for which the `SMGEBTCHT` column contains the value `"TruSeq.v1"`. [TruSeq](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html) is a library preparation kit from Illumina.

Q4.  Plot the number of samples from each tissue (`SMTSD`) as a barplot. (Hint: if you do not specify a y-axis, ggplot will use `stat = count` as the default, so the y-axis will represent the number of occurrences of each value of your x variable). See this [webpage](https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2) for a code snippet for rotating axis labels, which will be relevant throughout this exercise. Always be sure to label your axes with informative names!

Q5.  The RNA integrity number is a measurement of the degree of RNA degradation based on characteristics of an electropherogram trace. It ranges from 1 to 10, with 10 being the least degraded. Plot the distribution of RNA integrity numbers across your samples. What type of plot is best for visualizing a single continuous distribution? Take a look at this "[cheat sheet](https://images.datacamp.com/image/upload/v1666806657/Marketing/Blog/ggplot2_cheat_sheet.pdf)" for hints.

What is the shape of the distribution? Is it unimodal?

Q6.  Copy your code from above, but now plot the distribution of RIN, stratified by tissue. Consider what type of plot is best for contrasting continuous distributions across multiple groups.

Do you notice any differences across tissues? Are certain tissues outliers? What are your hypotheses to explain these observations?

Q7.  Visualize the number of genes detected per sample, stratifying by tissue. Again consider what type of plot is best for contrasting continuous distributions across multiple groups.

Do you notice any differences across tissues? Which tissues are outliers? Look over the abstract of this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7891839/) for one hypothesis to explain these observations.

Q8.  Plot the relationship between ischemic time and RIN. Consider what type of plot is best for visualizing the relationship between two continuous variables. Create sub-panels that stratify the data by tissue using `facet_wrap()`. Resize points to `size = 0.5` and set the opacity to `alpha = 0.5`. Add linear trend lines to your plot with \`geom_smooth(method = "lm").

What relationships do you notice? Does the relationship depend on tissue?

Q9.  Copy your answer from question 6 above, but modify it to color your points by autolysis score (`SMATSSCR`). Note that if we place `aes(color = SMATSSCR)` within the `ggplot()` portion of the code, it will attempt to apply this mapping to all `geom_`s, including `geom_smooth`. To avoid this, place `aes(color = SMATSSCR)` within the `geom_point()` portion of the code.

What relationships do you notice? Does the relationship depend on tissue?

Q10. If you finished early, make some more plots! What else can you learn about these data? Consider which type of plot visualizes a particular relationship most effectively. Keep it simple. Ideally, each figure should convey one main point. The purpose of a figure is to convey that point to the audience as clearly and concisely as possible.
