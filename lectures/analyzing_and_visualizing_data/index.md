## Motivation (9:00-9:05; PowerPoint)

- State learning objectives

# Software Carpentry: Analyzing Patient Data (Numpy; 9:05-10:00)

- [https://swcarpentry.github.io/python-novice-inflammation/instructor/02-numpy.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/02-numpy.html)

# Software Carpentry: Visualizing Patient Data (Matplotlib; 10:00-11:00)

- [https://swcarpentry.github.io/python-novice-inflammation/instructor/03-matplotlib.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/03-matplotlib.html)

# Exercise (11:00-12:00)

## Intro (live-coding; 11:00-11:15)

- Describe goal: recreate Fig 3B from Lott et al 2011 PLoS Biology [https://pubmed.gov/21346796](https://pubmed.gov/21346796)
- Discuss `plot-sxl.py` (i.e., “Who can tell me what this line is doing?”)
- Have students independently run code and confirm initial plot created

```
#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Get dataset to recreate Fig 3B from Lott et al 2011 PLoS Biology https://pubmed.gov/21346796
# wget https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/bulk_RNA-seq/extra_data/all_annotated.csv

# Overview of genfromtxt() https://numpy.org/doc/stable/user/basics.io.genfromtxt.html
data = np.genfromtxt( "all_annotated.csv", delimiter=",", names=True, dtype=None, encoding=None )

# Find transcript of interest
transcript = []
for row in data:
    if row['t_name'] == "FBtr0331261":
        print( row )
        transcript = row

# Obtain x and y values
x = []
y = []
for col in transcript.dtype.names:
    if "female" in col:
        x.append( col )
        y.append( transcript[col] )

# Plot data
fig, ax = plt.subplots()
ax.set_title( "FBtr0331261" )
ax.plot( x, y )
fig.savefig( "FBtr0331261.png" )
plt.close( fig )
```

## Extensions for student exercise (11:15-12:00)

- Annotate Plot
- Rotate x-axis labels
- Rename title to: "Sxl (FBtr0331261)"
- Add x- and y-axis labaels
- Create plot-msl-2.py for a different gene (FBtr0077571)

## Additional exercises (optional / advanced)

- Obtain male data (HINT: .startswith( "male" ))
- Plot female as red, male as blue
- Add legend
