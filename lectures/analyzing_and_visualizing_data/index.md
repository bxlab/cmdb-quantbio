# Motivation (9:00-9:05; PowerPoint)

- State learning objectives

# Conditionals and comparisons (live-coding; 9:05-9:20)
- [https://swcarpentry.github.io/python-novice-inflammation/instructor/07-cond.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/07-cond.html)
- everything before "Checking our Data"

# Software Carpentry: Analyzing Patient Data (Numpy; live-coding; 9:20-10:10)

- [https://swcarpentry.github.io/python-novice-inflammation/instructor/02-numpy.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/02-numpy.html)

- introduce concept of libraries, which allow access to functions that are not built in to python (9:20-9:25)
- use `numpy.loadtxt()` to load the inflammation dataset, then view it, store it in a variable and print it (9:25-9:30)
- discuss the concept of parameters to a function and order-based versus explicit specification (9:30-9:35)
- examine and discuss the `type()` (an n-dimensional array) (9:20-9:22)
- examine and discuss the `.shape` and use this as an example of an "attribute" of a python object (9:35-9:38)
- view all the attributes using the built-in `dir()` function (9:38-9:40)
- index various elements from the numpy array `data[0,0]` and discuss that these start from the top left (9:40-9:43)

- discuss and demonstrate data slices - start at the first index and go up to but not including the second index (9:43-9:48)
- show what happens if you dont include bone of the bounds on the slice `data[:3, 36:]`

- discuss and demonstrate built-in functions such as `numpy.mean(data)` (9:48-9:52)
- demonstrate some other numpy functions and the concept of multiple assignment

```
maxval, minval, stdval = numpy.amax(data), numpy.amin(data), numpy.std(data)
```

- what if we want to compute a statistic on only a subset of the data (e.g., a specific patient?) (9:52-9:58)
```
patient_0 = data[0, :] # 0 on the first axis (rows), everything on the second (columns)
print('maximum inflammation for patient 0:', numpy.amax(patient_0))
```

- this can be done in one line instead
```
print('maximum inflammation for patient 2:', numpy.amax(data[2, :]))
```

- numpy allows you to use the `axis` parameter to apply a function to each row or column (9:58-10:00)
```
print(numpy.mean(data, axis=0))
print(numpy.mean(data, axis=1))
```

- conditionals, continued "Checking our data" (10:00-10:10)
- [https://swcarpentry.github.io/python-novice-inflammation/instructor/07-cond.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/07-cond.html)



# Software Carpentry: Visualizing Patient Data (Matplotlib; 10:10-11:00)

- [https://swcarpentry.github.io/python-novice-inflammation/instructor/03-matplotlib.html](https://swcarpentry.github.io/python-novice-inflammation/instructor/03-matplotlib.html)

- work thorugh SWC as written, but using the plotting conventions below:
```
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.plot(xdata, ydata)

plt.show()
```

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
