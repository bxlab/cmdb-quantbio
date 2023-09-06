# QBB2023 - Day 2 - Homework

## Instructions

Add each of your answers to `/Users/cmdb/qbb2023-answers/day2-homework`.  Please `git push` after each update and do not wait until the end of the session e.g.

```
git add code.py
git add image.png
git commit -m "Add answer for exercise 1, step 1"
git push
# Confirm at github.com
``` 

## Exercises

1. Starting with [plot-sxl.py](../../../webpages/plot-sxl.md), create `plot-sisA.py` to visualize sisA (FBtr0073461) in a fashion similar to [Lott et al 2011 PLoS Biology](https://pubmed.gov/21346796) Fig 3A by adding elements in the following order.  After each step, push your code **and** plot to your git repository and check your repository using the https://github.com web interface (i.e. results in four separate commits).

    - Plot female data
    - Add male data
    - Add 2*male data (HINT: 2 * np.array( y ))
    - Annotate plot (generalize x-axis to 10 not female_10, add title, add x- and y-axis labels)

1. Modify `plot-sisA.py` (do not create a new file) to load the transcripts information using `open()` and a `for` loop rather than `np.loadtxt()`.  Remember that the first line is a header and should not be stored in the transcripts list.  Push just your code to your git repository and confirm at https://github.com that your code no longer uses `np.loadtxt()`.

1. Recover your original `plot-sisA.py` using `git checkout <commit> plot-sisA.py`.  Push your code to your git repository and confirm at https://github.com that your code once again uses `np.loadtxt()`.

## Additional Exercises

1. Create `plot-fig3a.py` that contains a for loop to create individual plots for FBtr0073461, FBtr0070073, and FBtr0077279.  Push your code and three plots.

1. Create `plot-fig3a-subplots.py` to plot all three transcripts onto a single multipanel plot as in Fig 3A (HINT: plt.subplots(3)).  Push your code and new plot.

1. Complete [git-game](https://github.com/git-game/git-game) and [git-game-v2](https://github.com/git-game/git-game-v2).

