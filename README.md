This public repository hosts materials for the Quantitative Biology Bootcamp and Fall Lab course (2023). It contains all materials necessary to render the course website, as well as any materials students will need for lectures and assignments.

The website is located at [bxlab.github.io/cmdb-quantbio](https://bxlab.github.io/cmdb-quantbio/) and rendered from the root directory of the `bxlab/cmdb-quantbio` repository. Materials necessary for **website content rendering**:

* `index.md` (website homepage)
* `_layouts/default.html`
* `Gemfile`
* `webpages/` (Contains subdirectories for all course days/weeks. Each subdirectory contains an `index.md` that will be the rendered webpage linked in the navigation bar dropdown menus and within the homepage schedules)

**Assignment materials** are within subdirectories for the two following directories:

* `assignments/bootcamp/` for bootcamp sessions/assignments
* `assignments/lab/` for lab sessions/assignments

Each assignment subdirectory contains:

* an `assignment/` directory with an `index.md` file where the assignment itself should be written and stored
* an `slides_asynchronous_or_livecoding_resources/` directory with a `README.md` file where any resources should be listed and linked or added
* an `extra_data/` directory with a `README.md` file where data should be listed, including any data pre-loaded on the student laptops. Extra data not pre-loaded onto the student laptops may be added to these directories as well. 

**Student resource materials:**

* Scripts and Python modules: `resources/code/`
* Syllabi: `resources/syllabi/`
* Reference sheets for syntax, etc: `resources/references/`
* Plotting gallery: `resources/gallery/`
