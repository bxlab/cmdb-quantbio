# Project Work Guide

Quantitative Biology Lab consists of two main parts: 
Lab Topics that introduce you to fundamental concepts and applications in biological data science and 
Project Work where you take a computational project of your own choosing from start to finish.
For the Project Work you will work throughout the semester in teams of two to identify a topic, develop a strategy, and implement your strategy.
Instructors, TAs, and even your peers will provide feedback and suggestions along the way.

## Grading and Timeline

Project Work represents 20% of your final grade and is comprised of four main assessments of equal weight

- Proposal -- 1-2 pages in README.md
- Check-in #1 -- GitHub repository review
- Check-in #2 -- GitHub repository review
- Final Presentation -- Rotation Talk format

## Background

### Ideas

Some possible project ideas depending on your interests and prior experience

- Reproduce a figure in a manuscript of interest (prior research, current rotation)
- Deep dive into the contents of your-favorite-genome
- Genome assembly using short, long, Hi-C datasets
- Exploring variation across ancestry (human) or evolution (viral, immune)
- Comparing gene expression between conditions
- Selecting two or more tools and comparing their performance on a test dataset

### Past Group Projects

- [QB24 Group Projects](../resources/past_group_projects/group_projects_2024.html)

### Example CMDB Projects

Here are some example GitHub repositories for *large* projects undertaken by current and former CMDB students (who also TA'ed Quantitative Biology Lab!).
These are provided purely as inspiration; your project is not expected to have such a scope

- Dylan Taylor, Stephanie Yan -- https://github.com/mccoy-lab/MAGE
- Sara Carioscia, Kathryn Weaver, and Andrew Bortvin -- https://github.com/mccoy-lab/transmission-distortion
- Peter DeFord -- https://github.com/pdeford/strum_paper
- Gherman Uritskiy -- https://github.com/bxlab/metaWRAP

## Data Repositories

Here are some data repositories that may contain datasets of interest

- [SRA](https://www.ncbi.nlm.nih.gov/sra) -- NCBI Sequence Read Archive
- [GEO](https://www.ncbi.nlm.nih.gov/geo) -- NCBI Gene Expression Omnibus
- [1KGP](https://www.internationalgenome.org) -- The 1000 Genomes Project
- [CaeNDR](https://caendr.org) -- The Caenorhabditis Natural Diversity Resource
- [GTEx](https://gtexportal.org) -- The Genotype-Tissue Expression Project
- [recount3](https://rna.recount.bio) -- recount3: uniformly processed RNA-seq
- [cellxgene](https://cellxgene.cziscience.com) -- Single-cell datasets
- [IDR](https://idr.openmicroscopy.org) -- The Image Data Resource
- [Allen Cell](https://allencell.org) -- Allen Cell Collection
- [Napari List](https://github.com/napari/docs/blob/main/docs/further-resources/sample_data.md) -- Public sample image databases across various imaging modalities

## Project Organization

Review "Good enough practices in scientific computing" [[pubmed](https://pubmed.gov/28640806)] for guidance on how to manage your data, organize your code, set up your project, and more.
Each team should create a new GitHub repository to collaborate with an organization similar to https://carpentries-lab.github.io/good-enough-practices/05-project_organization.html#example

```
	...
	|-- README            <-- title, description, overview figure
	|-- requirements.txt  <-- software (and versions)
	|-- data/             <-- only small (<1 Mb), instructions on retrieving larger
	|-- doc/              <-- reports, slidedecks
	|-- results/          <-- only small (<1 Mb), instructions on retrieving larger
	|-- src/              <-- code
	...
```

## Deliverables

### Proposal

Submit your Project Proposal as a README.md document in a new GitHub repository:

- Pick a repository name that is short but descriptive of your topic e.g. ebola-genomes
- Create a single GitHub repository that your group will use to collaborate on this project
- Add a README.md document with your proposal e.g. github.com/username/project-short-title/README.md
    - NOTE: For this week, it's fine (recommended?) that only the owner edit the file

Proposals should be ~1-2 pages in length with the following information:

- Title (1 pt)
- Description of 50-100 words (2 pts)
- Example published figure (display via Markdown using `![]` syntax ) (1 pt)
- Datasets with IDs, links (1 pt)
- Software with versions, links (1 pt)
- 2-3 proposed goals, with the first doable in <10 hours and the other two being "stretch" (4 pt)

### Check-ins

Progress will be assessed through two Check-ins.
In addition to addressing prior feedback given through GitHub Issues, describe your progress and any struggles by adding a new Markdown document in your `doc/` directory e.g. project-name/doc/checkin-date.md

Your check-in should be ~1-2 pages in length with the following information

- How you've addressed prior feedback (2 pts)
- New progress since last submission (4 pts)
- Project Organization (2 pts)
- Struggles you are encountering and questions you would like advice on (2 pts)

### Peer Review (Extra credit)

Each person will work individually to peer review two separate projects.
Pick one project that is familiar and one project that is unfamiliar to you and spend ~10 min on each.
Create a GitHub Issue with your feedback that provides 2-3 comments on the following:

- Aspects of the project that you find fascinating
- Questions you have about the approach or results
- Suggestions on next steps, project organization, etc.

### Final Presentation

Present your work in 6-8 min, leaving 2 min for questions:

- Create a slidedeck in the Shared Drive by copying template or making your own
    - Be sure you see your slidedeck appear in the Shared Drive
- Post link to slidedeck in README.md document of your GitHub repository
- Split time equally among group members, focusing on your contributions
