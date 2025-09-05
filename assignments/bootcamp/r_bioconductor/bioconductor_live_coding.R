library(tidyverse)
library(SummarizedExperiment)
library(airway)      # provides the 'airway' SummarizedExperiment

## 1) What is a SummarizedExperiment?
##    - A container for: assays (matrices), rowData (features/genes),
##      colData (samples), and metadata. 

## 2) Load the airway dataset
data("airway")    # loads object named 'airway'
airway            # print a summary

# Quick orientation
class(airway)
dim(airway)       # rows = genes, cols = samples
assays(airway)    # names of assay matrices (usually "counts")
rownames(airway)[1:5]
colnames(airway)

## 3) Inspect assays, colData, rowData (6–8 min)

# 3a) assay: count matrix (genes x samples)
counts <- assay(airway, "counts")
dim(counts)
counts[1:5, 1:5]

# 3b) sample annotations: colData (DataFrame)
colData(airway)
as_tibble(colData(airway)) %>% 
  head()

# 3c) feature annotations: rowData (genes)
rowData(airway) %>% 
  head()

# 3d) metadata (unstructured list)
metadata(airway)

## 4) operate on the contents of the SummarizedExperiment

# 4a) shorten the sample IDs
colData(airway)$sample_id <- str_replace(colnames(airway), "^SRR", "S")

# 4b) Compute library sizes (column sums of counts) and store as colData
lib_sizes <- colSums(counts)
# colData has to be a data.frame and not a tibble, so mutate will not work on it
# can add a new column "manually" as follows
colData(airway)$lib_size <- lib_sizes

# meanwhile, we could convert it to a tibble for our own downstream use, but this
# does not modify the original object
as_tibble(colData(airway)) %>% 
  select(sample_id, dex, lib_size)

# 4c) Plot library sizes by treatment (tidyverse + ggplot2)
as_tibble(colData(airway)) %>%
  ggplot(aes(x = dex, y = lib_size, color = dex)) +
    geom_point(size = 3) +
    labs(x = "Dexamethasone treatment", y = "Library size (total counts)") +
    theme_minimal() +
    theme(legend.position = "none")

# 4d) Subset the SummarizedExperiment by samples or by genes
#     (square brackets: [rows, cols] = [genes, samples])
# demonstrate logical indexing

# Subset to treated samples only

names <- c("Rajiv", "Fred", "Mike", "Sadhana", "Sneha", "Lance")
is_ta <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
names[is_ta]

treated_subset <- as_tibble(colData(airway)) %>%
  filter(dex == "trt")
treated_samples <- treated_subset$Run
airway_trt <- airway[, (colnames(airway) %in% treated_samples)]
dim(airway_trt)

# Subset to the first 1000 genes (just as a demo)
airway_1k <- airway[1:1000, ]
dim(airway_1k)

## 5) Add a normalized assay (log-CPM)
##    - logCPM = log2( (counts / library_size) * 1e6 + 1 )

# Compute CPM using current counts and library sizes
lib_sizes <- colSums(assay(airway, "counts"))
# don't worry about the function below for now
cpm <- sweep(assay(airway, "counts"), 2, lib_sizes, FUN = "/") * 1e6
logCPM <- log2(cpm + 1)

# Store as a new assay
assay(airway, "logCPM") <- logCPM
assays(airway)

# Quick sanity check: values and range
assay(airway, "logCPM")[1:5, 1:5]

## 6) Create a SummarizedExperiment from scratch

# Small toy count matrix (genes x samples)
toy_counts <- matrix(
  c(10,  50,  5,  0,
    200, 150, 1,  2,
    0,   0,   25, 40),
  nrow = 3, byrow = TRUE,
  dimnames = list(
    c("GeneA","GeneB","GeneC"),
    c("S1","S2","S3","S4")
  )
)

# Sample annotations (colData)
toy_col <- tibble(
  sample = c("S1","S2","S3","S4"),
  condition = c("ctrl","ctrl","trt","trt")
) %>% 
  column_to_rownames("sample")

# Feature annotations (rowData)
toy_row <- tibble(
  gene = c("GeneA","GeneB","GeneC"),
  biotype = c("protein_coding","protein_coding","lncRNA")
) %>% 
  column_to_rownames("gene")

toy_se <- SummarizedExperiment(
  assays = list(counts = toy_counts),
  colData = S4Vectors::DataFrame(toy_col),
  rowData = S4Vectors::DataFrame(toy_row),
  metadata = list(description = "Demo SummarizedExperiment")
)

toy_se
assays(toy_se)
colData(toy_se)
rowData(toy_se)
metadata(toy_se)

## 8) Save & reload

# Save your processed object to reuse later
saveRDS(toy_se, file = "~/toy_se.rds")

# use the broom to clear workspace

# Reload (e.g., next session)
toy_se_reloaded <- readRDS("~/toy_se.rds")
