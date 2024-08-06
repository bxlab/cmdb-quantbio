library(readr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

data = readr::read_delim(args[1])

# Add log-transformed data column
data = data %>%
    dplyr::mutate(Log2_Expr=log2(Expr + 1))

# Add tissue-gene data column
data = data %>%
    dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", GeneID))

# Create plot
p = ggplot(data, aes(x=Tissue_Gene, y=Log2_Expr)) +
    geom_violin() +
    coord_flip() +
    xlab("Tissue + Gene") +
    ylab("Log2 Expression")

# Save plot
pdf(args[2])
print(p)
dev.off()

