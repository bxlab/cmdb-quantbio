library(readr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

data = readr::read_delim(args[1])

# Create plot
p = ggplot(data, aes(x=Min, y=Max)) +
    geom_point(alpha=0.2) +
    xlab("Min Log2 Expression") +
    ylab("Max Log2 Expression")

# Save plot
pdf(args[2])
print(p)
dev.off()

