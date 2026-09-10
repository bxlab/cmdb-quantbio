# =============================================================================
# BLOCK 2 (LIVE CODING) -- Data visualization with ggplot2
#
# You met ggplot yesterday: ggplot(), aes(), geom_point(), and modifying a plot.
# Today we build up the rest of the grammar, and then write functions that
# produce our figures.
# =============================================================================

# start just like before to load in the data

library(tidyverse)

penguins <- read_delim(
  "https://gist.githubusercontent.com/slopp/ce3b90b9168f2f921784de84fa445651/raw/4ecf3041f0ed4913e7c230758733948bc561f434/penguins.csv",
  delim = ","
)

# Every ggplot is the same set of decisions, added together with `+`:
#
#   ggplot(DATA, aes(MAPPINGS)) +   # which columns map to which visual dimensions
#     geom_*() +                    # what to draw
#     scale_*() +                   # how values become positions / colors
#     facet_*() +                   # how to split into panels
#     theme_*()                     # everything that isn't data
#
# One consequence worth noting now: a ggplot is an object. You can save it to a
# variable, add to it later, and return it from a function.

p <- ggplot(penguins, aes(x = bill_length_mm, y = flipper_length_mm)) +
  geom_point()

p
p + theme_minimal()


## histogram

# build up from parts, varying bins, etc.
ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_histogram(bins = 30, position = "dodge") +
  labs(
    x = "Body mass (g)",
    y = "Number of penguins"
  )

# read the messages R prints under the plot. Two things it tells you:
#   - it picked 30 bins for you, which is a choice you might want to control
#   - it removed rows with missing values, which is data you threw away

# the number of bins changes the story; try it
ggplot(penguins, aes(x = body_mass_g)) + geom_histogram(bins = 5)
ggplot(penguins, aes(x = body_mass_g)) + geom_histogram(bins = 30)
ggplot(penguins, aes(x = body_mass_g)) + geom_histogram(bins = 100)

# binwidth is an alternative way of controlling this
ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram(binwidth = 250) +
  labs(
    x = "Body mass (g)",
    y = "Number of penguins"
  )

# the `position` argument matters once you add a fill
ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_histogram(binwidth = 250, position = "stack")        # the default

ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_histogram(binwidth = 250, position = "identity", alpha = 0.5)   # overlaid

ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_histogram(binwidth = 250, position = "dodge")
# "dodge" shifts bars sideways within a bin, so a bar's position on the x axis
# no longer means exactly what the axis says. This could be misleading
# for a continuous variable.


## density plot
ggplot(penguins, aes(x = body_mass_g, color = species)) +
  geom_density() +
  labs(
    x = "Body mass (g)",
    y = "Density"
  ) +
  theme_minimal()

# each curve integrates to 1 separately, so this hides how many penguins are in
# each group, meaning that you need to check the counts before you interpret it
penguins %>% count(species)


## scatter
# also demonstrate scale_x_log10() and geom_abline() here
ggplot(penguins, aes(x = bill_length_mm, y = flipper_length_mm, color = species)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Bill length (mm, log10)",
    y = "Flipper length (mm, log10)",
    caption = "Gorman et al., 2014"
  ) +
  scale_color_brewer(palette = "Set2", name = "Species") +
  theme_minimal()

# One thing to know about log scales before the assignment: log10(0) is -Inf and
# log10 of a negative number is undefined, so ggplot drops those rows. There are
# no zeros in the penguin measurements, but there are a great many zeros in gene
# expression data. Watch what happens:

toy_expression <- tibble(
  gene = paste0("gene", 1:10),
  tpm  = c(0, 0, 0.5, 2, 10, 50, 0, 120, 3, 0)
)

ggplot(toy_expression, aes(x = tpm)) +
  geom_histogram(bins = 10) +
  scale_x_log10()
# "Removed 4 rows containing non-finite values" -- 4 of our 10 genes are gone,
# and they are exactly the genes that are switched off, which is often the
# interesting group. One imperfect way to address this is to add a pseudocount, 
# which at least makes the choice visible in the code:

ggplot(toy_expression, aes(x = tpm + 1)) +
  geom_histogram(bins = 10) +
  scale_x_log10() +
  labs(x = "TPM + 1 (log10 scale)")

# with tens of thousands of points, alpha is not enough and you want to bin
ggplot(penguins, aes(x = bill_length_mm, y = flipper_length_mm)) +
  geom_hex(bins = 20) +
  scale_fill_viridis_c() +
  theme_minimal()


## boxplot and violin plot
ggplot(penguins, aes(x = species, y = body_mass_g, fill = species)) +
  geom_boxplot() +
  labs(
    x = "Species",
    y = "Body mass (g)"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2", name = "Species") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(penguins, aes(x = species, y = body_mass_g, fill = species)) +
  geom_violin() +
  labs(
    x = "Species",
    y = "Body mass (g)"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2", name = "Species") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# a boxplot shows five numbers, so it cannot show that a distribution has two
# peaks. The violin can. Body mass is bimodal within species here, and the
# next plot shows why.

library(ggbeeswarm)

ggplot(penguins, aes(x = species, y = body_mass_g, color = species)) +
  geom_beeswarm() +
  labs(
    x = "Species",
    y = "Body mass (g)"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2", name = "Species") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# when the sample size allows it, showing the observations on top of a summary
# is usually the right call
ggplot(penguins, aes(x = species, y = body_mass_g)) +
  geom_violin(fill = "grey90", color = NA) +
  geom_beeswarm(aes(color = sex), alpha = 0.7, size = 1) +
  labs(
    x = "Species",
    y = "Body mass (g)",
    color = "Sex"
  ) +
  theme_minimal()


## facets

ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram() +
  facet_grid(species ~ .)

ggplot(penguins, aes(x = body_mass_g, fill = sex)) +
  geom_histogram() +
  facet_grid(species ~ island)

# note that the empty panels are informative: Gentoo only occur on Biscoe

# facet_wrap() when you are splitting on a single variable
ggplot(penguins, aes(x = bill_length_mm, y = flipper_length_mm)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ species, ncol = 3)

# scales = "free_y" lets each panel choose its own y range. Useful when groups
# differ a lot in size, but it means panel heights can no longer be compared by
# eye.
ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram(binwidth = 250) +
  facet_grid(species ~ ., scales = "free_y")


## bar plot
# note that here I am piping the output of previous tidyverse commands (a tibble)
# into ggplot rather than providing the data argument
penguins %>%
  group_by(species) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = species, y = n, fill = species)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    labs(
      x = "Species",
      y = "Number of penguins"
    ) +
  theme(legend.position = "none")


## line plot

penguins %>%
  group_by(island, species, year) %>%
  summarize(mean_body_mass = mean(body_mass_g, na.rm = TRUE),
            sd_body_mass = sd(body_mass_g, na.rm = TRUE)) %>%
  ggplot(aes(x = year,
             y = mean_body_mass,
             ymin = mean_body_mass - sd_body_mass,
             ymax = mean_body_mass + sd_body_mass,
             color = species)) +
    geom_point() +
    geom_line() +
    geom_errorbar(width = 0.3) +
    facet_grid(. ~ island) +
    labs(
      x = "Year",
      y = "Body mass (g) +/- Std. Dev."
    ) +
  theme_bw()


# save a plot
p1 <- penguins %>%
  group_by(species) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = species, y = n, fill = species)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    x = "Species",
    y = "Number of penguins"
  ) +
  theme(legend.position = "none")

ggsave("~/Downloads/my_plot.pdf", plot = p1, width = 4, height = 4)

## 15 min: writing functions that make plots

# Look back at the plots above and count how many times we typed
# theme_minimal(), scale_fill_brewer(palette = "Set2"), and
# theme(legend.position = "none").

# Since a ggplot is an object, a function can build one and return it -- exactly
# the same idea as the functions from this morning.

plot_trait_histogram <- function(values, trait_name, binwidth = NULL) {
  plot_data <- tibble(value = values)

  p <- ggplot(plot_data, aes(x = value)) +
    geom_histogram(binwidth = binwidth, fill = "steelblue", color = "white") +
    labs(
      x = trait_name,
      y = "Count"
    ) +
    theme_minimal()

  return(p)
}

# note that this takes a plain vector, which we pull out with $, so there is no
# new syntax to learn -- it is the same kind of function as count_above()

plot_trait_histogram(penguins$body_mass_g, "Body mass (g)", binwidth = 250)
plot_trait_histogram(penguins$flipper_length_mm, "Flipper length (mm)", binwidth = 2)
plot_trait_histogram(penguins$bill_length_mm, "Bill length (mm)", binwidth = 1)

# because the function returns a plot object, the caller can still add to it
plot_trait_histogram(penguins$body_mass_g, "Body mass (g)", binwidth = 250) +
  theme_classic() +
  labs(title = "All penguins")


# a function can also take a value to filter on
plot_one_species <- function(data, species_name) {
  species_data <- data %>%
    filter(species == species_name)

  p <- ggplot(species_data, aes(x = bill_length_mm, y = flipper_length_mm)) +
    geom_point(alpha = 0.7, color = "steelblue") +
    labs(
      title = species_name,
      x = "Bill length (mm)",
      y = "Flipper length (mm)"
    ) +
    theme_minimal()

  return(p)
}

plot_one_species(penguins, "Gentoo")
plot_one_species(penguins, "Adelie")
plot_one_species(penguins, "Chinstrap")


# and saving can be a function too, so that every figure gets the same treatment
save_figure <- function(plot, filename, width = 5, height = 4) {
  ggsave(filename, plot = plot, width = width, height = height, dpi = 300)
  return(filename)
}

p2 <- plot_one_species(penguins, "Gentoo")
save_figure(p2, "~/Downloads/gentoo.pdf")


# if you have several to make, a for loop could work
for (this_species in unique(penguins$species)) {
  p <- plot_one_species(penguins, this_species)
  save_figure(p, paste0("~/Downloads/", this_species, ".pdf"))
}
