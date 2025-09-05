# start just like before to load in the data

library(tidyverse)  

penguins <- read_delim(
  "https://gist.githubusercontent.com/slopp/ce3b90b9168f2f921784de84fa445651/raw/4ecf3041f0ed4913e7c230758733948bc561f434/penguins.csv",
  delim = ","
) 

## histogram

# build up from parts, varying bins, etc.
ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_histogram(bins = 30, position = "dodge") +
  labs(
    x = "Body mass (g)",
    y = "Number of penguins"
  )

## density plot
ggplot(penguins, aes(x = body_mass_g, color = species)) +
  geom_density() +
  labs(
    x = "Body mass (g)",
    y = "Density"
  ) +
  theme_minimal()

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

## facets

ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram() +
  facet_grid(species ~ .)

ggplot(penguins, aes(x = body_mass_g, fill = sex)) +
  geom_histogram() +
  facet_grid(species ~ island)

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
  
