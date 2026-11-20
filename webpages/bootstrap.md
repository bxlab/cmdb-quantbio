```
library(tidyverse)
library(palmerpenguins)

# Three uses of simulation:
#   1. Bootstrap: how uncertain is our estimate?
#   2. Permutation test: how surprising is our result if there is no effect?
#   3. Power analysis: how often would we detect an effect of a given size?


### Part 1: Bootstrap ##########################################################

# 1a) Prepare data
adelie_dt <- penguins %>%
  filter(species == "Adelie") %>%
  select(sex, body_mass_g) %>%
  drop_na()

# 1b) Point estimate: difference in median body mass (male - female)
obs_diff <- median(filter(adelie_dt, sex == "male")$body_mass_g) -
  median(filter(adelie_dt, sex == "female")$body_mass_g)
obs_diff

# 1c) Bootstrap: resample penguins with replacement, within each sex
set.seed(1)
n_boot <- 5000
bootstrap_diffs <- rep(0, n_boot)

for (i in 1:n_boot) {
  male_b <- filter(adelie_dt, sex == "male") %>%
    slice_sample(prop = 1, replace = TRUE)
  female_b <- filter(adelie_dt, sex == "female") %>%
    slice_sample(prop = 1, replace = TRUE)

  bootstrap_diffs[i] <- median(male_b$body_mass_g) - median(female_b$body_mass_g)
}

# 1d) 95% percentile confidence interval
ci <- quantile(bootstrap_diffs, c(0.025, 0.975))
ci

# 1e) Check against theory: bootstrap the regression coefficient
#     and compare with the OLS confidence interval
set.seed(1)
boot_coefs <- rep(0, n_boot)

for (i in 1:n_boot) {
  boot_dt <- slice_sample(adelie_dt, prop = 1, replace = TRUE)
  boot_coefs[i] <- coef(lm(body_mass_g ~ sex, data = boot_dt))["sexmale"]
}

quantile(boot_coefs, c(0.025, 0.975))
confint(lm(body_mass_g ~ sex, data = adelie_dt))


### Part 2: Permutation test ###################################################

# Under the null hypothesis, sex labels are unrelated to body mass.
# Shuffle the labels (without replacement) and recompute the difference in medians.
set.seed(2)
n_perm <- 10000
perm_diffs <- rep(0, n_perm)

for (i in 1:n_perm) {
  adelie_perm <- adelie_dt %>%
    mutate(sex = sample(sex, replace = FALSE))

  perm_diffs[i] <- median(filter(adelie_perm, sex == "male")$body_mass_g) -
    median(filter(adelie_perm, sex == "female")$body_mass_g)
}

# One-sided p-value (male > female)
p_one_sided <- mean(perm_diffs >= obs_diff)
p_one_sided

# Two-sided p-value
p_two_sided <- mean(abs(perm_diffs) >= abs(obs_diff))
p_two_sided

# Null distribution, with the observed difference in red
ggplot(tibble(perm_diffs), aes(x = perm_diffs)) +
  geom_histogram() +
  geom_vline(xintercept = obs_diff, linewidth = 1, color = "red") +
  theme_classic() +
  labs(x = "Difference in medians (g)", y = "Count")


### Part 3: Power analysis #####################################################

# Example: allele-specific expression. At a heterozygous site, are reads from the
# reference (R) and alternate (A) alleles balanced 50:50?
# Background figure:
# https://www.researchgate.net/profile/Rui-Zhang-377/publication/262018045/figure/fig1/AS:296950378319872@1447809895486/Schematic-of-allele-specific-expression-A-The-two-chromosomal-copies-alleles-of-a.png

# Simulating reads with sample()
sample(c("R", "A"), size = 1, replace = TRUE)
sample(c("R", "A"), size = 5, replace = TRUE)
sample(c("R", "A"), size = 100, replace = TRUE) %>% table()
sample(c("R", "A"), size = 100, replace = TRUE, prob = c(0.9, 0.1)) %>% table()

# Binomial test: are 55 ref reads out of 100 consistent with a 50:50 ratio?
binom.test(x = 55, n = 100, p = 0.5)

# rbinom() simulates the number of ref reads directly
hist(rbinom(100000, size = 100, prob = 0.3), breaks = 0:100)
table(rbinom(100000, size = 100, prob = 0.3))

# Power for one scenario: allelic ratio = 0.3, depth = 100 reads.
# Simulate many sites, test each against 0.5, and count how often we reject the null.
set.seed(1)
sim <- tibble(ref_reads = rbinom(100000, size = 100, prob = 0.3), total_reads = 100) %>%
  rowwise() %>%
  mutate(pval = binom.test(ref_reads, total_reads, p = 0.5)$p.value)

power <- mean(sim$pval < 0.05)  # proportion of simulated sites where we reject the null
power

# Generalize to any allelic ratio and read depth
# (more simulations = more precise estimate of power)
get_power <- function(p_ref, n_reads, n_sims = 1000) {
  sim <- tibble(ref_reads = rbinom(n_sims, size = n_reads, prob = p_ref)) %>%
    rowwise() %>%
    mutate(pval = binom.test(ref_reads, n_reads, p = 0.5)$p.value)

  mean(sim$pval < 0.05)
}

power_sim <- expand_grid(prob_ref = c(0.5, 0.51, 0.55, 0.6, 0.7, 0.8, 0.9),
                         num_reads = c(5, 10, 20, 50, 100, 200, 500, 1000)) %>%
  rowwise() %>%
  mutate(power = get_power(prob_ref, num_reads))

ggplot(power_sim, aes(x = prob_ref, y = power, color = factor(num_reads))) +
  geom_line() +
  labs(x = "Simulated allelic ratio (effect size)",
       y = "Power to detect allele-specific expression") +
  scale_color_brewer(palette = "Dark2", name = "Depth of coverage\n(total aligned reads)")

```