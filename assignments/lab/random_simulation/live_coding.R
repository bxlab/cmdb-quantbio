library(tidyverse)
library(palmerpenguins)

### Part 1: demonstration of bootstrap
# goal: compute the difference in median body mass between male and female Adelie penguins,
# as well as a confidence interval around this estimated difference in medians
# question: does the confidence interval overlap zero?

# 1a) prepare data

adelie_dt <- penguins %>%
  filter(species == "Adelie") %>%
  select(c("sex", "body_mass_g")) %>%
  drop_na()

# 1b) Point estimate (difference in medians)

obs_diff <- median(filter(adelie_dt, sex == "male")$body_mass_g) - median(filter(adelie_dt, sex == "female")$body_mass_g)


# 1c) Bootstrap (resample within each sex, with replacement)

bootstrap_replicates <- 5000
bootstrap_diffs <- rep(0, bootstrap_replicates) # create empty vector

set.seed(1)

for (i in 1:bootstrap_replicates) {
  male_b <- filter(adelie_dt, sex == "male") %>%
    slice_sample(n = nrow(.), replace = TRUE)
    
  female_b <- filter(adelie_dt, sex == "female") %>%
    slice_sample(n = nrow(.), replace = TRUE)

  d <- median(male_b$body_mass_g) - median(female_b$body_mass_g)

  bootstrap_diffs[i] <- d
}

# 1d) Compute 95% percentile CI

ci <- quantile(bootstrap_diffs, c(0.025, 0.975))
ci


### Part 2: demonstration of permutation test
# goal: compute a p-value for the difference in medians between male and female Adelie penguins

set.seed(2)

perm_replicates <- 10000
perm_diffs <- rep(0, perm_replicates)

for (i in 1:perm_replicates) {
  
  # shuffle the sex labels within the tibble
  adelie_perm <- adelie_dt %>%
    mutate(sex = sample(sex, replace = FALSE))
  
  # compute difference in medians using the same code as before
  d <- median(filter(adelie_perm, sex == "male")$body_mass_g) -
       median(filter(adelie_perm, sex == "female")$body_mass_g)
  
  perm_diffs[i] <- d
}

# One-sided p-value (male > female)
p_one_sided <- mean(perm_diffs >= obs_diff)
p_one_sided

# Two-sided p-value (optional)
p_two_sided <- mean(abs(perm_diffs) >= abs(obs_diff))
p_two_sided

# corresponding visualization
ggplot(aes(x = perm_diffs), data = tibble(perm_diffs)) +
  geom_histogram() +
  geom_vline(xintercept = obs_diff, linewidth = 1, color = "red") +
  theme_classic() + 
  xlab("Difference in medians (g)") +
  ylab("Count")

### Part 3: demonstration of power analysis

# https://www.researchgate.net/profile/Rui-Zhang-377/publication/262018045/figure/fig1/AS:296950378319872@1447809895486/Schematic-of-allele-specific-expression-A-The-two-chromosomal-copies-alleles-of-a.png


sample(c("R", "A"), size = 1, replace = TRUE)

sample(c("R", "A"), size = 5, replace = TRUE)

sample(c("R", "A"), size = 100, replace = TRUE) %>%
  table()

sample(c("R", "A"), size = 100, replace = TRUE, prob = c(0.9, 0.1)) %>%
  table()

binom.test(x = 55, n = 100, p = 0.5)

####

# power analysis

hist(rbinom(100000, size = 100, prob = 0.3), breaks = 0:100)
table(rbinom(100000, size = 100, prob = 0.3))

# simulate allele-specific expression (allelic ratio [prob] = 0.3)
# for 100,000 simulated sites with depth of 100 reads each, count how many times you see the ref allele
set.seed(1)
sim <- tibble(ref_reads = rbinom(100000, size = 100, prob = 0.3), total_reads = 100, prob_ref = 0.3) %>%
  rowwise() %>%
  mutate(pval = binom.test(ref_reads, total_reads, p = 0.5)$p.value) # compute p-values from binomial tests on each trial

power <- mean(sim$pval < 0.05) # power is the proportion of trials in which we rejected the null hypothesis
# note that we basically just need to set the number of trials high enough to ensure that our calculation of power
# doesn't depend on the number of trials that we ran

# now generalize to consider an arbitrary "fairness" of the coin and an arbitrary number of flips per trial
get_power <- function(p_ref, n_reads) {
  
  sim <- tibble(ref_reads = rbinom(1000, size = n_reads, prob = p_ref), 
                tot_reads = n_reads, prob_heads = p_ref) %>%
    rowwise() %>%
    mutate(pval = binom.test(ref_reads, tot_reads, p = 0.5)$p.value)
  
  power <- mean(sim$pval < 0.05)
  
  return(power)
  
}

power_sim <- expand_grid(prob_ref = c(0.5, 0.51, 0.55, 0.6, 0.7, 0.8, 0.9), 
                         num_reads = c(5, 10, 20, 50, 100, 200, 500, 1000)) %>%
  rowwise() %>%
  mutate(power = get_power(prob_ref, num_reads))

ggplot(data = power_sim, aes(x = prob_ref, y = power, color = factor(num_reads))) +
  geom_line() +
  ylab("Power to detect allele-specific expression") +
  xlab("Simulated allelic ratio (effect size)") +
  scale_color_brewer(palette = "Dark2", name = "Depth of coverage\n(total aligned reads)")

