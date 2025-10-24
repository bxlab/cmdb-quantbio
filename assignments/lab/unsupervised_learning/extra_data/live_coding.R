remotes::install_github("allisonhorst/palmerpenguins")
library(tidyr)
library(dplyr)
library(matrixStats)
library(palmerpenguins)
data(package="palmerpenguins")

# Let's see how many NAs are in the data
penguins %>% summarize_all(funs(sum(is.na(.))))

# Let's see how the continuous variables relate to each other
penguins %>%
  select_if(is.numeric) %>%
  drop_na() %>%
  cor()

# Now let's remove the NAs and get rid of the non-numeric variables
valid_penguins = penguins[complete.cases(penguins),]
species = valid_penguins$species

# Since PCA doesn't handlle non-numeric data, let's convert our data to a usable format
summary(penguins$island)
num_penguins = valid_penguins %>%
  mutate(biscoe=as.numeric(island=="Biscoe")) %>%
  mutate(dream=as.numeric(island=="Dream")) %>%
  mutate(torgersen=as.numeric(island=="Torgersen")) %>%
  mutate(sex=ifelse(sex=="male", 0, 1)) %>%
  select_if(is.numeric)


# Now we will normalize the data
norm_penguins = num_penguins %>% scale()
colMeans(norm_penguins)
colSds(norm_penguins)

# Finally, let's run the PCA analysis
pca_results = prcomp(norm_penguins)
summary(pca_results)
pca_summary = tibble(PC=seq(1,ncol(norm_penguins),1),
                     sd=summary(pca_results)$sdev) %>%
  mutate(norm_var=sd^2/sum(sd^2)) %>%
  mutate(cum_var=cumsum(norm_var))

# Let's create a scree plot to see the variance explained
pca_summary %>% ggplot(aes(x=PC, y=norm_var)) +
  geom_line() +
  geom_point() +
  labs(y="Percent variance explained")
pca_summary %>% ggplot(aes(x=PC, y=cum_var)) +
  geom_line() +
  geom_point() +
  geom_text(aes(y=cum_var, label=round(cum_var, 2)), vjust=-1) +
  labs(y="Cumulative variance explained")

# Finally, let's see the PCA plot
PC_data = tibble(PC1=pca_results$x[,1],
                 PC2=pca_results$x[,2],
                 species=species,
                 sex=valid_penguins$sex)
PC_data %>% ggplot(aes(PC1, PC2, color=species, shape=sex)) +
  geom_point(size=3)

# Let's also look at the PC weightings
weights = tibble(x=pca_results$rotation[,1],
                 y=pca_results$rotation[,2],
                 name=rownames(pca_results$rotation))
r = max((weights$x^2 + weights$y^2)^0.5)
weights %>% ggplot(aes(x, y)) +
  geom_point() +
  geom_text(aes(x, y, label=name), vjust=-1) +
  annotate("path", x=r*cos(seq(0,2*pi,length.out=100)), y=r*sin(seq(0,2*pi,length.out=100))) +
  annotate("path", x=c(-r, r), y=c(0, 0)) +
  annotate("path", y=c(-r, r), x=c(0, 0))

# Time to move on to K-means clustering. We need to remove categorical variables for k-means
kmeans_penguins = scale(num_penguins %>% select(bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g))
kmeans_results = kmeans(as.matrix(kmeans_penguins), centers=6, nstart=100)

# Let's visualize the clusters using a heatmap
ordering = order(kmeans_results$cluster)
heatmap(as.matrix(kmeans_penguins)[ordering,], Rowv=NA, Colv=NA,
        RowSideColors=RColorBrewer::brewer.pal(12, name="Paired")[kmeans_results$cluster[ordering]], scale='none')

# And let's visualize the clusters using PCA
pca_results = prcomp(kmeans_penguins)
pca_tib = tibble(PC1=pca_results$x[,1],
                 PC2=pca_results$x[,2],
                 sex=valid_penguins$sex,
                 island=valid_penguins$island,
                 species=valid_penguins$species,
                 cluster=as.factor(kmeans_results$cluster))

# By species
ggplot(pca_tib, aes(PC1, PC2, color=cluster, shape=species)) +
  geom_point(size=3) +
  scale_color_discrete()
# By sex
ggplot(pca_tib, aes(PC1, PC2, color=cluster, shape=sex)) +
  geom_point(size=3) +
  scale_color_discrete()
# By sex
ggplot(pca_tib, aes(PC1, PC2, color=cluster, shape=island)) +
  geom_point(size=3) +
  scale_color_discrete()

