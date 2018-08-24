# QE Practicle 
# Jesse Smith

# datasets include:
# species data (spe)
# spetial data (spa)
# environmental data (env)
# species mass (spe_mass)


# load libraries  ---------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(vegan)
library(gclus)

# Load additionnal functions
# (files must be in the working directory)
source("evplot.R")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")
source("coldiss.R")
source("panelutils.R")


# read data ---------------------------------------------------------------

env <- as.tibble(read.csv("env.csv"))
spe <- as.tibble(read.csv("spe.csv"))
spa <- as.tibble(read.csv("spa.csv"))
spe_mass <- as.tibble(read.csv("spe_mass.csv"))

# question 1 --------------------------------------------------------------

# dimensions given as rows x columns
dim(env)
dim(spe)
dim(spa)
dim(spe_mass)

# question 2 --------------------------------------------------------------

# site location maps 
mut_spa <- spa %>% 
  mutate(N = 1:28)

fig1 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_text(aes(label = N), colour = "blue4") +
  labs(title = "Site Location", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_classic()

fig1

# Compute the number of species at each site
sit.pres <- apply(spe > 0, 1, sum)

site_df <- data.frame(spa = sit.pres)

# Sort the results in increasing order
sit.pres_df <- data.frame(spe_rich = sit.pres)

# a) Species Richness at the site positions 
spe_richness <- sit.pres_df %>% 
  mutate(N = 1:28)

fig2 <- ggplot(spe_richness, aes(x = N, y = spe_rich)) +
  geom_step() +
  geom_text(aes(label = N), colour = "red") +
  labs( x = "Positions of sites", y = "Species Richness", 
        title = "Species Richness") +
  theme_classic()

fig2

# species richness on Cartesian coordinates 

fig3 <-  ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=spe_richness$spe_rich), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Map of Species Richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_classic()

fig3 

# b) Shannon-Weaver Diversity Index
H <- as.tibble(diversity(spe, index = "shannon"))

fig4 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=H$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Shannon-Weaver Diversity Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_classic()

fig4

# c) Simpson index
simp <- as.tibble(diversity(spe, index = "simpson"))

fig5 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=simp$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Simpson Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_classic()

fig5 

ggarrange(fig2, fig3, fig4, fig5)


# ii)
# Species richness: Majority of the sites have 3 or fewer species. 
# Although sites 17 and 18 differ in there environments, i.e. wood and cardboard compared to carpet, these sites had more than 3 species each.
# Site 5 was the point of species dispersal and thus had the most species present. Sites 27 and 28 were the furthest away from the point of dispersal and therefore had no species present.
# Wind was aplied from between sites 1 and 8 and therefore, species richness was not as abundant at these sites either. Once the dispersal position was changed, the effect of the wind helped species distribution to vary.
# 
# Shannon Weaver: A high species diversity was prevalent around the middles quadrats of each transect. 
# The first point of dispersal showed a high species diversity, whereas the second point of dispersal had a lower species diversity due to species having a greater mass (gums and beans), and being unaffected by wind. 
# Wind affected colonised species at the second point of dispersal.
# 
# Simpsons index: a high Simpson index value implies that there is a high species richness as well as species evenness. 
# Overall species diversity was fairly high at most of the sites,  with higher values of diversity was found at the sites of dispersal. 
# Lower values occur at sites 1,2,7 and 8, where wind was a major infuence. 


# question 3 --------------------------------------------------------------
# i)
spe.t <- t(spe)
spe.jacc <- vegdist(spe.t)
matrix <- as.tibble(round(as.matrix(spe.jacc), 2))

# ii)
dev.new(title="Percentage difference (Jaccard), raw data", width=10, height=5)
coldiss(spe.jacc, byrank=FALSE, diag=TRUE)

#iii)
# The shades of magenta range from values of between 0.5 - 0 indicating high levels of dissimilarity, i.e. very similar. 
# Shades of blue range from 0.5 - 1, showing low simmilarity between species (maximum dissimilarity). 
# Gums, having a higher body mass, are poor dispersers and therefore are poorly associated with any of the other species, and have high dissimilairty values. 
# Squares, triangles, crumbled squares and crumbled triangles are strongly correlated with eachother.
# Small triangles are strongly correlated with big triangles and big half circles are  correlated with big pink elephants.
# species mass affects their dispersal ability especially since wind was the main dispersal factor, therefore species are more strongly correlated with species of similar mass, and lighter species are more easily dispersed than heavier species such as gum and beans. 
# Individuals of the same species are similar to one another and have dissimilarities close to 1 and species in the same genus (same shapes of big and small sizes) also share dissimilarites close to 0 and 
# individuals of different species share dissimilarities close to 1 (minimum similarity).

# question 4  -------------------------------------------------------------

# i)
spe.pca <- rda(spe, scale = TRUE)

summary(spe.pca)

# ii)
# Plot eigenvalues and % of variance for each axis
ev <- spe.pca$CA$eig
dev.new(title="PCA eigenvalues")
evplot(ev)

# PCA biplots
par(mfrow = c(2, 2))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") 
cleanplot.pca(spe.pca, scaling = 1, mar.percent = 0.08) 
cleanplot.pca(spe.pca, scaling = 2, mar.percent = 0.04) 

# The principle components help identify which variables are most strongly correlated and influential with each component. 
# the first 6 axes shuld be retained as they explain 75.96% of the variation. 
# the first 7 axes on the other hand, all have means above the average mean of all the axes.
# Species are randomly distributed among sites, with only certain species present at certain sites. Many sites have no species present, and many species have significant abundances at sites 5,11 and 18 where arrows extend further than the circle in the cleanplot. 
# The total iniertia is the sum of the diagonal of the correlation matrix that feeds into the PCA the sum of the eigenvalues. The value for total inertia is 17 for this dataset. 
# Only 22% varience is explained by PC1. (eigenvalue/ total inertia *100). 
# The cumulative proportion: upto and including PC6 only explains 75.96% of the variance. 
# Species scores: makes reference to the position (direction and length) of the arrow heads. 
# The most important original variables that make up PC1 are b_c_square, b_c_triangle, b_c_half_circle, b_triangle, s_c_half_circle, s_square and s_c_circle. 
# Site scores: sites 4, 5, 11, 12, 19 are most similar in PC1. 
# The figure on the left makes use of sites scaling (scaling 1) so the sites are scaled by eigenvalues. Here the distances in multidimensional space are more accurately reflected on the graph plane so it shows the relationships between sites better. 
# Sites that share similar environmental characteristics plots closer together and the ones that are further apart on the graph are also more dissimilar in their actual multidimensional space.
# The circle of equilibrium on the bottom left panel shows that species (gums, big half circles and pink elephant) strongly influence the spread of sites. 
# The presence of species are generally associated with the spread of sites towards the left of the vertical zero line, although their influence is felt rather obliquely (not perfectly alligned with the direction of spread in sites). 
# This might indicate that another important variable is lacking that can account for this direction of spread,
# but certainly the presence of certain species can be used to explain some of the spread of sites in this direction.

# question 5 --------------------------------------------------------------
# i)
# removing rows without any data 
spe1 <- spe[-16,]
spe2 <- spe1[-26:-27, ]

species <- vegdist(spe2, method = "bray")
spe.nmds <- metaMDS(spe2, distance = "bray")
spe.nmds
spe.nmds$stress
dev.new(title = "NMDS on Species - Percentage difference")
plot(spe.nmds, type = "t", main = paste("NMDS/Percentage difference - Stress =", round(spe.nmds$stress,3)))

par(mfrow = c(1,2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main = "Goodness of Fit")
points(spe.nmds, display = "sites", cex=gof*300)

# ii)
# In the goodness of fit model on the right, poorly fitted sites have larger bubbles, and well fitted sites have smaller bubbles. Well fitted sites include sites 2,8,9,10,12, and 22. 
# Beans are closely related to sites 2,8,9, and 10. b_C_triangle and b_c_square are associated with site 12. 
# The high scattering along the the NMDS1 axis implies a non linear relatonship.
# This NMDS produces a low stress value of 0.135 indicating a fairly good representation in reduced dimensions, and we can interpret the graph with confidence.
# It also means that the pairwise dissimilarities are well represented in ordination space (there is little scatter around the red line in the Shepard diagram), 
# with a non-metric r-squared value of 0.982 and a linear fit r-squared value of 0.906.


# question 6 --------------------------------------------------------------
# i)
# Hellinger-transform the species dataset
spe.hel <- decostand(spe, "hellinger")
(spe.rda <- rda(spe.hel ~ ., env))
summary(spe.rda)	# Scaling 2 (default)
sum(spe.rda$CCA$eig) # To calculate sum of inertia of CCA and eigenvalues 

# The variance is partitioned into constrained and unconstrained proportions. Where 7.27% of the variance is explained by constrained analysis and 92.73% of the vaiation is explained by unconstrained analysis.
# The canonical (RDAx) eigenvalues measures amounts of variance explained by the RDA model, where RDA1 and RDA2 cumulatively explains 7.23% of the variance in the RDA model. 
# The residual (PCx) eigenvalues measure amounts of variance represented 
# by the residual axes, but not explained by the RDA model l. The first 6 PC axes explain 79.20% of the residual variance. 
# Accumulated constrained eigenvalues: these are cumulative amounts of variance expressed as
# proportions of the total explained variance, where RDA1 explains 20.20% and RDA2 explains 79.80% of the variance.
# Species scores are the coordinates of the tips of the vectors representing the response variables in the bi- or triplots. Influencing species in RDA1 include s_rect, b_half_circles, pink_elephants, beans and gums.
# Site scores (weighted sums of species scores): coordinates of the sites as expressed in the space of the response variables Y. similar sites include sites 2, 8, 9, 11, 14, 17 and 18. 

# ii)
# Global test of the RDA result
anova(spe.rda, permutations=how(nperm=999))
# Tests of environmental variables
anova(spe.rda, by="term", permutations=how(nperm=999))
# Tests of all canonical axes
anova(spe.rda, by="axis", permutations=how(nperm=999))
# Variance inflation factors (VIF)
vif.cca(spe.rda)

# The global test of the RDA result showed no significant difference in the entire RDA with a p-value of 0.616. 
# The tests of all the environmental variables also showed no signoficant difference in environment type, with carpet having a p-value of 0.081 and cardboard with a p-value of 0.960 and wood not producing a result. 
# The anova of the axes are also insignificant with p-values of 0.203 and 0.924 for RDA1 and RDA2 respectively.

# iii)
## Triplots of the rda results (wa scores)
## Site scores as weighted averages (vegan's default)
# Scaling 1: distance triplot
dev.new(title="RDA scaling 1 + wa")
par(mfrow = c(2, 2))
plot(spe.rda, scaling=1, 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2 (default): correlation triplot
# dev.new(title="RDA scaling 2 + wa")
plot(spe.rda, main="Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
# dev.new(title="RDA scaling 1 + lc")
plot(spe.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2
# dev.new(title="RDA scaling 2 + lc")
plot(spe.rda, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length=0, lty=1, col="red")

# the RDA triplot is displayed in two forms: the first being "wa" in which the site scores are weighted sum of species and the second,"lc" being the fitted site scores. The fitted site scores (lc), are strictly orthogonal linear combinations of the explanatory variables and the site scores that are weighted sums of species appear more robust to noise in the environmental variables. All replicate sites with the same combination
# of factor levels are represented on top of one another in the fitted site scores
# (lc) triplot. The "wa" triplot is preferable in that case because
# the sites are separated in the plot and their labels can be read.
# Eigen values of the first 3 axes are 0.04, 0.01 and 0.15.
# From Scaling 1 we see that cardboard is mostly associated with gums, s_half_circles, pink_elephants and b_half_circles. Whereas carpet is closely associated with  the remainder of the species except b_square and s_rect. 
# Scaling 2: the angle between cardboard and carpet is greater than 90 edgrees, indicating that these two environments are strongly uncorrelated. 


# Question 7 --------------------------------------------------------------

# Landscape ecology refers to the study of the relationships found among the species, ecological 
# processess in the environment and ecosystems. The landscape in this study consisted of 
# cardboard, wood and carpet. The carpet covered the entire landscape, whereas the wood and cardboard were only
# found at specific sites of the landscape. The wood and cardboard represented mountains and inclines within the study area. This was calculated as percentage cover.
# species were randomly distributed amoung sites and between 
# carpet, wood and cardboard environments (most species were present on carpet).
# Sites either contained no species, few species or were greatly abundant in their
# species composition. Sites with similar species frequencies were close to one another, 
# implying that there may be an underlying environmental variable controlling their abundance 
# at these sites. Sites with no species present, also grouped close to one another, therefore,
# another environmental factor may be controlling the absence of species at these sites. 
# Species abundances on wood and cardboard environments were low, indicating that these
# environments may be unfavourable.


# bonus question  ---------------------------------------------------------

mass.hel <- decostand(spe_mass, "hellinger")
(mass.rda <- rda(mass.hel ~ ., env))
summary(mass.rda)	# Scaling 2 (default)
sum(mass.rda$CCA$eig) # To calculate sum of inertia of CCA and eigenvalues 

# Global test of the RDA result
anova(mass.rda, permutations=how(nperm=999))
# Tests of environmental variables
anova(mass.rda, by="term", permutations=how(nperm=999))
# Tests of all canonical axes
anova(mass.rda, by="axis", permutations=how(nperm=999))
# Variance inflation factors (VIF)
vif.cca(mass.rda)

## Triplots of the rda results (wa scores)
## Site scores as weighted averages (vegan's default)
# Scaling 1: distance triplot
dev.new(title="RDA scaling 1 + wa")
par(mfrow = c(2, 2))
plot(mass.rda, scaling=1, 
     main="Triplot RDA mass.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(mass.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2 (default): correlation triplot
# dev.new(title="RDA scaling 2 + wa")
plot(mass.rda, main="Triplot RDA mass.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(mass.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
# dev.new(title="RDA scaling 1 + lc")
plot(mass.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="Triplot RDA mass.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2
# dev.new(title="RDA scaling 2 + lc")
plot(mass.rda, display=c("sp", "lc", "cn"), 
     main="Triplot RDA mass.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length=0, lty=1, col="red")

# RDA scaling 1: on the carpet environment, bean species are most likely to be found, compared to to the cardboard environment where mostly gums are found. These are the two species with the heaviest mass and have poor dispersal ability. The lighter species are not correlated with any environment as they have a good dispersal ability. 
# RDA scaling 2: as seen in the species data, carpet and cardboard environments are uncorrelated. more sites are correlated with the carpet environment where only site 18 is associated with cardboard.
