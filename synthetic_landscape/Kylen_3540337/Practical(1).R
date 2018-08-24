# Kylen Brown 
# 18 August 2018
# Practical Assignment  

# Load required packages 
library(tidyverse)
library(vegan)
library(ggpubr)
library(gclus)
library(cluster)
source("coldiss.R")
source("cleanplot.pca.R")

# Load data 
env_prac <- as.tibble(read.csv("env.csv"))
spa_prac <- as.tibble(read.csv("spa.csv"))
spe_prac <- as.tibble(read.csv("spe.csv"))
spec_mass_prac <- as.tibble(read.csv("spe_mass.csv"))

# Question 1 
# dimensions given as rows x columns
dim(env_prac)
dim(spa_prac)
dim(spe_prac)

# Question 2 
# Location of each site 
mut_spa <- spa_prac %>% 
  mutate(N = 1:28)

plot1 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "darkblue") +
  geom_text(aes(label = N), colour = "red") +
  labs(title = "Location of each site", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
plot1

# ii
# a) Species Richness
# Compute the number of species at each site # baseR way
# To sum by rows, the second argument of apply(), MARGIN, is set to 1
site.pres <- apply(spe_prac > 0, 1, sum)
site_df <- data.frame(spa_prac = site.pres)
site.pres_df <- data.frame(spe_rich = site.pres)

spe_richness <- site.pres_df %>% 
  mutate(N = 1:28)

plot2 <-  ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "darkblue") +
  geom_point(aes(cex=spe_richness$spe_rich), shape = 1, colour = "red", show.legend = FALSE) +
  labs(title = "Species Richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
plot2

# b) Shannon-Weaver 
?diversity
spe_prac_divH <- as.tibble(diversity(spe_prac, index = "shannon"))

plot3 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "darkblue") +
  geom_point(aes(cex=spe_prac_divH$value), shape = 1, colour = "red", show.legend = FALSE) +
  labs(title = "Shannon-Weaver diversity index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
plot3

# c) Simpson's index
spe_prac_divsimp <- as.tibble(diversity(spe_prac, index = "simpson"))

plot4 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "darkblue") +
  geom_point(aes(cex=spe_prac_divsimp$value), shape = 1, colour = "red", show.legend = FALSE) +
  labs(title = "Simpson index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
plot4

combined1 <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2 )
combined1

# ii
# Species richness - Within the study area, the number of different species 
# differed in each quadrat due to the random dispersal of species from two 
# different points in the study area and the effects the of wind on dispersal.
# The number of different species, therefore, varied between sites (quadrats), 
# with many sites in the centre of the transects(4 and 5, 11 and 12 and 18 and 
# 19), having the greatest species richness. Site 7 and Site 1 , were the points 
# of dispersal for paper and gummy and jelly species respectively. This can 
# possibly explain the great species richness within the centre of the study 
# area as the effects of the wind caused species to be carried further away from
# their points of release into the centre of the study area. In the case of the 
# gummy and jelly species,they had a greater mass than the paper species and 
# naturally dispersed further away from their point of release than the paper 
# spcies, irrespective of the effects of the wind.Sites 27 and 28 have the 
# lowest species richness as species were dispersed towards the centre of the 
# study area from the two points of release, making it difficult for species to 
# end up at these sites (being in the outskirts of the study area) as well as 
# at other sites bordering the study area.

# Shannon - Weaver - The Shannon index emphasises evenness (how close in 
# numbers each species in an environment is). The eveness is high in the 
# quadrats in the centre of the study area due to these areas having the 
# greatest species richness and diversity. The high evenness in the centre 
# of the study area can be attributed to species being dispersed from their 
# points of release towards the centre of the study area. Sites on the 
# outskirts (sites 1 and 2, 27 and 29 and 22 and 23) have a lower species 
# richness and diversity and ,therefore, a lower evenness.  

# Simpson's index - The Simpson's index incorporates species richness and the 
# abundance of the relative species. The diversity of species closely 
# correlates to the species richness as seen in the "combined1" figure. The 
# diversity of species is greatest in the centre of the study area, the areas 
# possessing the greatest species richness. The species diversity is low at 
# sites 1 and 2 (site 1 being the release point of the gummy and jelly species), 
# whereas there is a higher species diversity at sites 6 and 7 (site 7 being 
# the release point of the paper species). As the gummy and jelly species 
# possess a greater mass, they dispersed further away from their point of 
# release than the paper species, irrespective of the wind. The wind also played
# a role in further dispersing the colonized paper species through the release 
# of the gummy and jelly species. As species richness and evenness increased, 
# so did the diversity. 

# Question 3 
# i 
spe_prac_t <- t(spe_prac)
spe_prac_t
spe_prac_t_S7 <- vegdist(spe_prac_t)
round(as.matrix(spe_prac_t_S7), 2)
# OR
cor(spe_prac)

# ii
coldiss(spe_prac_t_S7, byrank=FALSE, diag=TRUE)
# magenta = dissimilarity close to 0 = max similarity 
# cyan = dissimilarity close to 1 = min similarity 

# iii
# Individuals of the same species (same exact shapes) have dissimilarities of 0
# (maximum similarity). Big and small individuals of the same species have 
# varying dissimilarities, for example a small and big triangle have a 
# dissimilarity of 47 % whereas a small and big rectangle have a dissimilarity 
# of 82 %. Crumbled and uncrumbled species of the same shape have greater
# dissimilarities, for example a small crumbled triangle and a small triangle
# have a dissimilarity of 87 % and a small crumbled half circle and a small 
# half circle have a dissimilarity of 100 %. The crumbled shapes are greatly
# dissimilar to the gums and have dissimilarities of 1  (minimum similarity),
# for example a small crumbled triangle and the gums have a dissimilarity of 1 
# and a big crumbled square and the gums have a dissimilarity of 1 
# (minimum similarity). From an ecological context, individuals of the same
# species are similar to one another whereas individulas of different species
# are greatly dissimilar to one another. 

# Question 4  
# i 
spe_prac.pca <- rda(spe_prac, scale=TRUE) # scale = TRUE , mean = 0 , sd = 1
spe_prac.pca 
summary(spe_prac.pca) # Default scaling 2
# summary(spep.pca, scaling=1)

# Inertia explains the "variation" in the dataset and is calculated as the sum 
# of the diagonal of the correlation matrix that feeds into the PCA. It 
# represents the sum of all correlations of variables with themselves,
# correspnding to the number of variables (17 in this example). 
# The eigenvalues indicate the importance of the PC axes, the first PC axis has 
# the highest eigenvalue (3.791 in this example) and, therefore, possesses 
# the greatest explanatory power. The importance of the Pc axes diminishes until
# the smallest one is found for the last PC axis (0.029192 in this example). The
# sum of all the eigenvalues equals the total inertia. The eigenvalues can be 
# expressed as the Proportions explained (eigenvalue/total inertia). The 
# proportion of the variation that is explained for PC1 = 0.223, for 
# PC2 = 0.1456 and for PC3 = 0.1250 for example. The proportion of the variation
# that is cumulatively explained by the first three PC axes is 0.4936.
# Scaling is the way in which the ordination results are projected in 
# reduced space for graphical display. Scaling 1 should be used if the main 
# interest is to determine the relationships amoung objects (sites) whereas
# scaling 2 should be used if the main interest is to determine relationships 
# amoung descriptors (species).
# Species scores (response variables are always called species) refer to the 
# position of the arrow heads on the PCA biplot; they indicate the direction 
# and relative magnitude of species it is associated with as a linear gradient 
# across the 'landscape', i.e. the reduced space represented by the 
# corresponding PC axes. 
# Site scores (objects are always called sites) are the coordinates of the sites
# (rows).

# ii
par(mfrow = c(2, 2))
biplot(spe_prac.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe_prac.pca, main = "PCA - scaling 2 (species)") # default scaling 2
cleanplot.pca(spe_prac.pca, scaling = 1, mar.percent = 0.08) # sites scaling
cleanplot.pca(spe_prac.pca, scaling = 2, mar.percent = 0.04)

# As many axes maintained as there are principal components (dimensions/number 
# of variables). 

# Scaling 1 biplot - The first group of sites, site 6, 14, 17, 24 and 25 have 
# the greatest number of small half rectangles and gums with the numbers of 
# these species decreasing as you move further away from these sites. These 
# sites have the lowest numbers of all of the other species shapes. The third 
# group of sites, site 11, 12 and 13 have high numbers for a range of different
# species, however, have the lowest number of gummy species. The big and small
# crumbled half circles as well as the big square are greatest at sites 4 and 5
# with the values decreasing the further you move away from these sites. A range
# of other shapes are also correlated with both of these sites. 

# Scaling 2 biplot - It can be seen that the sites are organised in groups, with
# no descriptors being present in the top left quadrant of this ordination due 
# to no species being present at these sites. The lower left quadrant of the 
# biplot shows that small half circles and gums are highly positively 
# correlated with one another and are highly negatively correlated with the 
# species in the top right side of the biplot. The top right portion of the 
# biplot shows species which are highly correlated with each other.
# Species present in the top right portion of the biplot is highly correlated 
# with species in the bottom right quadrant of the biplot. The species 
# distributed on the right portion of the biplot can be observed as orthogonal
# arrows, indicating a correlation close to 0. For example, when looking at the 
# big half circles and the big crumbled half circles. Small half circles, big 
# squares and small trangles have shorther arrows and, therefore, are less 
# important for the ordination of sites in the ordination plane. These plots 
# demonstrate the importance of biplot representation as it summarises the main
# features of a data set and clusters of species are evident at specific sites. 

# Question 5 
# i
spe_prac_1 <- spe_prac[-16,]
remove_spe_prac <- spe_prac_1[-26:-27,]

spe_prac.nmds <- metaMDS(remove_spe_prac, distance="bray")
spe_prac.nmds
spe_prac.nmds$stress
dev.new(title="NMDS on species - Percentage difference")
plot(spe_prac.nmds, type = "t", main = paste("NMDS/Percentage difference", round(spe_prac.nmds$stress, 3)))

# Shepard plot and goodness of fit 
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe_prac.nmds, main="Shepard plot")
gof <- goodness(spe_prac.nmds)
plot(spe_prac.nmds, type="t", main="Goodness of fit")
points(spe_prac.nmds, display="sites", cex=gof*300)

# ii
# The Shepard diagram shows that the observed dissimilarities and the 
# ordination distances were up to 98% correlated (P,0.0001). The red step line
# reprents the fit of the data set.
# The Goodness of fit shows sites that are poorly fitted with larger bubbles. 
# The ordination plot shows that most sites are poorly fitted with site 10
# possessing the smallest circle (better fitted but not strongly fitted). 
# Crumbled paper species occur close to sites 4, 5, 12 and 20.
# No species are present close to sites 3, 13, 19, 21 and 22. The gummy species
# are not associated with any of the other species. Sites 7 and 8  have smaller 
# bubbles than the other sites, therefore, are positively fitted, however no 
# species occur at these sites(bordering the landscape). Sites 15 and 23 are 
# highly associated with small trangles. 

# Question 6 
#i
# Hellinger-transform the species dataset
spe_prac.hel <- decostand(spe_prac, "hellinger")
(spe_prac.rda <- rda(spe_prac.hel ~ ., env_prac))
summary(spe_prac.rda)	# Scaling 2 (default)
sum(spe_prac.rda$CCA$eig) # To calculate sum of inertia of CCA and eigenvalues 

# The overall variance is partitioned into constrained and unconstrained 
# fractions.The constrained fraction (0.05134) is the amount of variance of the 
# Y matrix and the unconstrained fraction yielded a value of 0.65514.
# The analysis yielded 2 canonical axes, eigenvalues labelled RDA1 to RDA2 and 
# 16 additional, unconstrained axes(labelled PC1 to PC16) for the residuals.
# The proportion constarined is 0.07267 in this example. 
# The canonical (RDAx) eigenvalues measures amounts of variance explained by the
# RDA model, whereas the residual (PCx) eigenvalues measure amounts of variance
# represented by the residual axes, but not explained by the RDA model l.
# The canonical eigenvalues RDA1 to RDA2 increase in value. The first
# residual eigenvalue (PC1) is, however, larger than the last canonical
# eigenvalue and is  larger than most RDA eigenvalues. The first residual 
# structure (axis),therefore, has the most variance. 
# Accumulated constrained eigenvalues: these are cumulative amounts of variance
# expressed as proportions of the total explained variance.The last cumulative 
# value is therefore 1.
# Species scoress are the coordinates of the tips of the vectors representing 
# the response variables in the bi- or triplots.
# Site scores (weighted sums of species scores): coordinates of the sites as 
# expressed in the space of the response variables Y.
# Site constraints (linear combinations of constraining variables):coordinates 
# of the sites in the space of the explanatory variables X.
# Biplot scores for constraining variables:coordinates of the tips of the 
# vectors representing the explanatory variables.
# Centroids for factor constraints:coordinates of centroids of levels of factor
# variables. 

# ii
# Global test of the RDA result
anova(spe_prac.rda, permutations=how(nperm=999)) 
# Tests of all canonical axes
anova(spe_prac.rda, by="term", permutations=how(nperm=999)) # p value for each location
anova(spe_prac.rda, by="axis", permutations=how(nperm=999)) # p value for each axis
# Variance inflation factors (VIF)
vif.cca(spe_prac.rda)

# The global test of the RDA results indicate that the model had a P-value of 
# 0.598, df of 2 and variance of 0.05134. The overall model fit relative to a 
# randomized or permutated matrix of data is examined by the first test. The 
# partial effects of the individual variables in the model is examined in the 
# second test. The P value indicates that there is no significant difference.
# After 999 permutations, an overall test of significance showed that canonical
# axes between matrices X and Y are insignificant. The canonical axes have p 
# values greater than 0.05 and are,therefore, all insignificant.  

# iii
# Scaling 1: distance triplot
dev.new(title="RDA scaling 1 + wa") 
par(mfrow = c(2, 2))
plot(spe_prac.rda, scaling=1, 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe_prac.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2 (default): correlation triplot
# dev.new(title="RDA scaling 2 + wa")
plot(spe_prac.rda, main="Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe_prac.rda, choices=1:2, display="sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length=0, lty=1, col="red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
# dev.new(title="RDA scaling 1 + lc")
plot(spe_prac.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2
# dev.new(title="RDA scaling 2 + lc")
plot(spe_prac.rda, display=c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length=0, lty=1, col="red")

# When species plot closely together, they share similar frequencies in their 
# abundance. In the triplot, only the carpet and cardboard environment is 
# present and,therefore, the wood environment is insignificant in explaining the 
# variation amoung species in the community. The rate of change of the carpet
# environmental factor is greater than the rate of change of the carpet
# environmental variable. The abundance of species is greater
# for the carpet environment than for the cardboard environment as seen in the 
# triplot. The bean species were the most abundant on the carpet environment and 
# the gums and the small half circles were the most abundant on the 
# cardboard environmet.   

# Question 7
# The landscape that was sampled consisted mainly of a flat carpet environment
# with wood and cardboard structures forming the steep incline and the peak of
# a mountainous region respectively. A stainless steel bowl was present in 
# the landscape and created and uninhabitable environment for species 
# colonisation. The distribution of species across the landscape was infuenced 
# by their points of dispersal, the force at which the species were dispersed, 
# the wind acting on species dispersal, the weight of the species in addition 
# to the suitability of each environment. Most species were present on the 
# carpet environment, with this environment being the most extensive. Few 
# species were present on the incline and peak of the mountainous region and no 
# species were present on the stainless steel bowl. As seen in this 
# investigation, most species were present at sites situated in the centre of 
# the landscape, this can be explained by the fact that species were dispersed 
# from their respective points of dispersal (two points of dispersal bordering 
# the landscape) towards the centre of the landscape. The wind acting on 
# species dispersal in addition to the force at which species were dispersed 
#allowed species to move further away from their points of release towards the 
#centre of the study area.
# The mass of the gummy and jelly species present in this landsacpe was greater
# than that of the paper species and, therefore, the gummy and jelly species 
# were dispersed futher away from their point of release towards the centre of 
# the study area than the paper species. As a result of this dispersal, there 
# were fewer species in the sites bordering the landscape. The number of species
# across the landscape also differed across the different environments due to 
# the suitability of each environment. There were few species present on the 
# mountain-like structure, indicating that this environment may be difficult
# to colonise and that only species possessing favourable traits are able to 
# exist in this region. Few species were found beneath the steep wood incline
# of the mountainous region, creating a microhabitat (differed from the 
# extensive surrounding carpet habitat) with specific environmental conditions,
# conducive for the species that were found there. No species were 
# present on the stainless steel bowl as this environment was unsuitable for 
# species colonisation. As explained above, many environmental variables are 
# responsible for the distribution of species across this landscape.     



