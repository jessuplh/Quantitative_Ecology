# QuantAssign.R
# Adnaan Emandien
# 3540482
# 20 August 2018
# Quantitative Ecology Assignment

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(vegan)
library(cluster)
library(FD)
library(gclus)
library(ade4)

source("C:/Users/Adnaan/My Documents/NEwR_updated_material_R332-NEwR_ed1/coldiss.R")
source("C:/Users/Adnaan/My Documents//NEwR_updated_material_R332-NEwR_ed1/cleanplot.pca.R")

# Load datasets -----------------------------------------------------------

spa <- read.csv("Qspa.csv")
env <- read.csv("Qenv.csv")
spe <- read.csv("Qspe.csv")
mass <- read_csv("Qspe_mass.csv")

# Question 1 --------------------------------------------------------------

dim(spa)
dim(env)
dim(spe)
dim(mass)

# Question 2 --------------------------------------------------------------

# mutate data to create a row with site numbers
mut_spa <- spa %>% 
  mutate(ID = 1:28)

# first create a site map
site_map <- ggplot(data = mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "lightblue") +
  geom_text(aes(label = ID), colour = "red") +
  #  annotate("text", label = "Transect 1", x = 50, y = 0, colour = "red") +
  labs(title = "Study site", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
site_map


#i
#a) Species Richness

sit.pres <- apply(spe > 0,1,sum) # Sum by rows
sit.pres
aa <- as.tibble(sit.pres)

sort(sit.pres) # Sort into increasing numbers

ggplot(data = spa, aes(x = as.factor(x), y = as.numeric(y)))+
  geom_line(colour = "salmon") +
  geom_point(cex=5*sit.pres/max(sit.pres), shape = 1, colour = "red", size = 5) +
  labs(title = "Species richness", x = "x coordinate (cm)", y = "y coordinate (cm)")+
  theme_bw()

sit.pres.tib <- as.tibble(sit.pres)

sit.pres.mut <- sit.pres.tib %>% 
  mutate(N = 1:28)
sit_df <- as.data.frame(sit.pres.mut)

ggplot(sit_df, aes(x = N, y = value)) +
  geom_step(colour = "lightblue") +
  geom_text(aes(label = N), colour = "red") +
  labs( x = "Positions of sites", y = "Species Richness", 
        title = "Species Richness") +
  theme_bw()

# species richness on spatial
Fig1 <-  ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=sit_df$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Map of Species Richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()

# b) Shannon-Weaver

H <- as.tibble(diversity(spe, index = "shannon"))

Fig2 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=H$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Shannon-Weaver Diversity Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()
# c) Simpsons index

simp <- as.tibble(diversity(spe, index = "simpson"))

Fig3 <- ggplot(mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "deeppink4") +
  geom_point(aes(cex=simp$value), shape = 1, colour = "blue4", show.legend = FALSE) +
  labs(title = "Simpson Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_minimal()

# Question 3 --------------------------------------------------------------

# i)
spe.t <- t(spe) #t = transpose
spe.t[1:15, 1:15]
spe.t.S7 <- vegdist(spe.t, binary = TRUE)
round(as.matrix(spe.t.S7)[1:15, 1:15], 2)

bray_spe <- vegdist(spe.t)
bray_spe

# OR (No standardizing or transposing) 
cor(spe)

coldiss(bray_spe, byrank = FALSE, diag = TRUE)
# Pink = more similar (dissim close to 0)
# Blue = least similar (dissim close to 1)

# Question 4 --------------------------------------------------------------

# i)
cor(spe)
spe.pca <- rda(spe, scale = TRUE)
summary(spe.pca)

# ii)
par(mfrow = c(1, 2))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") 

# cleanplots
par(mfrow = c(1, 2))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") # default scaling 2
cleanplot.pca(spe.pca, scaling = 1, mar.percent = 0.08) # sites scaling
cleanplot.pca(spe.pca, scaling = 2, mar.percent = 0.04)


# Question 5 --------------------------------------------------------------

#i)
# remove the rows with NAs
sp.1 <- spe[-16,] 
rm_sp <- sp.1[-26:-27,]

# run bray-curtis on new dataset with removed sites
spe.bray <- vegdist(rm_sp, method = "bray", binary = TRUE)

# run nmds on bray data
spe.nmds1 <- metaMDS(spe.bray)

dev.new(title="NMDS on fish species - Percentage difference")
plot(spe.nmds1, type = "t", main = paste("NMDS/Percentage difference", round(spe.nmds1$stress, 3)))

par(mfrow = c(1, 2))
stressplot(spe.nmds1, spe.bray)
ordiplot(spe.nmds1, type = "t", display = c("sites"),
         main = "NMDS with site scores")
abline(v = 0, h = 0, lty = 3)

# dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe.nmds1, main = "Shepard plot")
gof <- goodness(spe.nmds1)
plot(spe.nmds1, type = "t", main = "Goodness of fit")
points(spe.nmds1, display = "sites", cex = gof*300)


# Question 6 --------------------------------------------------------------

spe.hel <- as.tibble(decostand(spe, "hellinger"))
spe.hel

spe.rda <- rda(spe.hel ~ ., env)
spe.rda

summary(spe.rda) 

# Canonical coefficients from the rda object
coef(spe.rda)

# Unadjusted R^2 retrieved from the rda object
R2 <- RsquareAdj(spe.rda)$r.squared
R2

# Adjusted R^2 retrieved from the rda object
R2adj <- RsquareAdj(spe.rda)$adj.r.squared
R2adj

# Triplots of the rda results (wa scores)
dev.new(title = "RDA scaling 1 + wa")
par(mfrow = c(2, 2))
plot(spe.rda, scaling = 1, 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length = 0, lty = 1, col = "red")

# Scaling 2 (default): correlation triplot
# dev.new(title = "RDA scaling 2 + wa")
plot(spe.rda, main = "Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices = 1:2, display = "sp")
arrows(0, 0, spe.sc2[, 1]*0.92, spe.sc2[, 2]*0.92, length = 0, lty = 1, col = "red")

# Triplots of the rda results (lc scores)
# dev.new(title = "RDA scaling 1 + lc")
# par(mfrow = c(1, 2))
plot(spe.rda, scaling = 1, display = c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1]*0.92, spe.sc1[, 2]*0.92, length = 0, lty = 1, col = "red")

# Scaling 2
# dev.new(title = "RDA scaling 2 + lc")
plot(spe.rda, display = c("sp", "lc", "cn"), 
     main="Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1]*0.92, spe.sc2[,2]*0.92, length = 0, lty = 1, col = "red") 

# permutations

# Global test of the RDA result
anova(spe.rda, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(spe.rda, by = "axis", permutations = how(nperm = 999))

anova(spe.rda, by = "term", permutations = how(nperm = 999))
# Variance inflation factors (VIF)
vif.cca(spe.rda)

# Apply Kaiser-Guttman criterion to residual axes
spe.rda$CA$eig[spe.rda$CA$eig > mean(spe.rda$CA$eig)]


# Answers -----------------------------------------------------------------


# Question 1
Dimensions:
# Species data - 28 rows, 17 columns.
# Environmental data - 28 rows, 3 columns.
# Spatial data - 28 rows, 2 columns.
# Species mass data - 28 rows, 17 columns.

# Question 2
# Species richness figure shows that the highest species richness is found closer to the centre of the landscape. This makes sense due to the fact that the 'species' was thrown towards the centre of the landscape when they first 
# invaded the land. 
# Shannon-Weaver index figure shows the species richness incorporated with the relative abundances, with the highest indices being towards the centre, with the highest Shannon-Weaver indices located between 100 and 300 cm along the 
# y-axis of the 'landscape'. 
# Simpson's index figure shows that the species evenness increases as you increase you go towards the centre of the environment and as you go up vertically along the different x-coordinates (50, 100, 150, 200 cm). It makes sense 
# for the species to be more even towards the centre as again the species 'immigrated' to this landscape with majority of the species numbers located in the centre of the landscape.

#Question 3
# According to the similarity and dissimilarity plots [where pink indicates more similar (dissimilarity closer to 0) and blue indicates more similar (dissimilarity closer to 1)], the 'gums' species displays more blue blocks than 
# any other species showing that they are less similar (more dissimilar) to any of the other species. This could be because the gum species was spread out more across the landscape and found on its own when sampling most of the time.

# Question 4
# i) Since the association matrix is a correlation matrix, the sum of the eigenvalues along the diagonals equals the number of 'species'. The total and unconstrained inertia both equals 17, which equals the sum of the eigenvalues and 
# number of 'species'. The eigenvalue, proportion and cumulative proportion of the principle components are shown, showing that the eigenvalues of PC1 and PC2 are 3.79 and 2.48 respectively and explains 22.3% and 14.6% of the variation 
# showing the highest percentages and thus the most variation. From looking at PC1 through to PC6, the total amount of explanation offered is 76%. Scaling 2 for species and site scores allows for species to be scaled proportional to eigenvalues. 
# The general scaling constant of scores is 4.63. The species and sites scores are plotted. With regards to the species scores, the loadings (scores/scaled eigenvectors) indicate the degree pf correlation between the original and new variables 
# (PC1, PC2, etc.), they indicate how much each of the original variables contribute to PC1, PC2, etc. PC1 is made up of uneven contributions from most of the original variables, 'gums' and 's_half' are less important. The sign of the eigenvector 
# indicates the polarity of the correlation between the original variable and new variable. Site scores: the (scaled) coordinates of the objects (sites).
# ii) 6 axes must be retained, because the cumulative proportion adds up to 0.76, explaining 76% of the variation. 
# In biplot - Scaling 1, Sites 1,2,3,8,9,10,13,15,16,21,22,23,26,27 and 28 are clustered very close together. This means that they're closely related when it comes to their distances. Their distances are approximations of their Euclidian distances in 
# multidimensional space. Sites 4, 5, 11, 12, 18, and 19 are quite far from the cluster in terms of distances on the biplot, which means that they're more different to the group. Sites 6, 14, 17, and 25 and quite close to the cluster with sites 7, 20 
# and 24 are quite close to the cluster. The sites in the cluster all have very low species composition. Site 5 has the highest values of s_c_half, while Site 18 has the highest value of gums. 
# In biplot - Scaling 2, in the bottom left part of the biplot, s_half, gums and s_rect are highly positively correlated and are very highly negatively correlated with b_c_rect and beans. Majority of the species were found on the right and side of the 
# biplot and this due to the wind vector (fanning the species) influencing the spread and distribution of the species as they immigrated to the landscapes. Species b_rect, beans and b_c_rect very highly correlated, and these are very negatively correlated with s_half, gums and s_rect. The species b_triangle and s_square positively correlated with b_half and pink_el. Species b_c_half and s_rect have nearly orthogonal arrows indicating correlation close to 0. Species b_c_strong and s-c_triangle; b_half and b_square; and s_c_half and s_triangle are the three sets of species closely correlated to one another. B_suqare and s_triangle has the shorter arrows, showing that it's less important for the ordination of the sites in the ordination plan. 

# Question 5 
# ii) A Shephard stress plot showing the relationship between observed dissimilarity and ordination distances (i.e. the distances on the final plot). The data is not well correlated has there are a few outliers with the data. The graph has a poor linear relationship with constant increasing and stasis periods of the graph. Majority of points on the graph are clustered between 0.3 and 0.8 observed dissimilarity between 0.1 and 0.5 ordination distance. 
# With regards to the Goodness of fit graph, sites 7, 8, 23, and 26 are all outliers, showing less similarity to the other sites. Groupings of sites 13 and 22, 11 and 5, 15 and 25, and 18 and 17 show that these groupings are quite similar to each other and species within those sites overlapping.

# Question 6
# i)Partitioning of variance: The overall variance (0.70648) is divided into constrained (0.05134) and unconstrained (0.65514) values, with constrained showing 7.267% and unconstrained showing 92.733% of the overall variance.
# Eigenvalues and their contributions to the variance: The analysis produces 2 canonical axes with eigenvalues labelled RDA1 and RDA2 and 16 additional, unconstrained axes for residuals (PC1 to PC16). The results give the eigenvalues, as well as the cumulative proportion of variance explained (for the RDA axes) or represented (for the residual axes PC1 to PC16). The last cumulative value is therefore 1. The cumulative contribution to the variance obtained by the two canonical axes is the proportion of the total variance of the response data explained by the RDA. It's the same value as the proportion constrained, 0.07267. RDA1 to RDA2 are decreasing in value, the first residual eigenvalue (PC1 = 0.1490), however, is larger than the last canonical eigenvalue (RDA2 = 0.01037). The canonical (RDAx) eigenvalues measure amounts of variance explained by the RDA model, whereas the residual (PCx) eigenvalues measure amounts of variance represented by the residual axes, but not explained by the RDA model.
# Accumulated constrained eigenvectors: These are the cumulative amounts of variance expressed as proportions of the total explained variance, as opposed to their contributions of the total variance described before.
# Scaling 2 for species and site scores:
#   Species are scaled proportional to eigenvalues. Sites are unscaled: weighted dispersion equal on all dimensions. General scaling constant of scores:  2.089857.
# The species scores of the coordinates of the tips of vectors representing the response variables in the variance. In RDA1, b_half had the highest score of 0.265460, and beans had the lowest of -0.148416. In RDA, b_square had the highest score of 0.1565875, and s_half had the lowest of -0.1229982. In the four residual axes (PC11 to PC14), in PC1, 2, 3, and 4 the highest scores were b_c_triangle (0.2849), b_c_traingle (0.1705), b_c_triangle (0.3005) and pink_el (0.2067) respectively, with the lowest scores belonging to beans (-0.7991), s_rect (-0.7495), b_square (-0.3088) and s_triangle (-0.3015) respectively.
# The site scores (weighted sum of species scores): The coordinates of the sites as expressed in the space of the response variables Y.
# Site constraints (linear combinations of constraining variables): The coordinates of the sites in the space of the explanatory variables X. These are fitted site scores.
# Biplot scores for constraining variables: The biplot scores for RDA, increases from -0.9987 to -0.05939 for carpet and the scores of PC1 to PC4 all equals 0. The scores for RDAs decreases from 0.7354 to -0.67760 for cardboard with the scores for PC1 to PC4 all equalling 0, the same as for carpet.
# ii)  When running the global permutation test of the RDA result, the Pr = 0.636 indicating that the model is not significant. When testing all canonical axes, the Pr for RDA1 and RDA2 are 0.178 and 0.940 respectively showing both canonical axes to be non-significant. The carpet and cardboard Pr scores of 0.069 and 0.966 are also both lower than 0.05 meaning that they are both non-significant. The variance inflation factors (VIF) measures the proportion by which the variance of a regression coefficient is inflated in the presence of other explanatory variables. Carpet has a VIF score of 1.928571 and cupboard has a VIF score of 1.928571 and Wood not being applicable. The Kaiser-Guttman criterion consists in computing the mean of all eigenvalues and interpreting only the axes whose eigenvalues are larger than the mean (in this case being PC1 to 7) with the PC's having the following scores: PC1 = 0.14898504, PC2 = 0.11246212, PC3 = 0.07949840, PC4 = 0.06204632, PC5 = 0.05426915, PC6 = 0.05094135 and PC7 = 0.0429224.
# iii) Triplot RDA spe.hel ~ env2 - scaling 1- wa scores: The angles between sites and species in the biplot reflect their correlations (but not the angles among sites). All sites (rows) are concentrated (clustered) around the point of origin (0,0). There are a few distinct outliers however as in the top right of the plot, row 17 and 23 are clearly seen along with b_square and s_rect, In the bottom right part, row 14 and 18 are clearly seen along with the gums, c_half and pink_half species, meaning that these species are a lot closer related than the rest of the others. Also located here is the vector for the variable cardboard. That means these sites and species and more closely related to this variable. Majority of the other species are found in the bottom left corner, strongly related to the variable carpet. The cluster of species here shows that they are close in terms of Euclidian distances. There are no species found in the top right corner. 
# Triplot RDA spe.hel ~ env2 - scaling 2- wa scores: The explanation for scaling 2 is similar to scaling 1 besides the fact that the angles in the biplot between sites and species, and between sites themselves or species themselves, reflect their correlations. Sites are more spread out with species being centred around origin once again. In the top right, rows (sites) 22, 23, 25 and 27 are clearly visible showing that these sites are positively correlated. Species b_square and s_rect are found in the top_right showing that they are positively correlated but they have really short arrows (as with most of the species on this plot) meaning they are less important for this plot. In the bottom rite, sites 5, 6, 14 and 18 are all visible but more spread out. Regardless of that, they're positively correlated with one another and with the variable cardboard. In the bottom left, site 7, 11, 26, 27 and 28 are visible with 26 being an outlier. These sites are positively correlated with carpet but further away than the sites in bottom right were to cardboard. In top left, sites 3, 4, 10 and 15 are positively correlated but with no species. Site 2 is very close to the origin, with sites 20 and 24 located on the negative y axis and site 13 located on the positive y axis.
# Triplot RDA spe.hel ~ env2 - lc scores - triplots show that the sites are clustered around the origin and based on the lc scores they're all closely related to each other. In the first lc triplot it shows clearly that site 17 plays an important role in locations of b_square and s_rect on the plot. It also shows that row 18 does the same for species gums, s_half, and pink half in the bottom right section of the plot where cardboard plays a big role. The rest of the species are influenced by the rest of the sites in the bottom left of the plot where carpet plays a big role. In scaling 2, however, sites 17 and 18 have less effect on the species as their arrow length decreases significantly.

# Question 7
# The landscape was 200cmx360cm comprised by carpet, cardboard and wood landscape. The 'species' 'immigrated' to the landscape from different angles with them moving towards the centre of the landscape. Since all the different species seemed to show an attraction to the centre, the species richness and relative abundance is very high towards the centre of the landscape. The landscape was sampled by using 50cmx50cm transects with 7 sections sampled across the y-coordinate and 4 sections across the x-axis and each species found in a transect was tallied and entered into an excel spreadsheet. The various analyses such as the species richness, Shannon-Weaver and Simpson's indices, all showed that the centre of the landscape contains the highest species richness, highest species richness and relative abundance index and highest evenness in the landscape. An association matrix for the species data was conducted and showed that the 'gums' species was the most dissimilar to the other species. A PCA analysis was conducted and that revealed that PC1 to PC6 should be taken into consideration as it shows more than 75% of the total variation. The biplots from the PCA analysis shows that the species are located more to the right and side of the biplot, showing that the wind vector (fanning) was the main force acting on the species and is the reason majority of the species immigrated to the centre of the landscape.