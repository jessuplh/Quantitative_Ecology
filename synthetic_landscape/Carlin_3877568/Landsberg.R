# Landsberg.R
# QE_prac.R
# Quantitative Ecology practical exercise 
# 18 August 2018
# Carlin Landsberg

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(vegan)

source("Desktop/NEwR_updated_material_R332-NEwR_ed1/coldiss.R")
source("Desktop/NEwR_updated_material_R332-NEwR_ed1/cleanplot.pca.R")
source("Desktop/NEwR_updated_material_R332-NEwR_ed1/evplot.R")

# Load data ---------------------------------------------------------------

env <- read_csv("Desktop/QE_prac_test/env-1.csv")
spa <- read_csv("Desktop/QE_prac_test/spa.csv")
spe <- read_csv("Desktop/QE_prac_test/spe.csv")
mass <- read_csv("Desktop/QE_prac_test/spe_mass.csv")

# Question 1: Dimensions ----------------------------------------------------------

dim(env)
dim(spa)
dim(spe)
dim(mass)

# ???????????????????????????????????????????????????????????????

# env: 28 rows, 3 columns 
# spa: 28 rows, 2 columns
# spe: 28 rows, 17 columns
# mass: 28 rows, 17 columns

# ????????????????????????????????????????????????????????????

# Q2 ----------------------------------------------------------------------

# i)

# mutate data to create a row with site numbers
mut_spa <- spa %>% 
  mutate(ID = 1:28)

# first create a site map
site_map <- ggplot(data = mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "lightblue") +
  geom_text(aes(label = ID), colour = "red") +
  labs(title = "Study site", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
site_map

# site map including first and second points of dispersal
site_map_full <- ggplot(data = mut_spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "lightblue") +
  geom_text(aes(label = ID), colour = "red") +
  annotate("text", label = "First point of dispersal", x = 0.75, y = 450, colour = "black", size = 3) +
  annotate("text", label = "Second point of dispersal", x = 0.75, y = -50, colour = "black", size = 3) +
  labs(title = "Study site", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
site_map_full

# a) Species richness
sit.sum <- apply(spe > 0, 1, sum)
sit.sum.tib <- as.tibble(sit.sum)
site_df <- data.frame(spa = sit.sum.tib)
sit.sum.df <- data.frame(spe_rich = sit.sum)
spe.rich <- sit.sum.df %>% 
  mutate(N = 1:28)
sit.sum.mut <- sit.sum.tib %>% 
  mutate(N = 1:28)
sit_df <- as.data.frame(sit.sum.mut)

# plot for species richness
sp_rich <- ggplot(data = spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_line(colour = "lightblue") +
  geom_point(aes(cex = spe.rich$spe_rich), shape = 1, colour = "red", show.legend = FALSE) +
  labs(title = "Species richness", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
sp_rich

# step plot for species richness
sp_rich_step <- ggplot(sit_df, aes(x = N, y = value)) +
  geom_step(colour = "lightblue") +
  geom_text(aes(label = N), colour = "red") +
  labs( x = "Positions of sites", y = "Species Richness", title = "Species Richness") +
  theme_bw()
sp_rich_step

# b) Shannon-weaver 
H <- as.tibble(diversity(spe, index = "shannon"))
H

# Shannon plot
shan <- ggplot(spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "lightblue") +
  geom_point(aes(cex=H$value), shape = 1, colour = "red", show.legend = FALSE) +
  labs(title = "Shannon-Weaver Diversity Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
shan

# c) Simpson index
simp <- as.tibble(diversity(spe, index = "simpson"))
simp

# Simpson plot
simpson <- ggplot(spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_path(colour = "lightblue") +
  geom_point(aes(cex=simp$value), shape = 1, colour = "red", show.legend = FALSE) +
  labs(title = "Simpson Index ", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
simpson

ggarrange(sp_rich, sp_rich_step, shan, simpson)

# ii) 

# ???????????????????????????????????????????????????????????????

# Species richness: Species richness is high between 100cm and 300cm latitude (y). 
#                   at the transect x=50, overall species richness is higher than the other transects
#                   thic could be explained by the fact that both points of entry or disperal
#                   for species were located verly close to this transect (see site_map_full).
# Species richness step plot: This plot shows that as sites reach the middle area of the study site, their species richness increases.
#                             sites closer to the edges of the transect AND further from the points of dispersal have fewer species/ individuals present. 
#                             site 5 had the highest species richness as this was in the middle of the two point of dispersal. 
#                             Site 27 and 28 have no species because they were the furthest away from the points of dispersal and thus species could not reach this area to colonise it. 
#                             Because of the wind applied at the two points of dispersal and the light mass of the species (i.e. easily blown by the wind),
#                             species were not very abundant at the points of dispersal as they were quickly blown further into the study site.. 
#                             In addition, species that colonised between 0cm and 100cm latitude (y) at the first disperal event may have later been blown from these areas further "inland",
#                             thus resulting in the high species richness we observe in the middle of the study site.
# Shanon weaver: This index looks at diversity of species present as well as eveness. W observe high species diversity in the 3 middle quadrats of each of the transects. 
#               At the first point of dispersal, species diversity is high but at the second point of dispersal species diversity is low.
#               Reason for this would be that the gums and beans released at this point had a large biomass and did not disperse far and other species were able to disperse much further away.
#               Additionally, wind may have blown lighter species inward leaving only gums and beans at this area.
# Simpsons index: This index takes into account species diversoty as well as abundnace. So, diversity and abundance is relatively high throughout the landscapae,
#                 however, at the edges of the study site (transect 1, 2, and 4, at >100cm (y)), diversoty and abundnace are low as species may not have initially reached the area or were blown away during the second dispersal event.

# ???????????????????????????????????????????????????????????????

# Q3 ----------------------------------------------------------------------

# i)

# transpose the species data
trans_spe <- t(spe)

# bray-curtis (default)
bray_spe <- vegdist(trans_spe)
bray_spe

# OR SIMPLY; no need to transpose or standardize etc. 
cor(spe)

# ii)

# plotting the association matrix using the coldiss function (loaded earlier)
coldiss(bray_spe, byrank = FALSE, diag = TRUE)
  # Pink = more similar (dissim close to 0)
  # Blue = least similar (dissim close to 1)

# iii) 

# ???????????????????????????????????????????????????????????????

# the gum species are poorly associated with all of the other species due to it not being able to disperse as well as a result of its large biomass.
# Majority of crumbled species are associated with each other. Dispersal ability was controlled by the wind and as a result, heavier species such as the gums were poorly correlated but crumbled species were had higher associations.
# Flat species are somewhat associated, however they are very easily moved and strongly influenced by the wind. They thus dispersed so widely and did not settle in close assocoiaions.
# At the second dispersal event, these flat species may have also been more easlily displaced due to its dimensions that are easilt susceptible and swept up by the wind. 

# ???????????????????????????????????????????????????????????????

# Q4 ----------------------------------------------------------------------

# i)

# RDA/ PCA run on species data, scale = 2 standardizes data
spe.pca <- rda(spe, scale = TRUE)
spe.pca

summary(spe.pca)

# ???????????????????????????????????????????????????????????????

# Inertia refers to the sum of the diagonals (1+1+1+1...) of the correlation matrix. Inertia refers to the variation in the data.
# The sum of all the correlations adds up to the total number of variables in the dataset.
# In this case, inertia = 17 (there are 17 species or variables present in the spe dataset)

# The eigenvalues are the measures of importance or variance of the PC axes. The eigenvalues are the tranformed datapoints (e.g. by Euclidean distance)
# The first PC axis will have the higest eigenvalue and therefore the highest explanatory power.
# Eigenvalues for each PC axis will decrease and the sum of all the eigenvalues for all the PC axis should be equal to the total inertia.
# In this case, PC1 has the higest explanatory power, explaining 22.3% proportion of the variance. 

# The eigenvectors show the direction of the line. They represent the importance of the original variables before transformation with regards to the corresponding eigenvalue.
# Eigenvectors and eigenvalues always exist in pairs (i.e. each eigenvector will have a corresponding eigenvalue). 

# Proportion explained refers to the proportion of variance explained by the PC axis in question.
# This value is calculated by dividing the eigenvalue of the PC axis by the total inertia. 

# Cumulative proportion is calculated by adding the proportion explained by the PC axes. The sum of all the PC axes should be equal to 1 (i.e. all PC axes should explain 100% of the variance).
# In this case, PC1 shows the hghest explanatory power as it explains 22.3% of variance and PC2 explains 14.45% of variance. Collectively, they explain 36.86% of the variance observed. 

# Species scores refer to the coordinates of the arrow heads of the  response variables (species) on the PCA plot
# Site scores refer to the coordinates of the sites on the PCA plot

# ???????????????????????????????????????????????????????????????

# ii)

# PCA biplot scaling 1 and scaling 2
par(mfrow = c(1, 2))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") 

# PCA biplots, including cleanplots
par(mfrow = c(1, 1))
biplot(spe.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(spe.pca, main = "PCA - scaling 2 ('species')") 
cleanplot.pca(spe.pca, scaling = 1, mar.percent = 0.08)
cleanplot.pca(spe.pca, scaling = 2, mar.percent = 0.04)

# ???????????????????????????????????????????????????????????????

# 6 axes must be retained as the cumulative proportion explained by the first 6 PC axes is equal to 75.96%
# to calculate the proportion explained by each pc axis
p1 <- spe.pca$CA$eig[1] / sum(spe.pca$CA$eig)
p2 <- spe.pca$CA$eig[2] / sum(spe.pca$CA$eig)
p3 <- spe.pca$CA$eig[3] / sum(spe.pca$CA$eig)
p4 <- spe.pca$CA$eig[4] / sum(spe.pca$CA$eig)
p5 <- spe.pca$CA$eig[5] / sum(spe.pca$CA$eig)
p6 <- spe.pca$CA$eig[6] / sum(spe.pca$CA$eig)
# then sum them together and multiply by 100 to determine the percentage of variance explained.
((p1 + p2 + p3 + p4 + p5 + p6) * 100)

# Another way to determine how many axes should be retained is to plot the eigenvalues
ev <- spe.pca$CA$eig
# To do this, we apply Kaiser-Guttman criterion to select axes
ev[ev > mean(ev)]
evplot(ev)
# therefore, this test tells us that we should retain the first 7 axes as their eigenvalues are higher than that of the mean of all the eigenvalues.
# on the plot, only PCA1 and PCA2 are shown. 

# SCALING 1: scaling one shows the grouping of the sites. Sites that share similar species are plotted more closely together.
#            In this case, majority of sites are grouped very close together on the left hand side of the graph.
#            Sites 11, 19, 5 and 18 are the more obvious outliers. 
#            We can thus deduce that site 5 is mostly associated with small crumbled half circles, site 11 and 19 with big rectangles and site 18 with gums.
#            At these sites, species richness were relatively high, however, the patterns observed seem rather random and the reasons sites and species are associated may not have been taken into account in the PCA test.

# SCALING 2: Scaling two shows the correlation between the species. Species that plot with smaller angles btween them are more strngly correlated with each other.
#            All variables, or in this case, species, are closely correlated with the species adjacent to it as all the angles between them are acute. 
#            Some of the more obvious and stronger correlations occur between the big and small crumbled triangles.
#            Additionally, the big triangles are strongly correlated with small squares and big crumbled squares are correlated with big crumbled half circles.
#            From this, one could infer that crumbled species are asscociate with one another and flat species are associated with one another. 

# No species are explicitly negatively correlated with another as most species are plotted close to each other and none are plotted directly opposite another. 
# Overall, most of the sites and species are randomly distributed (with no clear gradient) with only a few exceptions or otliers that are associated with the species rich areas.
# Roughly half of the species fall outside of the radius of the circle and these can be interpreted as the species that play more of a role in determining the spread of sites.

# ???????????????????????????????????????????????????????????????

# Q5 ----------------------------------------------------------------------

# i)

# remove the rows with NAs
sp.1 <- spe[-16,] 
rm_sp <- sp.1[-26:-27,]

spe.nmds <- metaMDS(rm_sp, distance = "bray")
spe.nmds

summary(spe.nmds)

spe.nmds$stress

dev.new(title="NMDS on species data - Percentage difference")
plot(spe.nmds, type = "t", main = paste("NMDS/Bray - Stress =", round(spe.nmds$stress, 3))) # percentage difference

# stressplot and goodness of fit plot
# dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe.nmds, main = "Shepard plot") 
gof <- goodness(spe.nmds) # goodness of fit 
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof*300)

# ii)

# ????????????????????????????????????????????????????????????

# Stress = 0.136 which indicates a fair fit. 

# The Shepard stress plot shows relationship between the actual dissimilarities between objects (from the original dissimilarity matrix) and the ordination distances (i.e. the distances on the final plot).
# There appears to be a large amount of scatter (i.e. a poor linear relationship), therefore, this ordination is not representative of the original distances.

# A number of objects (sites) are plotted close together meaning that they may be more similar in sharing species. These arent very strong correlations however as sites are still somewhat of a distance away from each other. 
# They form loose clusters or groupings and a few sites (7, 8, 5, 22) are plotted very far away from all the other sites. These sites are thus the most dissmiliar from the other sites present.
# Upon inspection of the study site map, we see that these sites are so dissimilar from the rest because site 5 has the most species/individuals present whereas sites 7, 8, and 22 are at the edges of the landscape which may have had an effect on the species occurring there.

# it is also important to note the size of the bubbles in the goodness of fit plot.
# in this case, majority of the bubbles are very larger which indicate that the sites are poorly fitted. 
# Sites 7 and 8 bubbles are a bit smaller thus indicating that these sites are slightly better fitted.

# ????????????????????????????????????????????????????????????

# Q6 ----------------------------------------------------------------------

# i)

# transform the data using hellinger-transform
spe.hel <- as.tibble(decostand(spe, "hellinger"))
spe.hel

# run rda on hellinger transformed data by ALL environmental variables (.,)
spe.rda <- rda(spe.hel ~ ., env) # or rda(spe.hel ~ Carpet + Cardboard + Wood, data = env))
spe.rda

summary(spe.rda)  
  # scaling 2 default

# Canonical coefficients from the rda object
coef(spe.rda)
  # RDA1 and RDA2

# ????????????????????????????????????????????????????????????

# In an RDA, variance is proportioned into constrained and unconstrained proportions.
# In this case, 92.7% of the variance is explained by the unconstrained fraction of this test and only 7.3% of the varainace can be exlained by the constrained fraction.

# This analysis yielded 2 canonical axes, namely RDA1 and RDA2. 
# RDA1 explains 5.8% of the varaince observed and RDA2 explains 1.5% of the variance observed.
# These values add up to the 7.3% of variance explained by the constrained canonical axes.
# The other 92.7% is made up of the unconstrained axes, i.e. PC1-PC16 where PC1 explains most of the variance (21.1%) followed by PC2 (15.9%)

# It is important to note that the canonical (RDA1, RDA2) eigenvalues measure amounts of variance explained by the RDA model, 
# whereas the residual (PC1-PC16) eigenvalues measure amounts of variance represented by the residual axes, but not explained by the RDA model.
# The constrained axes represent the amount of variance in the species data; i.e. how species spread in response to environmental variables. 
# The unconstrained axes represents the amount of variance that still remains after the consrained axes have been retained. This represents variance that cannot be explained by the environmental variables.
# This can be seen as either noise or something else, such as  another environmental variable, that has not been accounted for.

# A reason for these constrained and unconstrained values could simply be statistical noise however it is likely that the constrained values are so low because in the case of the landscape, wind and point of entry may have been important environmental variables that were not accounted for.
# It was clear from the study that wind was the main influencers of where species initially landed and where they moved to after subsequent gusts of wind caused displacement of species. 

# Accumulated constrained eigenvalues represents the cumulative amounts of variance expressed as proportions of the total EXPLAINED (constrained) variance, as opposed to their contribution to the TOTAL variance or inertia.
# Therefore, RDA1 explains 79.8% of the constrained variances and is thus has the most explanatory power of the two axes. RDA2 explains 20.9% of the constrained variance. 
# In total, RDA1 and RDA2 cumulatively explain 100% of the constrained varaince and this, in turn, explains 7.3% of the TOTAL variance. 

# Species scores refers to the coordinates of the tips of the vectors representing the RESPONSE variables. 

# The site scores or the weighted sums of species scores represent the coordinates of the sites as expressed in the space of the RESPONSE variables. 

# Site constraints or the linear combinations of constraining variables represents the coordinates of the sites in the space of the EXPLANATORY variables. These are the fitted site scores.

# Biplot scores for constraining variables represent the coordinates of the tips of the vectors representing the EXPLANATORY variables. 

# coef(spe.rda) gives us the canonical coefficients from the RDA project. This is the equivalent of regression coefficients for each explanatory value (i.e. carpet, wood, cardboard) on the canonical axes (RDA1, RDA2)

# ????????????????????????????????????????????????????????????

# ii)

# ***** PERMUTATIONS *****

# H0: environmental variables have no effect on species distributions
# H1: environmental variables have an effect on species distributions

# Global test of the RDA result
anova(spe.rda, permutations = how(nperm = 999))
  # entire RDA model

# Tests of all canonical axes
anova(spe.rda, by = "axis", permutations = how(nperm = 999))
  # RDA1 and RDA2

anova(spe.rda, by = "term", permutations = how(nperm = 999))
  # variables carpet and cardboard (wood; NA)

# Variance inflation factors (VIF)
vif.cca(spe.rda)
  # carpet and cardboard = 1.928571, wood = NA

# Apply Kaiser-Guttman criterion to residual axes (which axes should be  slected or retained)
spe.rda$CA$eig[spe.rda$CA$eig > mean(spe.rda$CA$eig)]
# PC1, PC2, PC3, PC4, PC5, PC6, PC7

# R2 adjusted value 
R2adj <- RsquareAdj(spe.rda)$adj.r.squared
R2adj

# ????????????????????????????????????????????????????????????

# anova on global RDA: Pr > 0.05, therefore we fail to reject the null hypothesis.
# The environmental variables have no significant effect on the distribution of species.
# The RDA model as a whole does not significantly explain the varaince observed.

# anova by axis: Pr > 0.05 for both RDA1 sand RDA2. The two axes do not significantly explain the variance observed.

# anova by term: Pr > 0.05 for cardboard and Pr > 0.05  for carpet. These two environmental variables do not significantly explain the variance observed. 

# Overall, the permutation results tell us that the model, axes, and environmental variables do not significantly explain the variance observed,
# it can thus be inferred that, from the results of the RDA, we fail to reject the null hypothesis and thus environmentak variables do not affect species distributions.  

# ????????????????????????????????????????????????????????????

# iii) 

# Triplot; scaling 1 (wa scores)
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

# Triplot scaling 1 (lc scores)
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

# ????????????????????????????????????????????????????????????

# Firstly, it is important to note that this is different to the PCA biplot as it is a triplot because it contains three entities mamely, sites, response variables (spe) and explanatory variables (env).
# wa scored refers to the weighted sums of species 
# lc scores refers to the fitted site scores

# wa scores (scaling 1): distance biplot: The angles between response and explanatory variables in the biplot reflect their correlations
# wa scores (scaling 2): correlation biplot: The angles in the biplot between response and explanatory variables, and between response variables themselves or explanatory variables themselves, reflect their correlations 

# because these two axes cumulatively explain only 7.3% of the total variance, we are not confident that the major trends have been well-modelled in this analysis. 
# Additionally, the adjusted R2 value is very low (-0.001512655), which is very low (although this may be expected because environmental variables are known to have a lot of noise)

# Distances between object points represent their  Euclidean distances. Thus, objects plotting closer together can be expected to have similar variable values.
# The angles between the response variable vectors and explanatory variable vectors reflect their linear correlation. 

# The angles between all vectors reflect their (linear) correlation, a vector pair describing an angle of 90?? are uncorrelated.
# smaller angles between vectors represent strong, positive correlation.
# the right angle drawn from explanatory vector to response variable approximates the environmental avriable's value realised for that object. 

# In scaling 1, the projection of ordinated objects onto a vector are represented, thus approximating the variable's value realised for that object.
# in this case, right angles drawn from the carpet vector show that majority of species are associated or influenced by the carpet.
# technically speaking, these species should occur where there are large amounts of the variabe (i.e. carpet). 
# The carpet covered majority of the landscape and that why majority of species are found there. 
# Pink elephants, small half circle big half circles and gums can be expected to be influenced by cardboard.
# i.e. more cardboard woudl be expected where these species occur compared to other species. 
# Upon interpreting this looking at rows, we can infer that less carpet would be found at rows 17 and 18. 
# This holds true as at row 17 and 18 (or site 17 ad 18), there was no carpet, but rather cardbaord and wood was present. 

# In this case, scaling 2 shows that the angle between cardboard and carpet is very large (> 90??) and thus the two vectors/variables are uncorrelated. 
# Additionally, row 17 and 18 are far displaced from the grouping or clustering of the rest fo the rows.

# ????????????????????????????????????????????????????????????

# Q7 ----------------------------------------------------------------------

# ????????????????????????????????????????????????????????????

# The landscape consisted of dimentions roughly 250cmx360cm. Species entered the landscape via several entry or dispersion points indicated on site_map_full.
# At these entry points, a strong wind was created in order to disperse the species across the landscape which consisted of mainly carpet, a wooden slope and a cardboard box. These were then referred to as the environmental variables. 

# At the first dispersal event at the first entry site, species colonised an area across the landscape. crumbled species were the first to colonise.
# At the next dispersal event at the first entry point again, flat species were then allowed to colonise the area. Because of the dimensions of the species and thus their susceptibility to wind, they may have been blown further by the wind.
# The rest of the species that were not made of paper were dispersed into the landscape from a different (again, a new dispersal event), still via wind.
# At each dispersal event, the species that were initailly colonised were displaced and in some cases, were displaced several times.

# The area was sampled using transects and quadrats at intervals along these transects. (transects at 50cm intervals, quadrats at 60cm intervals)
# Species at each "site" were determined and the number of each species at each site was recorded.
# Hereafter, the relative mass of all species recorded at each site was calculated and recorded.

# As a result, it seems that most individuals settled towards the middlemost quadrats of each transect.
# However, relative masses did not follow the exact same patterns as can be observed looking at the mass_spe_combined plot.

# Wind may have been the most important environmental factor influencing the spread of species on the landscape.
# However, wind was not quantified in this experiment and thus majority of analyses that were performed, above, did not significantly explain the reasons for spread of species.
# As a result, in most cases, we could not confidently accept the results as true fot this experiment. 

# Additionally, trampling may have occurred during the sampling period which may have also influenced spread of species and their displacement. This, along with wind, may have been important environmental factors to include in this study.

# ????????????????????????????????????????????????????????????

# EXTRA -------------------------------------------------------------------

mass.pres <- apply(mass > 0, 1, sum)
mass.pres
mass_tib <- as.tibble(mass.pres)

mass_rich <- ggplot(data = spa, aes(x = as.factor(x), y = as.numeric(y))) +
  geom_line(colour = "lightblue") +
  geom_point(cex = 5*mass.pres/max(mass.pres), shape = 1, colour = "red", size = 5) +
  labs(title = "Mass", x = "x coordinate (cm)", y = "y coordinate (cm)") +
  theme_bw()
mass_rich

mass_spe_combined <- ggarrange(sp_rich, mass_rich)
