#Dhiren Vanmari 
#3308987
#BCB743 Quantitative Ecology practical

#Load libraries
library(readr)
library(vegan)
library(ggplot2)

#load data
env <- read_csv("C:/Users/PC/Desktop/Quantitative ecology/quanti mock test/practical assignment/env.csv")
spa <- read_csv("C:/Users/PC/Desktop/Quantitative ecology/quanti mock test/practical assignment/spa.csv")
spe <- read_csv("C:/Users/PC/Desktop/Quantitative ecology/quanti mock test/practical assignment/spe.csv")
spe_mass <- read_csv("C:/Users/PC/Desktop/Quantitative ecology/quanti mock test/practical assignment/spe_mass.csv")


# Question 1 --------------------------------------------------------------
#dimenisions
dim(env) 
dim(spe) 
dim(spa)

# Question 2 --------------------------------------------------------------
#i.	Provide a plot on Cartesian (spatial) coordinates of 

#a)The species richness.
spe.rich <- apply(spe > 0, 1, sum) # sum the number of species > 0 in each row (site) using apply(), MARGIN, is set to 1.

# Bubble map of species richness
plot1 <- ggplot(data = spa, aes(x = x, y = y)) +
  geom_rect(mapping=aes(xmin=125, xmax=150, ymin=140, ymax=220, fill = "Brown"), color="black", alpha=1) + # this reflects the mountain region (box and wood)
  geom_point(aes(x=175, y=300), size=10, shape=1, color="gold4") + #this reflects the bowl in the landscape
  geom_point(aes(cex=spe.rich), fill = "red", colour = "Blue") +
  labs(x = "x coordinate (cm)/ Latitude ", y = "y coordinate (cm)/ Longitude", title = "Map of Species Richness") +
  theme_bw() 
  
plot1
#From the map of  species richness we infer along latitudes y=120 and y = 240 species richness is greater than 2
#site (x=50; y= 240) has the greatest species richness among all sites (spe richness > 10) and  may be as a result of dispersal from  the far North East and a south westerly wind
#along latitudes y= 0 and y= 60, species richness is low (species ricness < 2 ), with an exception a to site (x= 150; y= 0), and may be as a result of low wind dispersal
#quadrats (x=50 ; y =180) and (x= 150 ; y = 240) have a richness of >7  but <10. 
# the absence or low species divesity witnessed at sites (x=150 y=60), (x=200; y= 300), (x= 200 y= 360) may be as a result of low dispersal of species by wind.
#(x =150; y= 180) around the mountain (Brown rectangle) species richness is 2.5 


#b) Shannon-weaver index
spe.shnon <- round(diversity(spe, index = "shannon"),2)

#Bubble map of Shannon-Weaver.
plot2 <- ggplot(data = spa, aes(x = x, y = y)) +
  geom_rect(mapping=aes(xmin=125, xmax=150, ymin=140, ymax=220, fill = "Brown"), color="black", alpha=1) +
  geom_point(aes(x=175, y=300), size=10, shape=1, color="gold4")+
  geom_point(aes(cex=spe.shnon), fill = "red", colour = "Red") +
  labs(x = "x coordinate (cm)/ Latitude ", y = "y coordinate (cm)/ Longitude", title = "Shannon-weaver index of diversity") +
  theme_bw() 

plot2
# the shanon weaver index of diversity  plot indicates :
    #9 sites with an uncertainty index  of 0,  mostly along latitudes y = 0 and y = 60 aswell as along longitude x = 200
    #5 sites with an uncercertainty index of > 0.5 and < 1 
    #6 sites with an uncertainty of > 1 and < 1.5 mainly along longitude x=150 
    #5 sites with an uncertainty > 1.5 and < 2 including the mountainous region
    #3 sites with an uncertainty > 2 mostly along longitude x = 50
#In communities that contain only one species, the uncertainty equates zero,
#In communities with a high species richness, the level of uncertainty increases 

#c) Simpsons index
spe.simp <- round(diversity(spe, index = "simpson"),2)

#Bubble map of Simpson’s Index.
plot3 <- ggplot(data = spa, aes(x = x, y = y)) +
  geom_rect(mapping=aes(xmin=125, xmax=150, ymin=140, ymax=220, fill = "Brown"), color="black", alpha=1) +
  geom_point(aes(x=175, y=300), size=10, shape=1, color="gold4")+
  geom_point(aes(cex=spe.simp), fill = "red", colour = "Dark Green") +
  labs(x = "x coordinate (cm)/ Latitude ", y = "y coordinate (cm)/ Longitude", title = "Simpson index of diversity") +
  theme_bw() 

plot3
#the simpsons index of diversity  represents the probability that two randomly selected individuals will be of the same species:
    #6 sites have 0.00  probability of containing the sames  species. mainly at latitudes of y= 0 and y = 60
    #5 sites have a probability  > 25% and < 50% of  containing the same species
    #6 sites have a probabilty >50% and <75% containing the same species
    #9 sites have a probability >75% and <100% containing the same species mainly at y =180
    #3 sites have 100% probability of containing the same species
#the probability  of two randomly chosen species being the same decreases with increasing species richness, 




# Question 3 --------------------------------------------------------------
#i.	Calculate an association matrix for the species data. (1)
spe.trans <- t(spe)
spe.trans.bray <- vegdist(spe.trans, method = "bray")

spe.trans.bray

#ii.	Produce a plot that is suitable for showing the kinds of information contained in the association matrix. (3)

#load the required libraries and source required functions
library(gclus)
source("coldiss.R")
ass.matrix <- coldiss(spe.trans.bray, byrank = FALSE, nc = 2, diag = "True")

ass.matrix

#iii.	Briefly describe a few of the more ‘obvious’ positive and negative species associations. 
#If possible, add an ‘ecological’ explanation for why this is so. (4)


#magenta = 0 dissimilarity = positive association
#cyan = 1 dissimilarity = negativ association

#From the dissimilarity matrix a positive association exists with dissimilarity close to 0 between s_triangle and b_triangle, s_square, b_square, s_rect, b_rect and s_half. 
#This positive associations between species may be as a result of them both selecting the same
#habitat or requiring the same environmental conditions, some species may have symbiotic relationships or even have predatory, pray relations.

#A negative association is observed  with a dissimalarity close to 1 between s_ triangle and S_c_triangle, b_c_triangle, b_c_square, b_c_rect, b_half, S_c_half,
#this negative association betweeen these species may be as a result of species having differing ecological requirements with some required to 
#live on carpet and some live on wood, some may be easily dispersed by wind whilst others are have a greater resitance. 


# Question 4 --------------------------------------------------------------
#i.	Undertake a PCA on the species data. 
#Comprehensively explain the various (all) components of the summary() of the PCA object. (20)
spe.pca <- rda(spe)
summary(spe.pca, scale= T)

summary(spe.pca, scaling = 1)

#Total inertia
#the total variation of the whole dataset is 5.42
#this is calculated as the sum of all PC axis or sum of all diagonals

#Eigenvalues
# The eigenvalues indicate the importance of the PC axes
# The first PC(PC1) axis having the highest eigenvalue of 1.35  is of  the most importance, followed by the second axis having an eigenvalue of 0.89, 
# The eigenvalues from PC1 until PC17 diminish  from most important to least importance 

#Proportion explained
#indicates the percentage variation explained by each PC axis:
#PC1 explains 24% of the variation
#PC2 explains 16% of the variation 
#PC3 explains 14% of the variation 

#Eigenvectors
#The eigenvectors represent the 'loadings' of the original, untransformed environmental variables on the corresponding eigenvalue (PC axis); 
#they indicate the importance of the original variable along the new reduced dimension captured by the PC axes

  #species scores
      loadings.spe <- scores (spe.pca, display = 'species', scaling = 0)
      sort (abs (loadings.spe[,1]), decreasing = TRUE)
      sort (abs (loadings.spe[,2]), decreasing = TRUE)

  #PC1 species scores indicatee that b_c_squares, b_c_triangle and pink_el are important variables in the reduced PC1 axis
  #PC2 species scores indicate that s_rect, pink_el, b_half and b_c_ triangle are important variables in the  2nd (PC2) reduced axis


#ii.Provide plots of the PCA.

#Source additional functions
source("cleanplot.pca.R")

#Plot      
cleanplot.pca(spe.pca)
 cleanplot.pca(spe.pca, scaling = 1)


#How many axes must be retained, and why? 
 
#source the evplot function 
source ("http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=0")
ev <- spe.pca$CA$eig
evplot(ev)
#The evplot indicates that the first 6 axis must be retained 
#Because they collectively explain 78% of the variation in the dataset

#What patterns become visible from the biplot? Explain your findings in detail. (10)

#PCA-scaling 1
# A scaling 1 PCA biplot reveals the underlying gradients structuring the community; 
# The sites are ordered along the axes according to their positions along these gradients.
# sites 6, 17,18 and 25 have the highest values for s_rect, s_half, b_half and the lowest values for b_c_triangle, b_c_rect, b_rect, s_c_ rect and beans
#which is dissimilar to site 4, 11, 19, and 12  which have the highest values for b_c_triangle, b_c_rect, b_rect, s_c_ rect, beans and the lowest for values for s_rect, s_half, b_half
#site 5 and 6 has the high values for pink_el, s_half, s_triangle, b_half, s_c_half, b_square and gums, which are dissimilar to other sites (sites 2,3,7,8,9,20,26,27 and 28)
# Site 5, 4, 12, 11, 6, 18, 19 have high values for most species as opposed to all other sites.
# The circle of equilibrium contribution allows the identification of the species contributing most to the plotted pair of axes.
# 6 species contribute strongly to axes 1 and 2

#PCA-scaling 2
#A scaling 2 biplot reveals the relationships among species in a correlation-like fashion; 
#the biplot indicates that b_c_triangle, b_c_square, s_c_triangle, b_rect, beans, s_square and b_triangle are (negatively correlated with s_halg and positively corellated with 
#with gums, b_square, s_c_half,s_triangle, b_half, s_rect,and pink_el all of which are negatively correlated with b_c_rect.
#this positive and negative correlations indicates the association of species with one another relative to their sites.





# Question 5 --------------------------------------------------------------
#i.	Undertake an NMDS on the species data, and provide plots of the NMDS. (2)

#remove empty rows
spe1 <- spe[-16, ]
spe.rm <-spe1 [-26: -27,]

#NMDS analysis
spe.nmds <- metaMDS(spe.rm, distance = "bray")


#NMDS plot
plot(spe.nmds, type = "t", main = paste("NMDS/Percentage difference", round(spe.nmds$stress,3)))
abline(v = 0, h = 0, lty =3)

#NMDS with goodness of fit
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main ="Goodness of fit")
points(spe.nmds, display = "sites", cex = gof*100)
abline(v = 0, h = 0, lty =3)

#shepard plot and goodness of fit
stressplot(spe.nmds, main = "Shepard plot")


#ii.What patterns become visible from the ordination plot? Explain your findings in detail. (10)

#The ordination plot indicates:
    #sites 3,4,5,11,12,18 and 20, have higher values for b_triangle(concentrated at site 4), b_c_half and s_c_half(at site 5), b_c_triangle and b_c_square (at site 12), s_square and s_c triangle at sites 18 and 20, and b_c_rect at site 7
    #sites 6,14,16,17 and 19 have a higher abundances of pink_el, b_half (concentrated more at sites 16 and 17) , gums, bsquare (concentrated more at site 6)
    #sites 15, 24, 23 and 25 have higher values for b_rect (concentated more at site 24 and 15), s_triangle (concentrared more at sites 15 and 23),s_half and beans 
    #at sites 10 and 11 beans and s_rect have high values 

#the shepard stress plot indicates a linear R^2 of 0.904 and a non linear-metric fit R^2 value of 0.982.

# Question 6 --------------------------------------------------------------
#i.	Do a constrained analysis on the data, and comprehensively explain all the findings as per the of the summary() of the NMDS object. (20)

#redundancy analysis
spe.hell <- decostand(spe,method = "hellinger")

spe.rda <- rda (spe.hell ~ Carpet + Wood, data = env)
sum.rda <- summary(spe.rda)

sum.rda

#Partitioning of variance 
    #the total variance of all constrained axis and unconstrained axis = 0.706 explaining 100% of the total variation
    sum(eigenvals(spe.rda))
    # the constrained inertia equates to 0.051 and explains 7.2% of the variation
    #the unconstrained inertia equates to 0.66 and explains 92% of the variation

#proportion constrained
    spe.rda$CCA$eig
    #the two variables explain 7.3% of variance in the data
    #the 1st constrained axis RDA1 explains 5.8% of variance  with an eigenvalue of 0.040
    #the 2nd constrained axis RDA2 explains 1.5% of variance with an eigenvalue of 0.010
    #the eigenvalues of the two constrained axis decreases as does PC axis from PC1 to PC16
    #PC1 > RDA1 indicates that the first residual structure (axis) of the data has more variance than some of the structures that can be explained by the explanatory variables

#proportion unconstrained
    #the sum of all unconstrained  PC axis explain 92.7% of variance in the data set.
    spe.rda$CA$eig
    #PC1 explains 21% of the variation with eigenvalues 0.14
    #PC2 explains 16% of the variation with eigenvalues 0.11
    #PC3 explains 11% of the variation with eigenvalues 0.07

#Accumulated constrained eigenvalues: 
    #the cumulative amounts of variance expressed as proportions of the total explained variance, as opposed to their
    #contribution to the total variance described above

#Species scores
    #are the coordinates of the tips of the vectors representing the
    #response variables in the bi- or triplots. As in PCA, they depend on the scaling chosen.

#Eigenvalues for constrained axes
    #RDA1 is of more importance than RDA2 as RDa explains 80% of the accumulated constrianed eigenvalues, whilst RDA2 explains 20%



#ii.	Run the necessary permutation tests. Explain the outcomes. (6)
anova (spe.rda)

#The function of permutation test is to evaluate whether the variation explained by the model (variables)
#is higher than would be variation explained by the same number of randomly generated variables.
# the model for both constrained axis is insignificant as p-value is < 0.05  for 999 permutations  


#iii.	Produce the necessary biplots. Explain the findings. (14)
plot(spe.rda,
     scaling = 1,
     display = c("sp", "cn"),
     main = "Biplot RDA spe.hel scaling 1"
)
spe.sc1 <- 
  scores(spe.rda, 
         choices = 1:2, 
         scaling = 1, 
         display = "sp"
  )
arrows(0, 0, 
       spe.sc1[, 1] * 0.92,
       spe.sc1[, 2] * 0.92,
       length = 0, 
       lty = 1, 
       col = "red"
)

# Scaling 2
plot(spe.rda, 
     display = c("sp","cn"), 
     main = "Biplot RDA spe.hel scaling 2 "
)
spe.sc2 <- 
  scores(spe.rda, 
         choices = 1:2, 
         display = "sp"
  )
arrows(0, 0, 
       spe.sc2[, 1] * 0.92, 
       spe.sc2[, 2] * 0.92,
       length = 0,
       lty = 1,
       col = "red")



#the biplots indicate that carpet and  wood  do not play an important role in the dispersions of species along the 1st axis
#the scaling 2 biplot shows 6 species namely b_square, s_rect, b_half, pink_el, gums and s_half are associated with wood surface , whilst all other species are correlated to carpet surface.   




# Question 7 --------------------------------------------------------------
#Write down your understanding of the ‘ecology’ of the landscape that was sampled and analysed in the various steps, above.

#sampling 
# systematic sampling was undertaken accross an area of 275cm, x 325cm, with transects laid at 50 cm intervals
#along each of the 5 transects, quadrats (30cm x 30cm) was laid at 60cm intervals. this was done to obtain a consistent system of sampling and to reduce biasness
#sampling effort was considered to be at a maximised as most of the site was sampled

# the species distribution accross the landscape is somewhat, influenced by the dispersion of species from the far north east  of the map (x=50 y= 360)  and a southerly  blowing wind
plot1 #we therefore find that species richness at latitudes is low in the southern regions of the map (latitudes x=0 and x= 60), this may be as a result of diminishing wind strength from its source.
#it is also observed that species richness is considerably high between latitudes 120 and 300. this may be as a result of species weight realtive to the wind and thus contribute to their place of settlement.
sort(colMeans(spe_mass),decreasing = F) # which might mean that lighter species will disperse further than heavier species from poit of dispersion.
# at site (x=50 and y = 240) species richness is greater than 10, and may be as a result of  similar habitat selection.
# on the mountain (wood and box) region we see a species richness of 2.5 these species may be better adapted to either living at a altitude higher than other species or prefer mountain substrate (wood and box). 

plot2 #of the shannon indicates that sites with higher species richness have higher shannon indexes and sites with low species diversity have lower indexes, this may indicate that similar species may occupy similar regions as a result of habitat preference or requirements.


ass.matrix # the dissimilarity observed between species may ecologically described by the habitat specifications of species ranging from wind, weight and site of dispersion
#positive associations between species may be as a result of them requiring similar habitats  with similar environmental conditions, whilst some species may have symbiotic relationships or even have predatory, pray relations.
#the negative association betweeen species may be as a result of species having differing ecological requirements with some required to live on carpet and some live on wood, some may be easily dispersed by wind whilst others are have a greater resitance 
#the dissimilarity between species is what allows them to coexist otherwise competiton would cause them to outcompete one another.

#from the nmds biplots we infer 
#that sites between latitudes 120 and 240 are greatly structured by b_c_triangle, b_c_square, b_c_half, pink_el,b_half and s-rect.
#sites 3,4,5,11,12,18 and 20, have higher values for b_triangle(concentrated at site 4), b_c_half and s_c_half(at site 5), b_c_triangle and b_c_square (at site 12), s_square and s_c triangle at sites 18 and 20, and b_c_rect at site 7
#sites 6,14,16,17 and 19 have a higher abundances of pink_el, b_half (concentrated more at sites 16 and 17), gums and b_square (concentrated more at site 6)
#sites 15, 24, 23 and 25 have higher values for b_rect (concentated more at site 24 and 15), s_triangle (concentrared more at sites 15 and 23),s_half and beans 
#at sites 10 and 11 beans and s_rect have high values 

#from the rda biplots we infer 
#that carpet and wood do not play an important role in the dispersion of species hence the reason we see insignificant p-value.
#although most species were dispersed accross the carpet. 



###################################### End of assingment #############################################
