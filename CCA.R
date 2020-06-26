# http://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html

install.packages("vegan")
library(vegan)

# Download and extract the following zip file. 
# It has three matrices: spe (with species abundances per site), 
# env (with environmental variables per site), and spatial (with longitude and latitude values per site; 
# this matrix will be only used in partial CCA).

# http://www.mediafire.com/file/r54e4ymw5wz1t0a/cca.zip

# Import data into R:

spe <- read.csv(choose.files(), row.names=1, sep=";")
env <- read.csv(choose.files(), row.names=1, sep=";")
spatial <- read.csv(choose.files(), row.names=1, sep=";")

# Apply log+1 transformation to your species occurrences data (spe matrix) 
# in order to correct for possible statistical errors associated to rare or very common species:

spelog <- decostand(spe, "log")

# Perform CCA. We need to specify that spe (species distribution matrix) 
# is explained by env (environmental matrix).

ccamodel <- cca(spe ~ ., env)

## To perform a Partial CCA, we need to use a third matrix (conditioning matrix). 
# To do that, we must first combine variables from the environmental ("env") 
# and conditional ("spatial") matrices. We have to do that in order to later apply 
# another function ("ordistep" function). After we have combined all variables, 
# we will apply the Partial CCA using a formula where we specify the 
# response ("spe" matrix), the constraint variables (each of the variables from "env" matrix), 
# and the conditioning variables (variables from the conditioning matrix; in our case "spatial" matrix)

envspatial <- cbind(env, spatial)

nams <- names(envspatial)

partialccamodel <- formula(paste("spe ~", paste(nams[1: (length(envspatial)-(length(spatial)) )], 
             collapse = " + "),
             "+ Condition(", paste(nams[(length(envspatial)-(length(spatial)-1) ):length(envspatial)], collapse ="+"),")"))

partialccamodel <- cca(partialccamodel, envspatial)

# However, we also have to automatically select variables of "env" matrix 
# that best explain "spe" matrix. We can do that by using a stepwise model
# from "ordistep" function. Let us do that with our "ccamodel" (not the partial cca).

finalmodel <- ordistep(ccamodel, scope = formula(ccamodel))

# Then, we can calculate Variance Inflation Factors (VIF) for each 
# of the constraints (variables) from the "env" matrix (environmental matrix). 
# If we find an environmental variable with VIF>10, 
# we'll know that this variable presents colinearity with another or other variables. 
# In that case, we would have to delete the variable from our initial dataset 
# and redo all the analysis. In our example, no variable is redundant with each other 
# (all of them have VIF<10).

vif.cca(finalmodel)

# Let us now fit the stepwise model (ordistep function) using our 
# partial cca model ("partialccamodel"). Note that we will use X (longitude) 
# and Y (latitude) variables from the previously created "envspatial" 
# object as our conditioning variables. After "vif.cca", note that "X" (spatial) 
# variable has a value > 10, so one should consider to delete the variable before the analysis.

partialccamodel <- formula(paste("spe ~", paste(nams[1: (length(envspatial)-(length(spatial)) )], collapse = " + "),"+ Condition(", paste(nams[(length(envspatial)-(length(spatial)-1) ):length(envspatial)], collapse ="+"),")"))

simplemodel<-cca(partialccamodel, envspatial)

finalmodelpartial<- ordistep(simplemodel, scope=formula(partialccamodel))

vif.cca(finalmodelpartial)

### Ok, let us now call "ccamodel" object (not the "partialccamodel") to see how to interpret results.

finalmodelpartial

## Note that "Total Inertia" is the total variance in species (observations matrix) distributions. 
# "Constrained Inertia" is the variance explained by the environmental variables (gradients matrix).
# The "Proportion" values represent the percentages of variance of species distributions 
# explained by Constrained (environmental) and Unconstrained variables. Eigenvalues of 
# constrained and unconstrained axes represent the amount of variance explained 
# by each CCA axis (graphs usually present the first two constrained axes, so take a look at their values).


# This is a critical step when doing the CCA. When we have our final model, 
# we must use permutation tests to observe if our whole CCA model, 
# the CCA terms (environmental varibles), and CCA axes explain more 
# variance of "spe" (observations) matrix than expected by chance 
# (tests should be significant; p values lower or equal 0.05). 
# If the tests are not significant, there is no point in using the CCA. 
# In our example, we'll see that we can continue using the CCA results.

# Testing the significance of the CCA model:
anova.cca(finalmodel)

# Testing the significance of terms (environmental variables):
anova.cca(finalmodel, by = "terms")

# Testing the significance of CCA axes (at least the first two or three 
# should present a significant p value):

anova.cca(finalmodel, by = "axis")

### Finally, we may want to generate graphs in order to better understand results. 
# To do that, we have to take a look at "xlim" (limits of x axis), 
# "ylim" (limits of y axis), and "display" 
# (if we want to observe species, environmental gradients, 
# and/or sites in the graph) arguments. For example, 
# to show only species scores along CCA axes, use only "sp" inside display" argument.

plot(finalmodel, xlim = c(-1.5, 2), ylim = c(-1, 1.5), display = c("sp"))

# If you want to show species and environmental gradients in our plot, 
# use both "sp" and "cn" inside "display" argument,

plot(finalmodel, xlim = c(-3, 3), ylim = c(-3, 3), display = c("sp", "cn"))

# To show species, environmental gradients, and sites, use "sp", "cn", 
# and "wa" inside "display" argument,

plot(finalmodel, xlim = c(-3, 3), ylim = c(-3, 3), display = c("sp", "cn", "wa"))

### To interpret the graph, we need to see if species or sites are close to 
# environmental gradients. For example, Alopacce, Alopfabr and Arctperi 
# are very related to Bare Sand, while Pardiugu is very related to Fallen Twigs.