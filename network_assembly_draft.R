
## IDEAS
# Manipulate functional trait space by putting an upper bound on potential trait space between 0 and 1
# Set "max.trait.space" to 1. No reason to have it as a different number. It'll make things easier to interpret.
# Use K-means clustering to manipulate the number of plant species. And create species from the dissimilarity data.
# Use euclidean distance to characterize functional trait differences on the random uniform niche axis.
# Run one simulation whether the number of plant and insect species is the same, but I manipulate the potential for functional dissimilarity in the plants.

#### OLD IDEA TO SIMULATE PLANT INDIVIDUALS AND THEN GROUP THEM INTO SPECIES. KEEP ON THE BACKBURNER ----
## specify number of plant individuals
#n.plants <- 100

## use random uniform distribution to simulate plant individuals in functional trait space.
#plant.trait.distrib <- list()
#for(i in 1:length(max.trait.space)){
#  plant.trait.distrib[[i]] <- data.frame(max.trait.space = max.trait.space[i],
#                                         plant.trait = sort(runif(n.plant.species, min = 0, max = max.trait.space[i]))
#  )
#}
#plant.trait.distrib.df <- plyr::ldply(plant.trait.distrib)
#functional.dist <- lapply(plant.trait.distrib, dist)

## use kmeans clustering to identify species in functional trait space
#n.plant.species <- 20 
#functional.species <- lapply(functional.dist, function(x) kmeans(x, centers = n.plant.species))
#str(functional.species)
#functional.species[[1]]$cluster

#### SETUP SIMULATION ----

## load libraries
library(tidyverse)
library(bipartite) # also loads vegan which is required 
library(visreg)
library(cowplot) # loads ggplot

## specify number of plant species
n.plant.species <- 20

## specify range of plant trait space to conduct simulations for
max.trait.space <- seq(0.1, 1, 0.05)

## specify the number of insect species
n.insect.species <- 20

## specify range of insect community structure (i.e. proportion of specialists)
prop.specialists <- seq(0, 1, 0.05)

## arbitrarily specify the range of plant trait space that specialist and generalist insects can interact with
specialist.degree <- 0.05
generalist.degree <- 0.5

## create all possible combinations of plant trait space and insect community structure
plant.insect.combos <- data.frame(
  expand.grid(plants.trait.space = max.trait.space, insect.community = prop.specialists)
  )

#### CREATE PLANT AND INSECT DISTRIBUTIONS IN FUNCTIONAL TRAIT SPACE ----
plant.insect.distrib <- list()
for(i in 1:dim(plant.insect.combos)[1]){
  tmp.trait.space <- plant.insect.combos$plants.trait.space[i]
  tmp.specialists <- round(n.insect.species*plant.insect.combos$insect.community[i],0)
  tmp.generalists <- n.insect.species - tmp.specialists
  
  plant.insect.distrib[[i]] <- data.frame(
    plant.species = factor(1:n.plant.species),
    trait.space = rep(tmp.trait.space, n.plant.species),
    plant.trait = sort(runif(n.plant.species, min = 0, max = tmp.trait.space)),
    insect.species = LETTERS[1:n.insect.species],
    insect.trait = c(sort(runif(tmp.specialists, min = 0, max = tmp.trait.space)),
                     sort(runif(tmp.generalists, min = 0, max = tmp.trait.space))),
    insect.niche.breadth = c(rep(specialist.degree, times = tmp.specialists), 
                             rep(generalist.degree, times = tmp.generalists))) %>%
    mutate(insect.max.trait = insect.trait + insect.niche.breadth) # not constrained within max.trait.space
}

#### DETERMINE WHETHER AN INTERACTION CAN OCCUR OR NOT ----
web.list <- list()
for(i in 1:length(plant.insect.distrib)){
  tmp.web <- plant.insect.distrib[[i]]
  
  plant.insect.links <- list()
  for(j in 1:dim(tmp.web)[1]){
    tmp.plant.species <- tmp.web[j, "plant.species"]
    tmp.plant.trait <- tmp.web[j, "plant.trait"]
    
    plant.insect.links[[j]] <- data.frame(
        plant.species = tmp.plant.species,
        plant.trait = tmp.plant.trait,
        insect.species = tmp.web$insect.species, 
        insect.trait = tmp.web$insect.trait,
        insect.max.trait = tmp.web$insect.max.trait) %>%
        mutate(link = ifelse(plant.trait > insect.trait & plant.trait < insect.max.trait, 1, 0))
  }
  web.list[[i]] <- plyr::ldply(plant.insect.links)
}

# function to create bipartite matrix
special.spread <- function(x) spread(select(x, -plant.trait, -insect.trait, -insect.max.trait), key = insect.species, value = link) %>% select(-plant.species)

# generate interaction matrices
make.webs <- lapply(web.list, FUN = special.spread)

#### CALCULATE NETWORK METRICS FOR EACH WEB ----
web.complexity <- lapply(make.webs, FUN = function(x) networklevel(x, index = c("connectance")))

web.NODF <- lapply(make.webs, FUN = function(x) nestednodf(x)$statistic[3])

web.modularity <- lapply(make.webs, FUN = function(x) computeModules(x)@likelihood)

## merge data together
model.predictions <- data.frame(plant.insect.combos,
                                Connectance = plyr::ldply(web.complexity)$connectance,
                                NODF = plyr::ldply(web.NODF)$NODF,
                                Modularity = plyr::ldply(web.modularity)$V1)

#model.ggplot <- select(model.predictions, -Connectance) %>% 
#  gather(key = Network_Metric, value = Network_Value, NODF, Modularity)

#plot(Modularity ~ NODF, model.predictions)
#plot(Modularity ~ Connectance, model.predictions)
#plot(NODF ~ Connectance, model.predictions)

## still a trade-off between nestedness and modularity after controlling for connectance
#network.lm <- lm(Modularity ~ log(Connectance) + log(NODF), model.predictions)
#visreg(network.lm)
#summary(network.lm)
#plot(network.lm)

## NEED TO THINK ABOUT SOMETYPE OF MODEL LIKE THIS IS WORTH DOING.
#predict.network <- lm(NODF ~ log(Connectance) + plants.trait.space*insect.community, model.predictions)
#visreg(predict.network)
#summary(predict.network)
#plot(predict.network)


NODF.plot <- ggplot(model.predictions, aes(x = plants.trait.space, y = insect.community, fill = NODF)) +
  geom_raster(interpolate = TRUE) + scale_fill_gradient2(midpoint = median(model.predictions$NODF), low = "#56B4E9", high = "#D55E00") +
  scale_x_continuous(breaks = c(0.1,0.25,0.5,0.75,1), labels = c(0.1, 0.25, 0.5, 0.75, 1.0)) +
  xlab("Plant Functional Diversity") +
  ylab("Proportion of Specialist Animals") 

NODF.Mod.cor <- ggplot(model.predictions, aes(x = NODF, y = Modularity)) + geom_point()

#ggplot(model.predictions, aes(x = plants.trait.space, y = insect.community, fill = Connectance)) +
#  geom_tile() + scale_fill_gradient2(low = "white", high = "red")

Mod.plot <- ggplot(model.predictions, aes(x = plants.trait.space, y = insect.community, fill = Modularity)) +
  geom_raster(interpolate = TRUE) + scale_fill_gradient2(midpoint = median(model.predictions$Modularity), low = "#56B4E9", high = "#D55E00") +
  scale_x_continuous(breaks = c(0.1,0.25,0.5,0.75,1), labels = c(0.1, 0.25, 0.5, 0.75, 1.0)) +
  xlab("Plant Functional Diversity") +
  ylab("")
  #ylab("Proportion of Specialist Animals")

fig_2 <- plot_grid(NODF.plot, Mod.plot, nrow = 1, labels = "AUTO", align = "hv")
#first <- plot_grid(NODF.plot, Mod.plot, nrow = 2, labels = c("A","C"), align = "hv")
#fig_1 <- plot_grid(first, NODF.Mod.cor, labels = c("","B"))
save_plot("fig_2.png", fig_2, base_height = 6, base_width = 8.5)

## Null Model Analysis ----

# nestedness
nested.test <- lapply(make.webs, FUN = function(x) oecosimu(x, nestednodf, method = "quasiswap", nsimul = 100))
nested.zscores <- plyr::ldply(lapply(nested.test, FUN = function(x) x$oecosimu$z["NODF"]))$NODF
nested.pvalues <- plyr::ldply(lapply(nested.test, FUN = function(x) x$oecosimu$pval[3]))$V1

# modularity
QuaBiMod <- function(x) computeModules(x)@likelihood
modularity.test <- lapply(make.webs, FUN = function(x) oecosimu(x, QuaBiMod, method = "quasiswap", nsimul = 100)) 

modularity.zscores <- plyr::ldply(lapply(modularity.test, FUN = function(x) x$oecosimu$z))$V1
modularity.pvalues <- plyr::ldply(lapply(modularity.test, FUN = function(x) x$oecosimu$pval))$V1

# merge z-scores and p-values into model predictions
model.results <- cbind(model.predictions, 
                       NODF.zscores = nested.zscores, NODF.pvalues = nested.pvalues,
                       mod.zscores = modularity.zscores, mod.pvalues = modularity.pvalues)

## INTERESTING, IT APPEARS THAT MODULARITY IS QUITE PROBABLE, BUT SIGNIFICANT NESTEDNESS DOES NOT HAPPEN IN THIS HEURISTIC MODEL...I WONDER WHY??? ESPECIALLY, SINCE IT CAN GENERATE HIGHLY NESTED NETWORKS.
## PERHAPS THIS IS AN ARTIFACT OF ALL OF THE GENERALISTS BEING THE SAME, WITH NO GRADIENT IN HOW GENERALIST SPECIES ARE, WHICH WOULD CREATE THE NESTED PATTERN...

## checking validity of the null models
#test <- simulate(vegan::nullmodel(make.webs[[1]], method = "quasiswap"), nsim = 1)

#library(metacom)
#NullMaker(make.webs[[1]], sims = 10, method = "tswap")

#nullmodel(make.webs[[1]], N = 2, method = "mgen", rep.cell = FALSE)

#nestednodf(mgen(as.matrix(make.webs[[1]]), rep.cell = FALSE))
#nestednodf(as.matrix(make.webs[[1]]))
