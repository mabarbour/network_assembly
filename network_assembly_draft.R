
## IDEAS
# Manipulate functional trait space by putting an upper bound on potential trait space between 0 and 1
# Set "max.trait.space" to 1. No reason to have it as a different number. It'll make things easier to interpret.
# Use K-means clustering to manipulate the number of plant species. And create species from the dissimilarity data.
# Use euclidean distance to characterize functional trait differences on the random uniform niche axis.
# Run one simulation whether the number of plant and insect species is the same, but I manipulate the potential for functional dissimilarity in the plants.

library(tidyverse)

## specify number of plant individuals
#n.plants <- 100

## specify number of plant species
n.plant.species <- 20

## specify range of plant trait space
max.trait.space <- seq(0.1, 1, 0.1)

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

## specify the number of insect species
n.insect.species <- 20

## proportion of specialists
prop.specialists <- seq(0, 1, 0.1)

## determine the range of trait space that specialist and generalist insects can interact with
specialist.degree <- 0.05
generalist.degree <- 0.5


plant.insect.combos <- data.frame(expand.grid(plants.trait.space = max.trait.space, insect.community = prop.specialists)) %>%
  filter(insect.community > 0.8 | insect.community < 0.6)

## WHY WON'T THIS WORK WITH 0.7 and 0.6!?!?!?!

## determine insect distributions in functional trait space
plant.insect.distrib <- list()
for(i in 1:dim(plant.insect.combos)[1]){
  tmp.trait.space <- plant.insect.combos$plants.trait.space[i]
  tmp.specialists <- n.insect.species*plant.insect.combos$insect.community[i]
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
#insect.distrib.df <- plyr::ldply(insect.distrib)

## determine whether an interaction can occur or not.

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

make.webs <- lapply(web.list, FUN = special.spread)

## CALCULATE NETWORK METRICS FOR EACH WEB
web.complexity <- lapply(make.webs, FUN = function(x) networklevel(x, index = c("connectance")))

web.NODF <- lapply(make.webs, FUN = function(x) nestednodf(x)$statistic[3])

web.modularity <- lapply(make.webs, FUN = function(x) computeModules(x)@likelihood)

## merge data together
model.predictions <- data.frame(plant.insect.combos,
                                Connectance = plyr::ldply(web.complexity)$connectance,
                                NODF = plyr::ldply(web.NODF)$NODF,
                                Modularity = plyr::ldply(web.modularity)$V1)

with(model.predictions, cor.test(Modularity, NODF)) # very strong negative correlation
### got up to here. NOW I NEED TO ADD THE INSECT TRAIT DISTRIBUTIONS TO DETERMINE WHO CAN INTERACT WITH WHOM.

#plant.insect.df <- cbind.data.frame(plant.trait.distrib.df, select(insect.distrib.df, -max.trait.space)) 

## determine the distribution and range of specialist trait space
#specialist.trait.distrib <- sort(runif(n.insect.species, min = 0, max = max.trait.space))
#specialist.max <- specialist.trait.distrib + specialist.degree # not constrained within max.trait.space

## determine the distribution and range of generalist trait space
#generalist.trait.distrib <- sort(runif(n.insect.species, min = 0, max = max.trait.space))
#generalist.max <- generalist.trait.distrib + generalist.degree

## generate data frame to look at the data
#data.frame(plant.trait.distrib, specialist.trait.distrib, specialist.max)

## determine whether an interaction can occur or not.

## specialists
#spec.links <- list()
#for(i in 1:length(plant.trait.distrib)){
#  spec.links[[i]] <- ifelse(plant.trait.distrib[i] > specialist.trait.distrib & plant.trait.distrib[i] < specialist.max, 1, 0)
#}
#library(bipartite)
#spec.web <- plyr::ldply(spec.links)

# plot web
#png(filename = "specialist.web.png")
#plotweb(spec.web)
#dev.off()

# complexity
#networklevel(spec.web, index = c("connectance"))

# nestedness
#nodf.spec <- nestednodf(spec.web)
#png(filename = "specialist.nestedness.png")
#visweb(spec.web, type = "nested")
#dev.off()

# modularity
#spec.mods <- computeModules(spec.web)
#png(filename = "specialist.modules.png")
#plotModuleWeb(spec.mods)
#dev.off()

## generalists
#gen.links <- list()
#for(i in 1:length(plant.trait.distrib)){
#  gen.links[[i]] <- ifelse(plant.trait.distrib[i] > generalist.trait.distrib & plant.trait.distrib[i] < generalist.max, 1, 0)
#}
#gen.web <- plyr::ldply(gen.links)

# plot web
#png(filename = "generalist.web.png")
#plotweb(gen.web)
#dev.off()

# connectance
#networklevel(gen.web, index = c("connectance"))

# nestedness
#nodf.gen <- nestednodf(gen.web)
#png(filename = "generalist.nestedness.png")
#visweb(gen.web, type = "nested")
#dev.off()

# modularity
#gen.mods <- computeModules(gen.web)
#png(filename = "generalist.modules.png")
#plotModuleWeb(gen.mods)
#dev.off()

## Null Model Analysis ----
library(vegan)

# nestedness
nested.test <- lapply(make.webs, FUN = function(x) oecosimu(x, nestednodf, method = "quasiswap", nsimul = 1000))

# modularity
QuaBiMod <- function(x) computeModules(x)@likelihood
modularity.test <- lapply(make.webs, FUN = function(x) oecosimu(x, QuaBiMod, method = "quasiswap", nsimul = 1000)) 
