# seasonal adaptation modeling -------------------------------------------------

require(stringr) 
require(tibble)
require(tidyr)
require(dplyr)

source("functions.R")

# Parameters -------------------------------------------------------------------
# number of loci in a genome
L <- 50
# population size 
pop_size <- 10
# dominance
d <- 0.5
# exponent of the fitness function? epistasis parameter or something?
y <- 1
# recombination rate
cross_prob <- 0.01
#mutation rate
mut_prob <- 0.05 #2*10^(-4)

#site in genome = 0 = summer, 1 = winter
#time of the year = 0 = summer, 1 = winter

# Initialize Population --------------------------------------------------------

genomes <- init_pop(L, pop_size)
fitness_all <- fitness_func(genomes)
season <- "summer"
selected_for_mating <- select_inds(genomes, fitness_all, season)
mutated_genome <- mutate_genome(genomes,mut_prob) 
#mutation
ind3 <- genomes[1:2,]
ind2 <- genomes[3:4,]
#crossover takes individuals in this format
ind_test<-cross_over(selected_for_mating[[1]], cross_prob)




