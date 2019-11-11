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
#duration of experiment in years
years <- 1
#number of generations in a season
generations <- 20

#site in genome = 0 = summer, 1 = winter
#time of the year = 0 = summer, 1 = winter

# Initialize Population --------------------------------------------------------

genomes <- init_pop(L, pop_size) # initial population
fitness_all <- fitness_func(genomes) # fitness at time 0
season <- "summer" # initial season


mutated_genome <- mutate_genome(genomes,mut_prob) 


# Start experiment -------------------------------------------------------------
for (year in 1:years){
  for (generation in 1:generations){
    for (i in 1:pop_size){
      selected_for_mating <- select_inds(genomes, fitness_all, season)
      crossed1<-cross_over(selected_for_mating[[1]], cross_prob)
      new_ind_chr1 <- crossed1[sample(seq_len(nrow(crossed1)), 1),]
      crossed2 <- cross_over(selected_for_mating[[2]], cross_prob)
      new_ind_chr2 <- crossed2[sample(seq_len(nrow(crossed2)), 1),]
      new_ind_name <- i
    }
    print("one generation has passed")
  }
  print("one year has passed")
}



