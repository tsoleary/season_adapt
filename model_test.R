# seasonal adaptation modeling -------------------------------------------------

require(stringr) 
require(tibble)
require(tidyr)
require(dplyr)
require(ggplot2)

source("functions.R")

# Parameter definitions --------------------------------------------------------
# L is the number of loci in a genome
# pop_size is the population size 
# d is the dominance parameter of the heterozygote
# y is the exponent of the fitness function
# cross_prob is the recombination rate
# mut_prob is the mutation rate
# years is the duration of experiment in years
# generations is the number of generations in a season
# seasonal_balance is the balance between seasons
#      (2 = even; < 2 more summer; > 2 more winter)

# Run a simulation -------------------------------------------------------------
sim_results <- run_simulation(L = 50, 
                              pop_size = 100,
                              d = 0.5,
                              y = 1, 
                              cross_prob = 0.03,
                              mut_prob = 1*10^(-4),
                              years = 2,
                              generations = 5,
                              seasonal_balance = 2)




numreps <- 1
d <- c(0.2, 0.5, 0.7)
sim_results <- sim()
for (i in 1:numreps){
  for (k in 1:d){
run_simulation(L = 5, 
                                   pop_size = 10,
                                   d = k,
                                   y = 1, 
                                   cross_prob = 0.03,
                                   mut_prob = 1*10^(-4),
                                   years = 2,
                                   generations = 5,
                                   seasonal_balance = 2)
  
  }
}
