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
tictoc::tic()
numreps <- 2
d <- c(0.2, 0.5, 0.8)
num_generations <- c(10,20,50)
cross_over_nums <- c(0.01, 0.05, 0.1)
pops <- c(100,1000,10000)
y <- c(0.5,1.0,2.0)
sesss <- c(1.25,1.5,2.0)

for (k in d){
  for (i in 1:numreps){
    run_simulation(L = 50, 
                   pop_size = 100,
                   d = k,
                   y = 1, 
                   cross_prob = 0.03,
                   mut_prob = 1*10^(-4),
                   years = 1,
                   generations = 20,
                   seasonal_balance = 2,
                   rep = i)
    
  }
}
tictoc::toc()

