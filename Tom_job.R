# seasonal adaptation modeling -------------------------------------------------

require(stringr) 
require(tibble)
require(tidyr)
require(dplyr)
require(ggplot2)
require(purrr)

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

numreps <- 2
d <- c(0.2, 0.5, 0.8) #Tom
num_generations <- c(10,20,50) #Tom

for (k in d){
  for (i in 1:numreps){
    run_simulation(L = 100, 
                   pop_size = 500,
                   d = k,
                   y = 1, 
                   cross_prob = 0.05,
                   mut_prob = 1*10^(-4),
                   years = 100,
                   generations = 20,
                   seasonal_balance = 2,
                   rep = i)
    
  }
}
for (k in num_generations){
  for (i in 1:numreps){
    run_simulation(L = 100, 
                   pop_size = 500,
                   d = 0.5,
                   y = 1, 
                   cross_prob = 0.05,
                   mut_prob = 1*10^(-4),
                   years = 100,
                   generations = k,
                   seasonal_balance = 2,
                   rep = i)
    
  }
}


