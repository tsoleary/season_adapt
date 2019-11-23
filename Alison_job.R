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

numreps <- 20
pops <- c(100,1000,10000) #Alison
y <- c(0.5,1.0,2.0) #Alison

for (k in pops){
  for (i in 1:numreps){
    run_simulation(L = 250, 
                   pop_size = k,
                   d = 0.5,
                   y = 1, 
                   cross_prob = 0.05,
                   mut_prob = 1*10^(-4),
                   years = 300,
                   generations = 20,
                   seasonal_balance = 2,
                   rep = i)
    
  }
}
for (k in y){
  for (i in 1:numreps){
    run_simulation(L = 250, 
                   pop_size = 1000,
                   d = 0.5,
                   y = k, 
                   cross_prob = 0.05,
                   mut_prob = 1*10^(-4),
                   years = 300,
                   generations = 20,
                   seasonal_balance = 2,
                   rep = i)
    
  }
}


#default parameters
# 
# (L = 250, 
#   pop_size = 1000,
#   d = 0.5,
#   y = 1, 
#   cross_prob = 0.05,
#   mut_prob = 1*10^(-4),
#   years = 300,
#   generations = 20,
#   seasonal_balance = 2,
#   rep = i)

