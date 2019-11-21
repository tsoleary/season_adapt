# seasonal adaptation modeling -------------------------------------------------

require(tidyverse)

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
sim_results <- run_simulation(L = 50, 
                              pop_size = 100,
                              d = 0.5,
                              y = 1, 
                              cross_prob = 0.03,
                              mut_prob = 1*10^(-4),
                              years = 1,
                              generations = 20,
                              seasonal_balance = 2)
tictoc::toc()
sim_clean <- sim_results %>%
  filter(genz %% 10 == 0)

plot_freq(sim_clean, figure_caption = "")
plot_freq(sim_clean, type = "avg", figure_caption = "")


