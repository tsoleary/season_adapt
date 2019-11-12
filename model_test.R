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
pop_size <- 100
# dominance
d <- 0.5
# exponent of the fitness function? epistasis parameter or something?
y <- 1
# recombination rate
cross_prob <- 0.03
# mutation rate
mut_prob <- 0.05 #2*10^(-4)
# duration of experiment in years
years <- 5
# number of generations in a season
generations <- 50
# balance between seasons (2 is even, less than 2 means more summer, etc.)
seasonal_balance <- 2

#site in genome = 1 = summer, 0 = winter

# Initialize Population --------------------------------------------------------

genomes <- init_pop(L, pop_size) # initial population
fitness_all <- fitness_func(genomes) # fitness at time 0
season <- "summer" # initial season

# Start experiment -------------------------------------------------------------
genomes_over_time <- vector()

new_overall_genome <- select(genomes, contains("locus"))
new_overall_genome <-  data.matrix(new_overall_genome, rownames.force = NA) #saving for comparing later
genomes_over_time[1]<-sum(new_overall_genome)

G <- 1

for (year in 1:years){
  for (generation in 1:generations){
    if (generation < generations / seasonal_balance){
      season <- "summer"
    } else {
      season <- "winter"
    }
    new_pop <- init_pop(L, pop_size) #prob_0 = 1, prob_1 = 0)
    for (i in 1:pop_size){
      selected_for_mating <- select_inds(genomes, fitness_all, season)
      
      crossed1 <- cross_over(selected_for_mating[[1]], cross_prob)
      new_ind_chr1 <- crossed1[sample(seq_len(nrow(crossed1)), 1),]
      crossed2 <- cross_over(selected_for_mating[[2]], cross_prob)
      new_ind_chr2 <- crossed2[sample(seq_len(nrow(crossed2)), 1),]
      
      loci1 <- select(new_ind_chr1, contains("locus"))
      loci2 <- select(new_ind_chr2, contains("locus"))
      
      new_pop[2*i,3:L] <- loci1 #this is not ideal, don't know how to improve
      new_pop[(2*i)-1,3:L] <- loci2
    }
    genomes <- new_pop
    new_overall_genome <- select(genomes, contains("locus"))
    new_overall_genome <-  data.matrix(new_overall_genome, rownames.force = NA) #saving for comparing later
    genomes_over_time[G]<-sum(new_overall_genome)
    fitness_all <- fitness_func(genomes)
    genomes <- mutate_genome(genomes,mut_prob) 
    print(paste("year", year, "generation", generation))
    print(season)
    G <- G + 1
    print(G)
  }
}


plot(genomes_over_time, ylab="Number of allele 1 in population", xlab="Number of generations")
lines(genomes_over_time)



