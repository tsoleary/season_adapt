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
# exponent of the fitness function
y <- 1
# recombination rate
cross_prob <- 0.03
# mutation rate
mut_prob <- 1*10^(-4)
# duration of experiment in years
years <- 20
# number of generations in a season
generations <- 20
# balance between seasons (2 is even, less than 2 means more summer, etc.)
seasonal_balance <- 2

#site in genome = 1 = summer, 0 = winter

start <- proc.time()

# Initialize Population --------------------------------------------------------

genomes <- init_pop(L, pop_size) 

# Start experiment -------------------------------------------------------------
genomes_over_time <- vector()

new_overall_genome <- select(genomes, contains("locus"))
new_overall_genome <-  data.matrix(new_overall_genome, rownames.force = NA)
genomes_over_time[1] <- sum(new_overall_genome)

G <- 1

for (year in 1:years){
  for (generation in 1:generations){
    if (generation < generations / seasonal_balance){
      season <- "summer"
    } else {
      season <- "winter"
    }
    # create an empty population data frame with all zeros
    new_pop <- init_pop(L, pop_size, prob_0 = 1, prob_1 = 0)
    for (i in 1:pop_size){
      # select parents based on fitness
      fitness_all <- fitness_func(genomes)
      selected_for_mating <- select_inds(genomes, fitness_all, season)
      # cross chromosomes of parents
      crossed1 <- cross_over(selected_for_mating[[1]], cross_prob)
      new_ind_chr1 <- crossed1[sample(seq_len(nrow(crossed1)), 1),]
      crossed2 <- cross_over(selected_for_mating[[2]], cross_prob)
      new_ind_chr2 <- crossed2[sample(seq_len(nrow(crossed2)), 1),]
      
      loci1 <- select(new_ind_chr1, contains("locus"))
      loci2 <- select(new_ind_chr2, contains("locus"))
      
      new_pop[2*i, which(grepl("locus", colnames(new_pop)))] <- loci1 
      new_pop[(2*i)-1, which(grepl("locus", colnames(new_pop)))] <- loci2
    }
    genomes <- new_pop
    
    # get a sum of all 1 alleles in the population for each generation
    new_overall_genome <- select(genomes, contains("locus"))
    new_overall_genome <-  data.matrix(new_overall_genome, rownames.force = NA) 
    genomes_over_time[G]<-sum(new_overall_genome)

    
    # mutate genomes for the next year
    genomes <- mutate_genome(genomes, mut_prob)
    
    # print information to keep track of simulation progress
    print(paste("year", year, "generation", generation, "season", season,
                "total generation", G))
    G <- G + 1
  }
}

time_elapsed <- proc.time() - start

plot(genomes_over_time, ylab = "Number of allele 1 in population", 
     xlab = "Number of generations")
lines(genomes_over_time)



