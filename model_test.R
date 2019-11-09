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

#site in genome = 0 = summer, 1 = winter
#time of the year = 0 = summer, 1 = winter

# Initialize Population --------------------------------------------------------

genomes <- init_pop(L, pop_size)

fitness_all <- fitness_func(genomes)

ind1 <- genomes[1:2,]
ind2 <- genomes[3:4,]
#crossover takes individuals in this format
ind_test<-cross_over(ind1, cross_prob)

# mutation on the entire population each generation before crossover
mut_prob <- 0.05 #1*10^(-4)

mut_test_mat <- matrix(runif(L*pop_size*2), 
                       nrow = pop_size * 2, 
                       ncol = L)
mut_pos <- mut_test_mat < mut_prob

mut_vals <- matrix(sample(0:1, 50, replace = TRUE), 
                   nrow = pop_size * 2, 
                   ncol = L)

mut_mat <- mut_vals * mut_pos

mat <- matrix(as.numeric(as.matrix(select(genome, contains("locus")))),
              nrow = pop_size * 2,
              ncol = L)

par_mat <- mat * !mut_pos

new_mat <- mut_mat + par_mat


# crossover 
# this cross over requires only one parent at a time
# would it be possible to do all at once, maybe with dplyr real quick??
cross_prob <- 0.2
cross <- sample(0:1, pop_size, prob = c(1 - cross_prob, cross_prob))
cross_loc <- sample(2:L, 1)
cross_genome <- select(genome, contains("locus"))[1:2, cross_loc:L]
cross_genome[c(1,2), ] <- cross_genome[c(2,1), ]
genome[1:2, (length(genome) - length(cross_loc:L) + 1):length(genome)] <- cross_genome





