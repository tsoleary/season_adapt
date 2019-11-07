# seasonal adaptation modeling -------------------------------------------------

require(tidyverse)
#couldn't install tidyverse so loaded stuff separately
require(stringr) 
require(tibble)
require(tidyr)
require(dplyr)

#Parameters
# number of loci in a genome
L <- 50
# population size 
pop_size <- 10
# dominance
d <- 0.5
# exponent of the fitness function? epistasis parameter or something?
y = 1
# recombination rate
r <- 1
#which season we are in, by default start in summer
season <- "summer"
#how many people are selected to reproduce to next generation
S<-5 # could cause bottlenecks if this value is low

#site in genome = 0 = summer, 1 = winter
#time of the year = 0 = summer, 1 = winter

#Function to initialise population

init_pop <- function(L, pop_size) {
  # initialize diploid chromosomes with a random population with equal allele freq
  genome <- matrix(sample(0:1, 2*L*pop_size, replace = TRUE, prob = c(0.5, 0.5)), 
                   nrow = pop_size*2, ncol = L)
  loci <- paste("locus", str_pad(1:L, 2, "0", side = "left"), sep = "_")
  colnames(genome) <- loci
  # name the individuals and their chromosome pairs for dplyr later
  individual <- paste("indiv", 
                      str_pad(rep(1:pop_size, each = 2), 6, "0", side = "left"), 
                      sep = "_")
  chr <- paste("chr", rep(1:2, times = pop_size), sep = "_")
  # bind the matrix and names together
  genome <- as.tibble(cbind(individual, chr, genome))
  return(genome)
}

genomes <- init_pop(L,pop_size) # gives you L-2 loci!

# Funtion to calculate fitness

fitness_func <- function(genomes) {
  # long format for dplyr
  gen_tidy <- genomes %>%
    pivot_longer(cols = contains("locus"), 
                 names_to = "loci", 
                 values_to = "allele")
  # sumarize each locus for each inidividual across the two chromosomes
  genotype <- gen_tidy %>%
    group_by(individual, loci) %>%
    summarize(geno = sum(as.numeric(allele)))
  
  # count of each genotype, genotype (0,1 or 2) for each loci in an individual
  geno_count <- genotype %>%
    group_by(individual, geno) %>%
    count()
  # seasonal score for each individual
  z_score <- geno_count %>%
    pivot_wider(names_from = geno, values_from = n) %>%
    rename("geno_0" = "0", "geno_1" = "1", "geno_2" = "2") %>%
    mutate(z_winter = geno_0 + geno_1 * d,
           z_summer = geno_2 + geno_1 * d)
  # fitness of each individual, based on equation in paper
  fitness <- z_score %>%
    mutate(f_winter = (1 + z_winter) ^ y,
           f_summer = (1 + z_summer) ^ y)
  return(fitness)
}

fitness=fitness_func(genomes)


# Selecting individuals for the next generation, depends on season
curr_season <- season
if (curr_season == "summer") { # select 2 parents at random weighted by fitness
  selected_idx <- sample(seq_len(nrow(fitness)), 2, prob=fitness$f_summer, replace = FALSE)
} else {
  selected_idx <- sample(seq_len(nrow(fitness)), 2, prob=fitness$f_winter, replace = FALSE)  
}

# selection of the parents 
# tournament
# totally random panmictic mating ???


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





