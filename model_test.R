# seasonal adaptation modeling -------------------------------------------------

require(stringr) 
require(tibble)
require(tidyr)
require(dplyr)
require(ggplot2)

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
years <- 2
# number of generations in a season
generations <- 50
# balance between seasons (2 is even, less than 2 means more summer, etc.)
seasonal_balance <- 2

#site in genome = 1 = summer, 0 = winter

# Initialize Population --------------------------------------------------------

genomes <- init_pop(L, pop_size) 
freq_df <- get_freqs(genomes)

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
    fitness_all <- fitness_func(genomes)
    for (i in 1:pop_size){
      # select parents based on fitness
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
    
    freq_temp <- get_freqs(genomes)
    freq_df <- full_join(freq_df, freq_temp, by = "loci")
    colnames(freq_df)[which(colnames(freq_df) == "freq_1")] <- 
      paste0("freq_G.", G)
    
    # # get a sum of all 1 alleles in the population for each generation
    # new_overall_genome <- select(genomes, contains("locus"))
    # new_overall_genome <-  data.matrix(new_overall_genome, rownames.force = NA) 
    # genomes_over_time[G]<-sum(new_overall_genome)

    
    # mutate genomes for the next year
    genomes <- mutate_genome(genomes, mut_prob)
    
    # print information to keep track of simulation progress
    print(paste("year", year, "generation", generation, "season", season,
                "total generation", G))
    G <- G + 1
  }
}

colnames(freq_df)[2:3] <- c("freq_G.0", "freq_G.1")

x <- freq_df %>% 
  pivot_longer(cols = contains("freq"), 
               names_to = "genz", 
               values_to = "freqs")
x$genz <- as.numeric(str_extract(x$genz, "[:digit:]+"))

# plotting
ggplot(x, mapping = aes(x = genz, y = freqs, color = loci)) + 
  geom_line() +
  labs(title = "Loci specific allele frequencies over time",
       caption = paste(paste("Generations per season", generations),
                        paste("Pop size", pop_size),
                        paste("Seasonal Balance", seasonal_balance),
                        paste("Number of Loci", L),
                        paste("Dominance", d),
                        paste("Epistasis", y), sep = "; ")) +
  xlab("Generations") +
  ylab("Freq of Summer Allele") + 
  theme_classic() +
  theme(legend.position = "none")




