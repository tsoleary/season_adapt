# seasonal adaptation modeling -------------------------------------------------

require(tidyverse)

# number of loci in a genome
L <- 50
# population size 
pop_size <- 1
# dominance
d <- 0.5
# exponent of the fitness function? epistasis parameter or something?
y = 1


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

# long format for dplyr
gen_tidy <- genome %>%
  pivot_longer(cols = contains("locus"), 
               names_to = "loci", 
               values_to = "allele")

# sumarize each locus for each inidividual 
genotype <- gen_tidy %>%
  group_by(individual, loci) %>%
  summarize(geno = sum(as.numeric(allele)))

# count of each genotype
geno_count <- genotype %>%
  group_by(individual, geno) %>%
  count()

# seasonal score for each individual
z_score <- geno_count %>%
  pivot_wider(names_from = geno, values_from = n) %>%
  rename("geno_0" = "0", "geno_1" = "1", "geno_2" = "2") %>%
  mutate(z_winter = geno_0 + geno_1 * d,
         z_summer = geno_2 + geno_1 * d)

# fitness of each individual
fitness <- z_score %>%
  mutate(f_winter = (1 + z_winter) ^ y,
         f_summer = (1 + z_summer) ^ y)

# selection of the parents 
# tournament
# totally random panmictic mating ???



# crossover 
cross_prob <- 0.2
runif(L) < cross_prob
select(genome, contains("locus"))


  
  


