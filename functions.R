# functions for seasonal adaptation model --------------------------------------

# Function to initialise population --------------------------------------------
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


# Funtion to calculate fitness -------------------------------------------------
fitness_func <- function(genomes) {
  # long format for dplyr
  gen_tidy <- genomes %>%
    pivot_longer(cols = contains("locus"), 
                 names_to = "loci", 
                 values_to = "allele")
  # sumarize each locus for each inidividual across the two chromosomes=
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

# Selecting individuals for the next generation, depends on season

select_inds <- function(genomes,fitness_all, season) {
  if (season == "summer") { # select 2 parents at random weighted by fitness
    selected_idx <- sample(seq_len(nrow(fitness_all)), 2, prob=fitness_all$f_summer, replace = FALSE)
  } else {
    selected_idx <- sample(seq_len(nrow(fitness_all)), 2, prob=fitness_all$f_winter, replace = FALSE)
  }
  temp1 <- selected_idx*2-1
  temp2 <-selected_idx*2
  ind1 <- genomes[temp1[1]:temp2[1],] #because we want both chrs
  ind2 <- genomes[temp1[2]:temp2[2],]
  selected_individuals <- list(ind1,ind2)
  return(selected_individuals)
}

# takes individual with its 2 chromosomes, does crossover ----------------------
cross_over <- function(ind, cross_prob) {
  stripped <- select(ind, contains("locus"))
  strip_names <- select(ind, -contains("locus"))
  cross <- sample(0:1, L, prob = c(1 - cross_prob, cross_prob), replace = TRUE)
  cross_locs <- which(cross == 1)
  for (cross_loc in cross_locs){
    cross_genome <- stripped[,cross_loc:L]
    cross_genome[c(1,2), ] <- cross_genome[c(2,1), ]
    stripped[,cross_loc:L] <- cross_genome
  }
  return(cbind(strip_names, stripped))
}

# mutation on the entire population each generation before crossover -----------
mutate_genome <- function(genomes,mut_prob) {
  mut_test_mat <- matrix(runif(L*pop_size*2), 
                         nrow = pop_size * 2, 
                         ncol = L)
  mut_pos <- mut_test_mat < mut_prob
  
  mut_vals <- matrix(sample(0:1, 50, replace = TRUE), 
                     nrow = pop_size * 2, 
                     ncol = L)
  
  mut_mat <- mut_vals * mut_pos
  
  mat <- matrix(as.numeric(as.matrix(select(genomes, contains("locus")))),
                nrow = pop_size * 2,
                ncol = L)

  par_mat <- mat * !mut_pos
  
  new_mat <- mut_mat + par_mat
  strip_names <- select(genomes, -contains("locus"))
  df <- cbind(strip_names, new_mat)
  colnames(df) <- colnames(genomes)
  return(df)
}




