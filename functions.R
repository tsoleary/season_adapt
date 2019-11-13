# Seasonal adaptation model functions ------------------------------------------

# Function to initialise population --------------------------------------------
init_pop <- function(L, pop_size, prob_0 = 0.5, prob_1 = 0.5) {
  # initialize diploid chromosomes with a random population
  genome <- matrix(sample(0:1, 2*L*pop_size, 
                          replace = TRUE, 
                          prob = c(prob_0, prob_1)), 
                   nrow = pop_size * 2, 
                   ncol = L)
  
  # name the loci, the individuals, and their chromosome pairs
  colnames(genome) <- paste("locus", 
                            str_pad(1:L, 2, "0", side = "left"), sep = "_")
  individual <- paste("indiv", 
                      str_pad(rep(1:pop_size, each = 2), 6, "0", side = "left"), 
                      sep = "_")
  chr <- paste("chr", rep(1:2, times = pop_size), sep = "_")
  
  # bind the matrix and names together
  return(as_tibble(cbind(individual, chr, genome)))
}


# Funtion to calculate fitness -------------------------------------------------
fitness_func <- function(genomes) {
  # creates the diploid genotypes: 0 (homo 00), 1 (het), or 2 (homo 11)
  genotype <- genomes %>%
    pivot_longer(cols = contains("locus"), 
                 names_to = "loci", 
                 values_to = "allele") %>%
    group_by(individual, loci) %>%
    summarize(geno = sum(as.numeric(allele)))
  
  # count of each genotype for each loci in an individual
  geno_count <- genotype %>%
    group_by(individual, geno, .drop = FALSE) %>%
    count()
  
  # seasonal score for each individual
  z_score <- geno_count %>%
    pivot_wider(names_from = geno, values_from = n) %>%
    rename("geno_0" = "0", "geno_1" = "1", "geno_2" = "2") %>%
    mutate(z_winter = geno_0, # + geno_1 * d,
           z_summer = geno_2) # + geno_1 * d)
  
  # heterozygote effect
  z_score[is.na(z_score)] <- 0
  z_score["z_winter"] <- z_score["z_winter"] + z_score["geno_1"]*d
  z_score["z_summer"] <- z_score["z_summer"] + z_score["geno_1"]*d
  
  # fitness of each individual
  fitness <- z_score %>%
    mutate(f_winter = (1 + z_winter) ^ y,
           f_summer = (1 + z_summer) ^ y)
  
  return(fitness)
}

# Selecting individuals for the next generation, depends on season -------------

select_inds <- function(genomes,fitness_all, season) {
  if (season == "summer") { # select 2 parents at random weighted by fitness
    selected_idx <- sample(seq_len(nrow(fitness_all)), 2, 
                           prob = fitness_all$f_summer, replace = FALSE)
  } else {
    selected_idx <- sample(seq_len(nrow(fitness_all)), 2, 
                           prob = fitness_all$f_winter, replace = FALSE)
  }
  temp1 <- selected_idx*2-1
  temp2 <-selected_idx*2
  ind1 <- genomes[temp1[1]:temp2[1],] #because we want both chrs
  ind2 <- genomes[temp1[2]:temp2[2],]
  selected_individuals <- list(ind1,ind2)
  return(selected_individuals)
}

# Crossover function that takes a single individual ----------------------------
cross_over <- function(ind, cross_prob) {
  # select only the loci or only the indiv and chr identifiers
  loci <- select(ind, contains("locus"))
  indiv_chr <- select(ind, -contains("locus"))
  
  # determine the location(s) for crossover
  cross <- sample(0:1, 
                  ncol(loci), 
                  prob = c(1 - cross_prob, cross_prob), 
                  replace = TRUE)
  cross_locs <- which(cross == 1)
  
  # do crossover at each location 
  for (cross_loc in cross_locs){
    cross_genome <- loci[,cross_loc:ncol(loci)]
    cross_genome[c(1,2), ] <- cross_genome[c(2,1), ]
    loci[,cross_loc:ncol(loci)] <- cross_genome
  }
  
  # add back the indiv and chr identifiers
  return(as_tibble(cbind(indiv_chr, loci)))
}

# Mutation on the entire population each generation before crossover -----------
mutate_genome <- function(genomes, mut_prob) {
  # multiply the mutation probabilty by two (because half the time 0's will 
  # be mutated to 0's and same for 1's)
  mut_prob_2 <- mut_prob * 2
  
  # total number of chromosomes and loci in the genome
  total_chr <- nrow(genomes)
  genome_length <- sum(grepl("locus", colnames(genomes)))
  
  # determine mutation positions with a random uniform distribution
  mut_pos <- matrix(runif(genome_length * total_chr), 
                         nrow = total_chr, 
                         ncol = genome_length) < mut_prob_2
  
  # mutation values (with zeros where the original values are located)
  mut_mat <- matrix(sample(0:1, genome_length, replace = TRUE), 
                     nrow = total_chr, 
                     ncol = genome_length) * mut_pos
  
  # original values (with zeros where the mutated values are located)
  org_mat <- matrix(as.numeric(as.matrix(select(genomes, contains("locus")))),
                nrow = total_chr, 
                ncol = genome_length) * !mut_pos

  # combine matrices by addition
  new_mat <- mut_mat + org_mat
  
  # add indiv and chr identifiers to dataframe
  indiv_chr <- select(genomes, -contains("locus"))
  df <- cbind(indiv_chr, new_mat)
  
  # keep column names the same for the loci
  colnames(df) <- colnames(genomes)
  
  return(as_tibble(df))
}

# Calculate loci specific frequencies 
get_freqs <- function(genomes, G){
  genomes %>%
    pivot_longer(cols = contains("locus"), 
                 names_to = "loci", 
                 values_to = "allele") %>%
    group_by(loci) %>%
    summarize(freq_1 = sum(as.numeric(allele)) / (pop_size * 2))
}


