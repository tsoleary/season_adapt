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
fitness_func <- function(genomes, d, y) {
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
select_inds <- function(genomes, fitness_all, season) {
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
cross_over <- function(parent) {
  loci <- genomes %>%
    filter(.$individual == as.character(parent)) %>%
    select(contains("locus"))
  
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
  
  pick_chr <- sample(c(1,2), 1)
  
  # add back the indiv and chr identifiers
  return(loci[pick_chr, ])
}

# cross over parents at the same time ------------------------------------------
cross_over_parents <- function(parent1, parent2) {
  
  p1_xover <- cross_over(parent1)#, genomes, cross_prob)
  p2_xover <- cross_over(parent2)#, genomes, cross_prob)
  
  new_chr <- rbind(p1_xover, p2_xover)
  bind_cols("chr" = c("chr_1", "chr_2"), new_chr)
  
}

# Parent selection -------------------------------------------------------------
parent_selection <- function(genomes, fitness_all, season){
  if (season == "summer"){
    fit <- fitness_all$f_summer
  } else {
    fit <- fitness_all$f_winter
  }
  # create an empty data frame for parent 1 and parent 2
  df <- tibble(p1 = character(),
               p2 = character()) 
  # sample pairs of potential parents based on fitness
  for (i in 1:nrow(fitness_all)){
    df[i, 1:2] <- sample(fitness_all$individual, 2, 
                         prob = fit,
                         replace = FALSE)
  }
  return(df)
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

# Calculate loci specific frequencies ------------------------------------------
get_freqs <- function(genomes, pop_size){
  genomes %>%
    pivot_longer(cols = contains("locus"), 
                 names_to = "loci", 
                 values_to = "allele") %>%
    group_by(loci) %>%
    summarize(freq_1 = sum(as.numeric(allele)) / (pop_size * 2))
}

# Plotting allele frequencies overtime -----------------------------------------
plot_freq <- function(loci_freq, type = "loci", figure_caption) {
  if (type == "avg"){
    avg_freq <- loci_freq %>%
      group_by(genz) %>%
      summarize(freqs  = mean(freqs))
    
    g <- ggplot(avg_freq, mapping = aes(x = genz, y = freqs)) + 
      geom_line() +
      labs(title = "Total summer allele frequency over time",
           caption = figure_caption) +
      xlab("Generations") +
      ylab("Freq of Summer Allele") + 
      ylim(0,1) +
      theme_classic() +
      theme(legend.position = "none")
  } else {
    plot_df <- loci_freq
    g <- ggplot(plot_df, mapping = aes(x = genz, y = freqs, color = loci)) + 
      geom_line() +
      labs(title = "Loci specific allele frequencies over time",
           caption = figure_caption)+
      xlab("Generations") +
      ylab("Freq of Summer Allele") + 
      ylim(0,1) +
      theme_classic() +
      theme(legend.position = "none")
  }
  return(g)
}

# Simulation -------------------------------------------------------------------
run_simulation <- function(L, pop_size, d, y, cross_prob, mut_prob, years, 
                           generations, seasonal_balance) {
  # Initialize Population 
  
  genomes <- init_pop(L, pop_size) 
  freq_df <- get_freqs(genomes, pop_size)
  
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
      fitness_all <- fitness_func(genomes, d, y)
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
      
      freq_temp <- get_freqs(genomes, pop_size)
      freq_df <- full_join(freq_df, freq_temp, by = "loci")
      colnames(freq_df)[which(colnames(freq_df) == "freq_1")] <- 
        paste0("freq_G.", G)
      
      # mutate genomes for the next year
      genomes <- mutate_genome(genomes, mut_prob)
      
      # print information to keep track of simulation progress
      print(paste("year", year, "generation", generation, "season", season,
                  "total generation", G))
      G <- G + 1
    }
  }
  
  colnames(freq_df)[2:3] <- c("freq_G.0", "freq_G.1")
  
  loci_freq <- freq_df %>% 
    pivot_longer(cols = contains("freq"), 
                 names_to = "genz", 
                 values_to = "freqs")
  loci_freq$genz <- as.numeric(str_extract(loci_freq$genz, "[:digit:]+"))
  
  caption <- paste(paste("Generations per year", generations),
                   paste("Pop size", pop_size),
                   paste("Seasonal Balance", seasonal_balance),
                   paste("Number of Loci", L),
                   paste("Dominance", d),
                   paste("Epistasis", y), sep = "; ")
  
  g1 <- plot_freq(loci_freq, figure_caption = caption)
  print(g1)
  g2 <- plot_freq(loci_freq, type = "avg", figure_caption = caption)
  print(g2)
  
  return(loci_freq)
}

# Simulation -------------------------------------------------------------------
run_simulation <- function(L, pop_size, d, y, cross_prob, mut_prob, years, 
                           generations, seasonal_balance, rep) {
  # Initialize Population 
  
  individual <- paste("indiv", 
                      str_pad(rep(1:pop_size, each = 2), 
                              6, "0", side = "left"), 
                      sep = "_")
  
  cross_prob <<- cross_prob
  genomes <<- init_pop(L, pop_size) 
  freq_df <- get_freqs(genomes, pop_size)
  
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
      fitness_all <- fitness_func(genomes, d, y)
      
      df <- parent_selection(genomes, fitness_all, season)
      
      new_pop <- map2_df(df[[1]], df[[2]], cross_over_parents)
      
      genomes <- cbind(individual, new_pop)
      
      freq_temp <- get_freqs(genomes, pop_size)
      freq_df <- full_join(freq_df, freq_temp, by = "loci")
      colnames(freq_df)[which(colnames(freq_df) == "freq_1")] <- 
        paste0("freq_G.", G)
      
      # mutate genomes for the next year
      genomes <<- mutate_genome(genomes, mut_prob)
      
      # print information to keep track of simulation progress
      print(paste("year", year, "generation", generation, "season", season,
                  "total generation", G))
      G <- G + 1
    }
  }
  
  colnames(freq_df)[2:3] <- c("freq_G.0", "freq_G.1")
  
  loci_freq <- freq_df %>% 
    pivot_longer(cols = contains("freq"), 
                 names_to = "genz", 
                 values_to = "freqs")
  loci_freq$genz <- as.numeric(str_extract(loci_freq$genz, "[:digit:]+"))
  
  caption <<- paste(paste("Generations per year", generations),
                    paste("Pop size", pop_size),
                    paste("Seasonal Balance", seasonal_balance),
                    paste("Number of Loci", L),
                    paste("Dominance", d),
                    paste("Epistasis", y), sep = "; ")
  
  g1 <- plot_freq(loci_freq, figure_caption = caption)
  print(g1)
  g2 <- plot_freq(loci_freq, type = "avg", figure_caption = caption)
  print(g2)
  
  file_names <- paste("results",
                      "G", generations, 
                      "Ps", pop_size,
                      "Sb", seasonal_balance,
                      "L", L,
                      "d", d,
                      "y", y,
                      "c", cross_prob,
                      "num_rep", rep,
                      sep = "_")
  
  write.csv(loci_freq, paste0(file_names, ".csv"), row.names = FALSE)
  
  return(loci_freq)
}


