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


# uniform crossover function ---------------------------------------------------
cross_over_uniform <- function(parent) {
  loci <- genomes %>%
    filter(.$individual == as.character(parent)) %>%
    select(contains("locus"))
  
  # determine the location(s) for crossover
  cross <- sample(1:2, 
                  ncol(loci), 
                  prob = c(0.5, 0.5), 
                  replace = TRUE)
  
  new_chr <- vector(mode = "double", length = ncol(loci))
  for (i in 1:length(cross)){
    new_chr[i] <- as.numeric(loci[cross[i], i])
  }
  
  names(new_chr) <- names(loci)
  
  return(new_chr)
}

# cross over parents at the same time ------------------------------------------
cross_over_parents_uniform <- function(parent1, parent2) {
  
  p1_xover <- cross_over_uniform(parent1)
  p2_xover <- cross_over_uniform(parent2)
  
  new_chr <- as_tibble(rbind(p1_xover, p2_xover))
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
plot_freq <- function(loci_freq, type = "loci", figure_caption,
                      fig_title = "Loci specific allele frequencies over time"){
  if (type == "avg"){
    avg_freq <- loci_freq %>%
      group_by(genz) %>%
      summarize(freqs  = mean(freqs))
    
    g <- ggplot(avg_freq, mapping = aes(x = genz, y = freqs)) + 
      geom_line() +
      labs(title = fig_title,
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
      labs(title = fig_title,
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
  
  # caption <<- paste(paste("Generations per year", generations),
  #                   paste("Pop size", pop_size),
  #                   paste("Seasonal Balance", seasonal_balance),
  #                   paste("Number of Loci", L),
  #                   paste("Dominance", d),
  #                   paste("Epistasis", y), sep = "; ")
  # 
  # g1 <- plot_freq(loci_freq, figure_caption = caption)
  # print(g1)
  # g2 <- plot_freq(loci_freq, type = "avg", figure_caption = caption)
  # print(g2)
  
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
  
  #return(loci_freq)
}



# simulation with uniform crossover --------------------------------------------
run_simulation_uniform <- function(L, pop_size, d, y, cross_prob, mut_prob, 
                                   years, generations, seasonal_balance, rep) {
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
      
      new_pop <- map2_df(df[[1]], df[[2]], cross_over_parents_uniform)
      
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
  
  # caption <<- paste(paste("Generations per year", generations),
  #                   paste("Pop size", pop_size),
  #                   paste("Seasonal Balance", seasonal_balance),
  #                   paste("Number of Loci", L),
  #                   paste("Dominance", d),
  #                   paste("Epistasis", y), sep = "; ")
  # 
  # g1 <- plot_freq(loci_freq, figure_caption = caption)
  # print(g1)
  # g2 <- plot_freq(loci_freq, type = "avg", figure_caption = caption)
  # print(g2)
  
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
  
  #return(loci_freq)
}

# analysis & ANOVA plotting function -------------------------------------------
do_analysis <- function(values, test_name, boxplot_n){
  
  files = c()
  for (value in values){
    file <- list.files()[grep(value, list.files())]
    files <- list.append(files, file) #contains all the files regarding your parameter in the correct order
  }
  
  sim_results <- vector(mode = "list", length = length(values))
  i = 1
  for (f in files){
    sim_result <- read.csv(f)
    sim_results[i] <- list(sim_result) #contains all dataframes of your parameter
    i <- i+1
  }
  
  maximums <- c()
  for (rep in sim_results){
    maximums <- list.append(maximums, max(rep$genz)) #this ensures that the minimum of all the maximum generations is used for analysis
  }
  
#------------------------------------------------------------------------------------------------------
  
  means <- c()
  standard_devs <- c()
  for (rep in sim_results){
    temp <- rep[rep$genz == min(maximums),] #check that all have same max generation!
    means <- list.append(means, mean(temp$freqs))
    standard_devs <- list.append(standard_devs, sd(temp$freqs))
  }
  
  lm <- length(means)/5 #number of parameters to compare
  t <- list(1:lm)
  
#Calculating means--------------------------------------------------------------------------------------
  
  x <- data.frame("Categories" = sort(rep.int(as.numeric(unlist(t)),5)),
                  "Means" = means)
  
  ANOVAofmeans <- aov(x$Means ~ x$Categories)
  print(summary(ANOVAofmeans))
  
  v_names = c()
  p_colors = c("gray", "blue", "green", "yellow", "red", "orange", "purple", "white", "black", "pink")
  v_colors = c()
  c <- 1
  for (v in values){
    v_names <- list.append(v_names, v)
    v_colors <- list.append(v_colors, p_colors[c])
    c <- c + 1
  }
  
  #b1 <- ggplot(x, aes(x=test_name, y="Mean loci frequency")) + 
    geom_boxplot()
  
  boxplot(x$Means ~ x$Categories,
          ylab="Mean loci frequency",
          xlab= test_name,
          names= v_names,
          col= v_colors
  )
  b1 <- recordPlot()
  
  
#Calculating standard deviation-------------------------------------------------------------------------
  
  y <- data.frame("Categories" = sort(rep.int(as.numeric(unlist(t)),5)), "SD" = standard_devs)
  
  ANOVAofsds <- aov(y$SD ~ y$Categories)
  print(summary(ANOVAofsds))
  
  boxplot(y$SD ~ y$Categories,
          ylab="Standard deviation of loci frequencies",
          xlab= test_name,
          names= v_names,
          col= v_colors
  )
  b2 <- recordPlot()
  
#Calculating fixed proportions---------------------------------------------------------------------------
  
  fixed <- c()
  for (rep in sim_results){
    temp <- rep[rep$genz == min(maximums),] #check that all have same max!
    temp <- round(temp$freqs, digits = 1) #round to 1 digit so that 0.99 would be considered fixed
    fixed <- list.append(fixed, length(which(temp==0 | temp == 1))/length(temp))
  }
  
  z <- data.frame("Categories" = sort(rep.int(as.numeric(unlist(t)),5)), "Fix" = fixed)
  
  ANOVAoffixs <- aov(z$Fix ~ z$Categories)
  print(summary(ANOVAoffixs))
  
  boxplot(z$Fix ~ z$Categories,
          ylab="Proportion of fixed loci",
          xlab= test_name,
          names= v_names,
          col= v_colors
  )
  b3 <- recordPlot()
  
  if (boxplot_n == "b1"){
    return(b1)
  }
  else if (boxplot_n == "b2"){
    return(b2)
  }
  else {
    return(b3)
  }
}

# read csv and save file name for tbl in map_df --------------------------------
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

