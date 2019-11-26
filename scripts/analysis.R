# results analysis -------------------------------------------------------------

require(tidyverse)

source(here::here("functions.R"))

# set working directory to the results folder of interest
setwd(here::here("results/dominance"))

# get a list of the files that match the variable value of interest
value <- "d_0.2"

# which rep do you want to look at? 
rep <- 1

# find files that match and create figure captions with parameter details
files <- list.files()[grep(value, list.files())]
names <- str_remove(files, ".csv")
param <- str_split(names, "_", simplify = TRUE)[,c(3,5,7,9,11,13,15,18)]
caption <- paste(paste("Generations per year", param[rep, 1]),
                paste("Pop size", param[rep, 2]),
                paste("Seasonal Balance", param[rep, 3]),
                paste("Number of Loci", paste0(param[rep, 4], "\n")),
                paste("Dominance", param[rep, 5]),
                paste("Epistasis", param[rep, 6]),
                paste("Crossover probability", param[rep, 7]),
                paste("Replicate number", param[rep, 8]), sep = "; ")


# read in file and plot results
sim_results <- read.csv(files[rep])

# plot results less than generation 2000
sim_results %>%
  filter(genz < 2000) %>%
  plot_freq(figure_caption = caption, fig_title = value)

# plot results between genz 5000 and 6000
sim_results %>%
  filter(genz > 5000 & genz < 6000) %>%
  plot_freq(figure_caption = caption, fig_title = value)








