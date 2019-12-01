# results analysis -------------------------------------------------------------

require(tidyverse)
require(rlist)

source(here::here("functions.R"))

# set working directory to the results folder of interest
setwd(here::here("results/generations"))

# get a list of the files that match the variable value of interest
value <- "G_50"

# which rep do you want to look at? 
rep <- 3

# find files that match and create figure captions with parameter details
files <- list.files()[grep(value, list.files())]
names <- str_remove(files, ".csv")
param <- str_split(names, "_", simplify = TRUE)[,c(3,5,7,9,11,13,15,18)]
caption <- paste(paste("Generations per year", param[rep, 1]),
                paste("Pop size", param[rep, 2]),
                paste("Seasonal Balance", param[rep, 3]),
                paste("Number of Loci", param[rep, 4]),
                paste("\nDominance", param[rep, 5]),
                paste("Epistasis", param[rep, 6]),
                paste("Crossover probability", param[rep, 7]),
                paste("Replicate number", param[rep, 8]), sep = "; ")


# read in file and plot results
sim_results <- read.csv(files[rep])

# plot results less than generation 2000
sim_results %>%
  filter(genz < 2000) %>%
  plot_freq(figure_caption = caption, fig_title = value)

# plot results less than generation 2000 and divisible by 10 to get only the 
# freq at the seasonal turns for gen/year = 20
sim_results %>%
  filter(genz < 2000 & genz %% 10 == 0) %>%
  plot_freq(figure_caption = caption, fig_title = value)

sim_results %>%
  filter(loci %in% as.character(sample(unique(sim_results$loci), 10, replace = FALSE))) %>%
  filter(genz > 2000 & genz < 2100 & genz %% 10 == 0) %>%
  plot_freq(figure_caption = caption, fig_title = value)

# plot results between genz 500 and 600
sim_results %>%
  filter(genz > 2800 & genz < 3000) %>%
  plot_freq(figure_caption = caption, fig_title = value)


# plot the proportion of loci at fixation
sim_results %>%
  filter(genz < 6000) %>%
  group_by(genz) %>%
  filter(freqs == 0 | freqs == 1) %>%
  count() %>%
  mutate(freqs = n / as.numeric(param[rep, 4])) %>%
  ggplot(mapping = aes(x = genz, y = freqs)) +
  geom_line() + 
  ylim(0, 1) +
  ylab("Propotion of fixed loci") +
  xlab("Generations")

# histogram at generation 6000 (if possible)
sim_results <- sim_results[sim_results$genz == max(sim_results$genz),]
hist(sim_results$freqs)
mean(sim_results$freqs)
sd(sim_results$freqs)

# COMPARE different parameters ----------------------------------------------------------------

setwd(here::here("results/cross_over"))
values = c("c_0_", "c_0.001_", "c_0.01","c_0.05", "c_0.1","d_0.5_y_1_c_0.5_") #this needs to be of length 3 and of increasing order
do_analysis(values, "cross_over") #first boxplot = means, second boxplot = standard deviation, categories appear
                    #in the order you put it in the values variable above
                    #also gives you 3 ANOVA: first means, second standard deviation, third ratio of fixed loci




# histogram and density plots --------------------------------------------------

# set working directory to the results folder of interest
setwd(here::here("results/generations"))

# load all files and filter only the final generation
all_tbl <- list.files() %>% 
  map_df(~read_plus(.)) 

tbl_max_gen <- all_tbl %>%
  split(.$filename) %>%
  map(function(.data){filter(.data, genz == max(genz))}) %>%
  do.call("rbind", .)

# plot the density distribution
ggplot(tbl_max_gen, aes(x = freqs, fill = filename)) +
  geom_density(alpha = 0.2) + 
  scale_fill_manual(values = c(rep("#e25dbb", 5), 
                               rep("#6cdf4c", 5),
                               rep("#258bd2", 5)),
                    breaks = unique(all_tbl$filename)[c(1,6,11)],
                    # legend name
                    name = "Generations\nper year",
                    # values tested
                    labels = c("10", "20", "50")) +
  ylab("Density") + 
  xlab("Allele Frequency") +
  theme_classic()

