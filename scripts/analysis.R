# results analysis -------------------------------------------------------------

require(tidyverse)
require(rlist)

source(here::here("functions.R"))

# set working directory to the results folder of interest
setwd(here::here("results/generations"))

# get a list of the files that match the variable value of interest
value <- "G_20"

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

# plot results between genz 500 and 600
sim_results %>%
  filter(genz > 900 & genz < 1000) %>%
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

setwd(here::here("results/generations"))
values = c("G_10","G_20","G_50") #this needs to be of length 3 and of increasing order
do_analysis(values) #first boxplot = means, second boxplot = standard deviation, categories appear
                    #in the order you put it in the values variable above
                    #also gives you 3 ANOVA: first means, second standard deviation, third ratio of fixed loci

#first run function below, then call it
do_analysis <- function(values){
  
  files = c()
  for (value in values){
    file <- list.files()[grep(value, list.files())]
    files <- list.append(files, file) #contains all the files regarding your parameter in the correct order
  }
  
  sim_results <- vector(mode = "list", length = 15)
  i = 1
  for (f in files){
    sim_result <- read.csv(f)
    sim_results[i] <- list(sim_result) #contains all dataframes of your parameter (length = 15)
    i <- i+1
  }
  
  maximums <- c()
  for (rep in sim_results){
    maximums <- list.append(maximums, max(rep$genz)) #this ensures that the minimum of all the maximum generations is used for analysis
  }
    
  means <- c()
  standard_devs <- c()
  for (rep in sim_results){
    temp <- rep[rep$genz == min(maximums),] #check that all have same max generation!
    means <- list.append(means, mean(temp$freqs))
    standard_devs <- list.append(standard_devs, sd(temp$freqs))
  }
  
  means1 <- means[1:5] #means of loci frequency when parameter = first input value
  means2 <- means[6:10] #means of loci frequency when parameter = second input value
  means3 <- means[11:15] #means of loci frequency when parameter = third input value
  
  x <- data.frame("Categories" = c("Means1","Means1","Means1","Means1","Means1","Means2","Means2","Means2","Means2","Means2","Means3","Means3","Means3","Means3","Means3"), "Means" = c(means1,means2,means3))
  
  ANOVAofmeans <- aov(x$Means ~ x$Categories)
  print(summary(ANOVAofmeans))
  
  boxplot(x$Means ~ x$Categories,
          ylab="Mean loci frequency",
          xlab= "Dominance values", #REPLACE THIS
          names=c("d = 0.2","d = 0.5","d = 0.8"), #REPLACE THIS
          col= c("gray", "blue", "green")
  )
  
  sd1 <- standard_devs[1:5] 
  sd2 <- standard_devs[6:10] 
  sd3 <- standard_devs[11:15]
  
  y <- data.frame("Categories" = c("SD1","SD1","SD1","SD1","SD1","SD2","SD2","SD2","SD2","SD2","SD3","SD3","SD3","SD3","SD3"), "SDs" = c(sd1,sd2,sd3))
  
  ANOVAofsds <- aov(y$SDs ~ y$Categories)
  print(summary(ANOVAofsds))
  
  boxplot(y$SDs ~ y$Categories,
          ylab="Standard deviation of loci frequencies",
          xlab= "Dominance values",
          names=c("d = 0.2","d = 0.5","d = 0.8"),
          col= c("gray", "blue", "green")
  ) 
  
  fixed <- c()
  for (rep in sim_results){
    temp <- rep[rep$genz == min(maximums),] #check that all have same max!
    temp <- round(temp$freqs, digits = 1) #round to 1 digit so that 0.99 would be considered fixed
    fixed <- list.append(fixed, length(which(temp==0 | temp == 1))/length(temp))
  }
  
  fix1 <- fixed[1:5]
  fix2 <- fixed[6:10]
  fix3 <- fixed[11:15]
  
  z <- data.frame("Categories" = c("fix1","fix1","fix1","fix1","fix1","fix2","fix2","fix2","fix2","fix2","fix3","fix3","fix3","fix3","fix3"), "Fixs" = c(fix1,fix2,fix3))
  
  ANOVAoffixs <- aov(z$Fixs ~ z$Categories)
  print(summary(ANOVAoffixs))
  
  boxplot(z$Fixs ~ z$Categories,
          ylab="Proportion of fixed loci",
          xlab= "Dominance values",
          names=c("d = 0.2","d = 0.5","d = 0.8"),
          col= c("gray", "blue", "green")
  ) 
}
