# real-data values
observed_test_stats <- read.csv("piBUSS/data/observed_test_statistics.csv")

# Compute posterior-predictive p-values
random_effects_pps <- lapply(list.files("piBUSS/randomEffects/",full.names=TRUE),read.csv)
random_effects_pps <- do.call(rbind,random_effects_pps)

random_effects_pppvals <- sapply(1:dim(observed_test_stats)[2],function(j){
  (sum(random_effects_pps[,j] < observed_test_stats[1,j]) + 0.5 * sum(random_effects_pps[,j] == observed_test_stats[1,j])) / dim(random_effects_pps)[1]
})

gtr_pps <- lapply(list.files("piBUSS/GTR/",full.names=TRUE),read.csv)
gtr_pps <- do.call(rbind,gtr_pps)

gtr_pppvals <- sapply(1:dim(observed_test_stats)[2],function(j){
  (sum(gtr_pps[,j] < observed_test_stats[1,j]) + 0.5 * sum(gtr_pps[,j] == observed_test_stats[1,j])) / dim(gtr_pps)[1]
})