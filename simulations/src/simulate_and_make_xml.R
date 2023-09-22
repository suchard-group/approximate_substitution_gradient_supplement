set.seed(42)
source("~/git_repos/approximate_substitution_gradient/xml/COVID/simulations/src/helper.R")

######
# Preliminaries
######

# Number of runs going into the posterior we're basing our simulations on
nrun <- 1

bridge_exponents <- c(1/8,1/4,1/2,1)

# Trees are logged every 50k, parameters every 500 (for ESS purposes)
log_thin_factor <- 50000/500

# Each tree file has 1001 trees (1000 post-0)
# If we remove the first 50% and take every 20th tree, this will leave us with 100 simulations total
keep <- seq(502,1000,20)

# If setting this to another directory, make sure it ends in a /
parent_dir <- ""

pibuss_script_refx <- scan(paste0(parent_dir,"piBUSS/src/piBUSS_randomEffects.xml"),what=character(),sep="\n",strip.white=FALSE)
mcmc_template <- scan(paste0(parent_dir,"simulations/src/mcmc_template.xml"),what=character(),sep="\n",strip.white=FALSE)
target_dir <- "SC2_HMC"
xml_dir <- paste0(parent_dir,"simulations/xml")

dir.create(xml_dir)

param_logs <- vector("list",nrun)
tree_logs <- vector("list",nrun)
for (i in 1:nrun) {
  param_logs[[i]] <- read.table(paste0(parent_dir,"output/",target_dir,"/job_",i,"/SC2_HMC.log"),header=TRUE,stringsAsFactors=FALSE)
  tree_logs[[i]] <- read.nexus(paste0(parent_dir,"output/",target_dir,"/job_",i,"/SC2_HMC.trees"))
  
  # Trim off state 0 for easier handling
  param_logs[[i]] <- param_logs[[i]][-1,]
  tree_logs[[i]] <- tree_logs[[i]][-1]
  
  nstates <- dim(param_logs[[i]])[1]
  if (nstates/log_thin_factor != length(tree_logs[[i]])) {
    stop("Tree or parameter log file dimensions do not match")
  }
  
  if (length(tree_logs[[i]]) != 1000) {
    stop("Unexpected number of trees.")
  }
  
  # Make logs match trees
  param_logs[[i]] <- param_logs[[i]][seq(log_thin_factor,nstates,log_thin_factor),]
  
  param_logs[[i]] <- param_logs[[i]][keep,]
  tree_logs[[i]] <- tree_logs[[i]][keep]
}

param_log <- do.call(rbind,param_logs)
trees <- do.call(c,tree_logs)

nsim <- length(trees)

refx_log <- param_log[,grepl("glmRandCoefficients",names(param_log))]
# We make the following rough classification of the inferred random-effects
shrunk <- c(1, 2, 4, 8, 9, 10, 11, 12)
# these have opposite signs, 3 is negative, 7 is positive
moderate <- c(3, 7)
strong <- c(5, 6)

# Fit Normals to the posteriors and use those to simulate from
# We use MAD because fat Bridge tails can make an SD-based fit very poorly representative in some regimes, and the answers should be equivalent in the other regimes
shrunk_mean <- 0
shrunk_sd <- mad(unlist(refx_log[,shrunk])) 

moderate_abs_mean <- mean(c(-refx_log[,moderate[1]],refx_log[,moderate[2]]))
moderate_sd <- mad(c(-refx_log[,moderate[1]],refx_log[,moderate[2]]))

strong_mean <- mean(unlist(refx_log[,strong]))
strong_sd <- mad(unlist(refx_log[,strong]))

normal_params <- list(strong_mean=strong_mean,
                      strong_sd=strong_sd,
                      moderate_abs_mean=moderate_abs_mean,
                      moderate_sd=moderate_sd,
                      shrunk_mean=shrunk_mean,
                      shrunk_sd=shrunk_sd)

refx_mean_vec <- rep(0,12)
refx_mean_vec[shrunk] <- shrunk_mean
refx_mean_vec[moderate] <- moderate_abs_mean * c(-1,1)
refx_mean_vec[strong] <- strong_mean

refx_sd_vec <- rep(0,12)
refx_sd_vec[shrunk] <- shrunk_sd
refx_sd_vec[moderate] <- moderate_sd
refx_sd_vec[strong] <- strong_sd

simRefx <- function() {
  rnorm(12, refx_mean_vec, refx_sd_vec)
}

######
# Simulate and make XML
######
simulating_values <- matrix(nrow=nsim,ncol=15)
colnames(simulating_values) <- c(paste0("glmRandCoefficients",1:12),"log.kappa","alpha","clock.rate")

set.seed(42)
for (i in 1:nsim) {
  suppressWarnings(rm(this_script))
  
  this_script <- pibuss_script_refx

  # Tree
  this_script <- gsub("NEWICKTREE",write.tree(trees[[i]]),this_script)
  
  # ASRV
  alpha <- param_log$alpha[i]
  this_script <- gsub("GAMMAALPHA",alpha,this_script)
  
  # Clock
  clock.rate <- param_log$clock.location[i]
  this_script <- gsub("CLOCKRATE",clock.rate,this_script)
  
  # Q matrix
  log_kappa <- param_log$log.kappa[i]
  this_script <- gsub("LOGKAPPA",log_kappa,this_script)
  refx <- simRefx()
  this_script <- gsub("LOGRANDOMEFFECTS",paste0(refx,collapse=" "),this_script)
  
  # Site count
  nsites <- 29903

  this_script <- gsub("FROMPARTITION1",1,this_script)
  this_script <- gsub("TOPARTITION1",nsites,this_script)

  tmpfile <- "tmp_pibuss.xml"
  beastjar <- "~/git_repos/beast-mcmc/build/dist/beast.jar"
  beaglepath <- "/usr/local/lib"
  tmpaln <- paste0(dirname(tmpfile),"/sequences.fasta")
  
  # DO THE EVOLUTION, BABY!
  cat(this_script,file=tmpfile,sep="\n")
  system(paste0("java  -Djava.library.path=",beaglepath," -jar ",beastjar," -working -overwrite ",tmpfile))
  
  sim_aln <- read.phyDat(tmpaln,"fasta")

  char_aln <- as.character(sim_aln)
  beast_aln <- formatAlnForXML(char_aln)
  beast_aln <- paste0(beast_aln,collapse="")
  
  for (bridge_exponent in bridge_exponents) {
    suppressWarnings(rm(this_mcmc))
    this_mcmc <- mcmc_template
    
    this_name <- paste0("SC2_simulation_",i,"_bridge_exponent_",bridge_exponent)
    
    this_mcmc <- gsub("ALIGNMENT",beast_aln,this_mcmc)
    this_mcmc <- gsub("BAYESIANBRIDGEEXPONENT",bridge_exponent,this_mcmc)
    this_mcmc <- gsub("NEWICKTREE",write.tree(trees[[i]]),this_mcmc)
    this_mcmc <- gsub("CLOCKRATE",clock.rate,this_mcmc)
    this_mcmc <- gsub("OUTPUTFILENAME",this_name,this_mcmc)
    
    cat(this_mcmc,file=paste0(xml_dir,"/",this_name,".xml"),sep="\n")  
  }
  
  # Record simulating values
  simulating_values[i,] <- c(refx,log_kappa,alpha,clock.rate)
}

write.csv(simulating_values,file=paste0(dirname(xml_dir),"/simulating_values.csv"),quote=FALSE,row.names=FALSE)
