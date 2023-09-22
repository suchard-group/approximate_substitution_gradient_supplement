library(phangorn)
set.seed(42)
source("piBUSS/src/model_adequacy.R")

######
# Preliminaries
######

nrun <- 1

# Trees are logged every 50k, parameters every 500 (for ESS purposes)
log_thin_factor <- 50000/500

# Each tree file has 1001 trees (1000 post-0)
# If we remove the first 50% and take every 2nd tree, this will leave us with 1000 simulations total
keep <- seq(502,1000,2)

parent_dir <- ""

pibuss_script_refx <- scan("piBUSS/src/piBUSS_randomEffects.xml",what=character(),sep="\n",strip.white=FALSE)

pibuss_script_gtr <- scan("piBUSS/src/piBUSS_GTR.xml",what=character(),sep="\n",strip.white=FALSE)


######
# Summary stats from real alignment
######
real_aln <- read.phyDat("piBUSS/data/aln.fasta","fasta")
real_aln_taxa <- names(real_aln)
real_aln_mask <- aln2mask(real_aln)

aln_taxa <- names(real_aln)
aln_as_num <- asNumAln(real_aln)
obs_test_stats <- numAlnTestStats(aln_as_num)

write.csv(rbind(obs_test_stats),"piBUSS/data/observed_test_statistics.csv",quote=FALSE,row.names=FALSE)

for (model in c("GTR","REfx")) {
  
  if ( model == "GTR" ) {
    target_dir <- "SC2_GTR"
    dir.create("piBUSS/GTR")
    log_name <- "SC2_GTR.log"
    tree_name <- "SC2_GTR.trees"
  } else {
    target_dir <- "SC2_HMC"
    dir.create("piBUSS/randomEffects")
    log_name <- "SC2_HMC.log"
    tree_name <- "SC2_HMC.trees"
  } 
  
  param_logs <- vector("list",nrun)
  tree_logs <- vector("list",nrun)
  for (i in 1:nrun) {
    param_logs[[i]] <- read.table(paste0(parent_dir,"output/",target_dir,"/job_",i,"/",log_name),header=TRUE,stringsAsFactors=FALSE)
    tree_logs[[i]] <- read.nexus(paste0(parent_dir,"output/",target_dir,"/job_",i,"/",tree_name))
    
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
  
  ######
  # Simulate
  ######
  
  for (i in 1:nsim) {
    suppressWarnings(rm(this_script))
    
    if ( model == "GTR" ) {
      this_script <- pibuss_script_gtr
    } else {
      this_script <- pibuss_script_refx
    }
    
    # Tree
    this_script <- gsub("NEWICKTREE",write.tree(trees[[i]]),this_script)
    
    # ASRV
    alpha <- param_log$alpha[i]
    this_script <- gsub("GAMMAALPHA",alpha,this_script)
    
    # Clock
    clock.rate <- param_log$clock.location[i]
    this_script <- gsub("CLOCKRATE",clock.rate,this_script)
    
    # Q matrix
    if ( model == "GTR" ) {
      gtr_params <- param_log[i,grepl("gtr.rates",names(param_log))]
      this_script <- gsub("GTRPARAMS",paste0(gtr_params,collapse=" "),this_script)
    } else {
      log_kappa <- param_log$log.kappa[i]
      this_script <- gsub("LOGKAPPA",log_kappa,this_script)
      refx <- as.numeric(param_log[i,grepl("glmRandCoefficients",names(param_log))])
      this_script <- gsub("LOGRANDOMEFFECTS",paste0(refx,collapse=" "),this_script)
    }
    
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
    sim_aln_taxa <- names(sim_aln)
    
    key <- match(real_aln_taxa,sim_aln_taxa)
    sim_aln <- sim_aln[key]
    num_sim_aln <- asNumAln(sim_aln)
    num_sim_aln <- maskNumAln(num_sim_aln,real_aln_mask)
    
    tmp_test_stats <- numAlnTestStats(num_sim_aln)
    
    if (model == "GTR") {
      statfile <- paste0("piBUSS/GTR/test_statistics_",i,".csv")
    } else {
      statfile <- paste0("piBUSS/randomEffects/test_statistics_",i,".csv")
    }
    write.csv(rbind(tmp_test_stats),statfile,quote=FALSE,row.names=FALSE)
    
  }
}
