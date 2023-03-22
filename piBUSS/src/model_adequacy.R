library(phangorn)

numAlnTestStats <- function(num.aln) {
  # recover()

  # raw counts of nucleotides at sites
  n <- apply(num.aln,2,function(x){
    sapply(1:4,function(i){
      sum(x == i)
    })
  })
  is_inv <- apply(n,2,function(x){sum(x > 0) == 1})
  
  # Calculate averages, variances, and covariances only on variable sites
  n <- n[,!is_inv]
  p <- apply(n,2,function(x){x/sum(x)})
  m <- cov(t(p))
  
  res_mean <- sapply(1:4,function(i){
    sum(n[i,])
  })
  res_mean <- res_mean/sum(res_mean)

  i <- c(1:4,1,1,1,2,2,3)
  j <- c(1:4,2,3,4,3,4,4)
  res_vcv <- sapply(1:10,function(idx){m[i[idx],j[idx]]})
  res <- c(res_mean,
           res_vcv)
  names(res) <- c("E(A)","E(C)","E(G)","E(T)","Var(A)","Var(C)","Var(G)","Var(T)","Cov(AC)","Cov(AG)","Cov(AT)","Cov(CG)","Cov(CT)","Cov(GT)")
  return(res)
}

maskNumAln <- function(num.aln,mask) {
  for (i in 1:dim(mask)[1]) {
    num.aln[mask[i,1],mask[i,2]] <- -1
  }
  return(num.aln)
}

# This is NOT sorted the way that arr.ind returns results!

aln2mask <- function(aln) {
  char_aln <- as.character(aln)
  per_taxon_mask <- lapply(1:dim(char_aln)[1],function(i){
    is_amb <- !(tolower(aln[i,]) %in% c("a","c","g","t"))
    if ( any(is_amb) ) {
      return(cbind(i,which(is_amb)))
    }
  })
  return(do.call(rbind,per_taxon_mask))
}

asNumAln <- function(aln) {
  char_aln <- as.character(aln)
  apply(char_aln,2,function(x){
    y <- rep(-1,length(x))
    y[x == "a"] <- 1
    y[x == "c"] <- 2
    y[x == "g"] <- 3
    y[x == "t"] <- 4
    return(y)
  })
}