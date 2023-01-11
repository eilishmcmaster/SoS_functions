
allele_counter <- function(x, out){ # counts the number of alleles present in a dms$gt 
  genos <- unique(x) # genotypes at this locus
  if(1 %in% genos | 0 %in% genos & 2 %in% genos){ # hetero loci
    out <- 2 #two alleles are present
  }else{
    if(0 %in% genos | 2 %in% genos){
      out <- 1 # one allele is present 
    }else{
      out <- 0 # if there are no alleles present dont add anything 
    }
  }
  return(out)
}

resample_analysis_function <- function(dms,schemes, min_maf, pop){
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  gt <- ds[,which(keepers>=min_maf)]
 # get loci with maf >= min_maf
  # dx <- which(apply(dms$gt,2,filter)<min_maf) # get loci with maf <0.05
  # dms_x <- remove.snps.from.dart.data(dms,dx) # remove the loci with maf<0.05
  # total_alleles <- 2*length(dms_x$locus_names)
  total_alleles <- 2*ncol(gt) # get total number of alleles (two alleles per loci) FOR THIS RESAMPLE GROUP -- not the original population
  cat("Total alleles: ", total_alleles,"\n")
  
  schemes_out <- mat.or.vec(length(schemes), 5) # make empty df
  # for each resampling scheme get the total number of common alleles present 
  for(i in 1:length(schemes)){
    s    <- schemes[[ i ]]
    ivec <- s$svec # locations of the resampled individuals
    gti  <- gt[ ivec, ] # get the altcount of the resampled individuals
    presence <- apply(gti,2,allele_counter) # for each loci, determine if there are one, two, or NA alleles found and write to a df
    schemes_out[i,1] <- s$nseed 
    schemes_out[i,2] <- s$nfam
    schemes_out[i,3] <- sum(presence) #total number of alleles found
    schemes_out[i,4] <- sum(presence)/total_alleles # proportion of all common alleles
    schemes_out[i,5] <- paste(pop)
  }
  colnames(schemes_out) <- c("nseed","nfam", "alleles","prop_total_alleles", "pop")
  return(schemes_out)
  
}

# scheming_function <- function(dmv){ # get the sampling schemes
#   
#   pop    <- as.vector(dmv$meta$site)
#   fam    <- as.vector(dmv$meta$analyses[,"families"])
#   tissue <- dmv$meta$analyses[,"tissue"]
#   tissue[which(dmv$meta$analyses[,"tissue"] == "mother")] <- "M"
#   tissue[which(dmv$meta$analyses[,"tissue"] == "seedling")] <- "P"
#   
#   # set up randomization schemes
#   
#   nr = 10 # number of resamples per scheme
#   scheme_list <- list()
#   scheme_list[[ 1 ]] <- list(nseed=1, nfam=9, npop=NULL, nr=nr)
#   scheme_list[[ 2 ]] <- list(nseed=2, nfam=5, npop=NULL, nr=nr)
#   scheme_list[[ 3 ]] <- list(nseed=5, nfam=2, npop=NULL, nr=nr)
#   scheme_list[[ 4 ]] <- list(nseed=10, nfam=1, npop=NULL, nr=nr)
#   
#   set.seed(9823984)
#   
#   schemes <- list()
#   cs <- 1
#   
#   # loop through the schemes
#   for (i in 1:length(scheme_list)) {
#     
#     s <- scheme_list[[ i ]]
#     # loop through replicates
#     for (i in 1:s$nr) { # schemes is a list of lists where each list is one resample (nseed and nfam specified) and the locations of individuals in the dataset are recorded. This allows the resampled individuals to be found in the proceeding simulations. 
#       svec <- resample_progeny(pop, fam, tissue, nseed=s$nseed, nfam=s$nfam)
#       sout <- list(nseed=s$nseed, nfam=s$nfam, svec=svec)
#       schemes[[ cs ]] <- sout
#       cs <- cs + 1
#     }
#   }
#   return(schemes)
# }


scheming_functionx <- function(dmv, seedvector, famvector, nr){ # get the sampling schemes
  
  pop    <- as.vector(dmv$meta$analyses[,"pop"])
  fam    <- as.vector(dmv$meta$analyses[,"families"])
  tissue <- dmv$meta$analyses[,"tissue"]
  tissue[which(dmv$meta$analyses[,"tissue"] == "mother")] <- "M"
  tissue[which(dmv$meta$analyses[,"tissue"] == "seedling")] <- "P"
  
  # set up randomization schemes
  if(isTRUE(length(seedvector)==length(famvector))){
    print("Seed and family vector are the same length, proceeding")
  }else{
    print("Seed and family vectors are different lengths, terminating process")
    stop()
  }
  
  scheme_list <- list()
  
  for(i in 1:length(seedvector)){
    scheme_list[[i]] <- list(nseed=seedvector[i], nfam=famvector[i], npop=NULL, nr=nr)
  }
  
  set.seed(9823984)
  schemes <- list()
  cs <- 1
  # loop through the schemes
  for (i in 1:length(scheme_list)) {
    
    s <- scheme_list[[ i ]]
    # loop through replicates
    for (i in 1:s$nr) { # schemes is a list of lists where each list is one resample (nseed and nfam specified) and the locations of individuals in the dataset are recorded. This allows the resampled individuals to be found in the proceeding simulations. 
      svec <- resample_progeny(pop, fam, tissue, nseed=s$nseed, nfam=s$nfam)
      sout <- list(nseed=s$nseed, nfam=s$nfam, svec=svec)
      schemes[[ cs ]] <- sout
      cs <- cs + 1
    }
  }
  return(schemes)
}



# jasons function

count_progeny <- function( pop, fam, tissue ) {
  
  families <- unique(fam)
  ifam     <- which(tissue=="M")
  ipro     <- which(tissue=="P")
  npro     <- length(ipro)
  nfam     <- length(families)
  
  out <- data.frame(population=as.character(rep("",nfam)), mother=as.character(rep("",nfam)), progeny=as.numeric(rep(0,nfam)))
  
  class(out$population) <- 'character'
  class(out$mother) <- 'character'
  
  c <- 1
  for (i in ifam) {
    
    n    <- length(which( fam==fam[i] & tissue == "P"))
    p    <- pop[ which( fam==fam[i] & tissue == "M") ]
    f    <- fam[ which( fam==fam[i] & tissue == "M") ]
    
    out$population[c] <- p
    out[c, 2] <- f
    out[c, 3] <- n
    
    c <- c + 1
  }
  return(out)
}

resample_progeny <- function( pop, fam, tissue, nseed, nfam, npop=NULL ) {
  
  cp <- count_progeny(pop, fam, tissue)
  
  class(cp$population) <- 'character'
  class(cp$mother) <- 'character'
  
  pops     <- unique(cp$population)
  
  if (is.null(npop)) {
    
    npop = length(pops)
    
  }
  
  sout <- NULL
  
  spop <- sample(pops)[1:npop]
  
  for (p in spop) {
    
    # eligible fams
    efams <- cp$mother[which( cp$population == p & cp$progeny >= nseed)]
    
    if (length(efams) == nfam) {
      sfam <- efams
    } 
    
    if (length(efams) > nfam) {
      sfam <- sample(efams)[ 1:nfam ]
    }
    
    if (length(efams) < nfam) {
      cat("Warning: could not find enough families of nominated size...\nPopulation", p, efams, nfam); stop();
    }
    
    # for fams in sfam, choose seeds
    
    for (f in sfam) {
      
      # eligible seeds
      eseeds <- which( fam==f & tissue == "P")
      
      if (length(eseeds) == nseed) {
        sseed <- eseeds
      } 
      
      if (length(eseeds) > nseed) {
        sseed <- sample(eseeds)[ 1:nseed ]
      }
      
      if (length(eseeds) < nseed) {
        cat("Warning: could not find enough seeds for nominated family...\n"); stop();
      }
      sout <- c(sout, sseed)
      
    }
    
  }
  
  return(sort(sout))
  
}

