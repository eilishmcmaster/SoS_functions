# these functions are for making a population and then simulating multiple generations 

# make population ##################################################################

generate_homogeneous_population <- function(n, num_loci, He, Ho){
  # # 2pq = He
  # # p=He/(2q)
  # # if : z= 0.5 * He
  # # p = z/q
  z <- He/2
  
  f2 <- function(q){ # get the function of q that would solve for = 0
    (q - 1) *q + z #=0
  }
  
  # get values of q between 0 and 1 where the function = 0 
  roots <- rootSolve::uniroot.all(f2,c(0,1)) # get q when function = 0 
  
  # x = 0:1
  # df = data.frame(x)
  # ggplot(df, aes(x))+stat_function(fun=f2)+geom_hline(yintercept = 0, color='blue')+
  #   theme_bw()+ylim(-0.5,0.5)+
  #   geom_point(data.frame(x=roots, y=c(0,0)), mapping=aes(x=roots, y=c(0,0)), color='red')
  
  ref_af <- sample(roots, num_loci, replace=TRUE) # get allele frequencies for every locus that result in specified He
  
  geno_matrix <- matrix(nrow=n, ncol=num_loci) # empty matrix 
  
  for(i in 1:num_loci){ # for each locus generate genotypes using allele frequency and desired Ho
    genotypes_1locus <- sample(c(0,1,2), n, replace=TRUE, 
                               prob=c((ref_af[i]-0.5*Ho), Ho, ((1-ref_af[i])-0.5*Ho)))
    geno_matrix[,i] <- genotypes_1locus
  }
  
  
  Ho_original <- sum(geno_matrix==1)/sum(!is.na(geno_matrix))
  He_original <- ref_allele_frequencies(geno_matrix) %>% # get current generation allele frequencies
    as.vector %>% sapply(.,He_from_q) %>% mean() #get He of each locus then average across all loci 
  
  print(paste0("The supplied matrix has ", n, " individuals, ", num_loci, " loci, expected heterozygosity (He) of ",
               round(He_original,3), " and observered heterozygosity (Ho) of ", round(Ho_original, 3)))
  
  return(geno_matrix)
}


# generation simulations ##################################################################

## small functions ##################################################################

### calculate heterozygosity from allele frequency #####
He_from_q <- function(x) {1-(x^2 + (1-x)^2)}

### get genotype when selfing ####
get_self_genotype <- function(x){
  ifelse(x==0, 'AA',
         ifelse(x==2, 'BB',
                sample(c('AA','AB','BB'), size=1, replace=TRUE, prob=c(0.25,0.5,0.25))))
}

# gamete (half genotype) from an individual####
get_gamete <- function(x){ # x is locus in a genotype vector
  ifelse(x==0, 'A',
         ifelse(x==2,'B',
                sample(c("A", "B"), size=1, replace=TRUE, prob=c(0.5, 0.5))))
}

# get reference allele frequency (I treat the 0 homozygous allele as reference)####
ref_allele_frequencies <- function( gt ) {
  
  alleles   <- 2*(colSums(!is.na(gt))) # total number of alleles for a locus (all samples in gt)
  alt_freq  <- colSums(gt,na.rm=TRUE) / (alleles) # get the allele frequencies (altcount data can be summed because 0=homo1, 1=het, 2=homo2)
  ref_freq  <- 1-alt_freq # get the alternative allele frequency
  return(ref_freq) # return minor allele frequency for each locus in gt
}

# run the generation simulations####
# makes a bottlenecked population with n families
# each generation the population stays the same size
# offspring replaces parent 
# when outcrossing, the offspring is a result of a cross within the population 
# when selfing, the offspring is a result of two gametes same parent 
# far out rate is the probability of crossing with the source population 

simulate_populations <- function(pop_size, outcross,far_out_rate, generations, num_loci, He_specified, Ho_specified) {
  He_list <- list()
  out_list <- list()
  name_vector <- character()
  
  for(j in 1:length(pop_size)){
    samples <- generate_homogeneous_population(pop_size[j], num_loci, He=He_specified, Ho=Ho_specified) %>% as.data.frame()
    original_allele_frequencies <- ref_allele_frequencies(samples)
    
    for(z in 1:length(outcross)){
      out <- outcross[z]
      current_generation_genotypes <- samples
      Ho_all <- vector()
      He_all <- vector()
      near_out_rate <- (1-far_out_rate)*outcross[z]
      self_or_out <- matrix(sample(c(0, 1, 2), size = generations * pop_size[j] ,
                                   replace = TRUE, prob = c((1-near_out_rate), near_out_rate, far_out_rate)),
                            nrow = generations)
      # matrix of self or out for every individual in every generation
      
      for(i in 1:generations){ #for each generation
        He <- ref_allele_frequencies(current_generation_genotypes) %>% # get current generation allele frequencies
          as.vector %>% sapply(.,He_from_q) %>% mean() #get He of each locus then average across all loci 
        He_all <- c(He_all, He)
        
        Ho <- sum(current_generation_genotypes==1)/(nrow(current_generation_genotypes)*ncol(current_generation_genotypes))
        Ho_all <- c(Ho_all, Ho) # append to list
        
        # make the next generation individual
        new_genotypes <- matrix(nrow=nrow(current_generation_genotypes), ncol=ncol(current_generation_genotypes))
        
        for(n in 1:pop_size[j]){
          if(self_or_out[i, n]==0){ # if selfing 
            
            new_individual <- sapply(current_generation_genotypes[n,], get_self_genotype) %>% as.vector(.)#alleles from another
          } 
          if(self_or_out[i,n]==1){ #if outcrossing
            alleles1 <- sapply(current_generation_genotypes[n,], get_gamete) %>% as.vector(.)#alleles from self
            alleles2 <- sapply(sample_n(current_generation_genotypes[-n,], 1), get_gamete)%>% as.vector(.) #alleles from another
            new_individual <- paste0(alleles1, alleles2)
          }
          if(self_or_out[i,n]==2){ #if crossing with source population 
            alleles1 <- sapply(current_generation_genotypes[n,], get_gamete) %>% as.vector(.)
            alleles2 <- sapply(original_allele_frequencies, get_outcross_allele) # pick allele from population based on frequency
            new_individual <- paste0(alleles1, alleles2)
          }
          
          new_genotypes[n,] <- ifelse(new_individual == 'AA', 0,
                                      ifelse(new_individual == 'BB', 2, 1))
        }
        current_generation_genotypes <- new_genotypes %>% as.data.frame(.)
      }
      He_list <- append(He_list, list(He_all))
      out_list<- append(out_list, list(Ho_all))
      # name_vector <- c(name_vector, paste0(pop_size[j],"_",outcross[z], "_", far_out_rate[z]))  
      name_vector <- c(name_vector, paste0(pop_size[j],"_",outcross[z]))  
      
    }
  }
  
  result <- list(He = He_list, Ho = out_list, name = name_vector)
  return(result)
}


get_outcross_allele <- function(x){ # x is locus in a genotype vector
  sample(c("A", "B"), size = 1, prob = c(x, 1 - x))
}


# cleanup the output from `simulate_populations`####
get_generation_df <- function(out_list, H) {
  df <- do.call(cbind, out_list[[paste(H)]]) %>% data.frame(.)
  colnames(df) <- out_list[['name']]
  df$generations <- c(1:nrow(df))
  df <- t(apply(df,1, function(x) tapply(x,colnames(df),mean))) %>% as.data.frame(.)
  
  df2 <- df %>% pivot_longer(cols = -generations, names_to = "popsize_outcrossrate", values_to = paste(H))
  df2[,c("popsize","outcross")] <- str_split_fixed(df2$popsize_outcrossrate,"_",2)
  
  return(df2)
}



# Example of how to use ####################################
# 
# devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/simulate_bottleneck.R?raw=TRUE")
# 
# num_loci <- 100
# He_specified <- 0.2
# Ho_specified <- 0.2
# n <- 10
# 
# generations <- 20
# pop_size <- rep(c(5,10),15) # more reps means a smoother outcome -- final data is mean 
# outcross <- c(1,0.5,0) # outcross rates (0=selfing only, 1=non-self only)
# 
# out_list <- simulate_populations(pop_size, outcross, generations, num_loci, He_specified, Ho_specified)
# 
# 
# He_df <- get_generation_df(out_list, 'He')
# Ho_df <- get_generation_df(out_list, 'Ho')
# 
# 
# ggplot(Ho_df, aes(x=generations, y=Ho, groups=factor(popsize_outcrossrate),color=factor(outcross)))+
#   geom_point()+theme_bw()+geom_line()+
#   geom_hline(yintercept = c(0.062, 0.038), linetype="dotted")+
#   # ylim(0)+
#   geom_vline(xintercept = 0, color="red")+
#   facet_wrap(popsize~., drop=FALSE)
# 
# 
# ggplot(He_df, aes(x=generations, y=He, groups=factor(popsize_outcrossrate),color=factor(outcross)))+
#   geom_point()+theme_bw()+geom_line()+
#   geom_hline(yintercept = c(0.062, 0.038), linetype="dotted")+
#   # ylim(0)+
#   geom_vline(xintercept = 0, color="red")+
#   facet_wrap(popsize~., drop=FALSE)
# 
