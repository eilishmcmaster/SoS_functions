species_site_stats <- function(dms, maf, pop_var, site_var){ 
  # This function allows you to calculate site stats for multiple genetic groups at the same time
  # dms has all of the samples youre interested in 
  # MAF is the threshold (0.05 usually)
  # pop_var is the genetic group -- use genetic groups to avoid biasing calculations
  # site_var is the sites within genetic groups to calculate stats for -- has to be a column name in the dms$meta$analyses dataframe 
  
  # example : x <- species_site_stats(dms, maf=0.05, pop_var="sp", site_var="site")
  # @dms Genind structure built by RRtools with the samples of interest
  # @maf Minor allele frequency
  # @pop_var Name of the variable containing populations/genetic neighbourhoods
  # @site_var Name of the variable containing sites

  
  if(isTRUE(site_var %in% names(dms[["meta"]])) &
     isFALSE(site_var %in% colnames(dms[["meta"]][["analyses"]]))){ #if "site" is in dms$meta and isnt in $analyses...
    dms[["meta"]][["analyses"]] <- cbind(dms[["meta"]][["analyses"]], site =unlist(dms[["meta"]][site_var], use.names=FALSE) )
    colnames(dms[["meta"]][["analyses"]])[ncol(dms[["meta"]][["analyses"]])] <- paste(site_var)
  } 
  
  if(isFALSE(site_var %in% colnames(dms[["meta"]][["analyses"]]))){
    stop("ERROR: cannot find site variable")
  }
  if(isFALSE(pop_var %in% colnames(dms[["meta"]][["analyses"]]))){
    stop("ERROR: cannot find population variable")
  }

  if(all(is.na(dms[["meta"]][["analyses"]][,site_var]))){
    stop("ERROR: site_var is empty")
  }
  if(all(is.na(dms[["meta"]][["analyses"]][,pop_var]))){
    stop("ERROR: pop_var is empty")
  }
  
    # removes samples with no site or sp classification
  samples_with_sp_and_site <- dms[["sample_names"]][-which(is.na(dms[["meta"]][["analyses"]][,site_var]) |
                                                              is.na(dms[["meta"]][["analyses"]][,pop_var]))]
  if(length(samples_with_sp_and_site)>=2){
    dms <- remove.by.list(dms, samples_with_sp_and_site)
  }
  if(length(dms[["sample_names"]])<=2){
    stop("ERROR: not enough samples to proceed (<=2)")
  }

  # removes samples with only one sample per site
  tab <- table(dms[["meta"]][["analyses"]][,pop_var], dms[["meta"]][["analyses"]][,site_var]) %>% as.data.table(.) 
  if(1 %in% tab[,N]){
    tab <- tab[N == 1, ]
    not_small <- dms[["sample_names"]][which(!(dms[["meta"]][["analyses"]][,pop_var] %in% tab$V1 &
                                                 dms[["meta"]][["analyses"]][,site_var] %in% tab$V2))]
    dms <- remove.by.list(dms, not_small)
  }

  # remove whitespaces
  dms[["meta"]][["analyses"]][,site_var] <- gsub("\\s", "_", dms[["meta"]][["analyses"]][,site_var])
  dms[["meta"]][["analyses"]][,site_var] <- gsub(",", "", dms[["meta"]][["analyses"]][,site_var])

  dms[["meta"]][["analyses"]][,pop_var] <- gsub("\\s", "_", dms[["meta"]][["analyses"]][,pop_var])
  dms[["meta"]][["analyses"]][,pop_var] <- gsub(",", "", dms[["meta"]][["analyses"]][,pop_var])
  
  
  out_list <- list()
  genetic_group <- unique(dms[["meta"]][["analyses"]][,pop_var])
  
  for(i in 1:length(genetic_group)){
    dmsx <- remove.by.list(dms, dms[["sample_names"]][(dms[["meta"]][["analyses"]][,pop_var] %in% paste(genetic_group[i]))]) %>% 
      remove.poor.quality.snps(., min_repro=0.96,max_missing=0.3) %>%
      remove.by.maf(., maf)
    sites <- dmsx[["meta"]][["analyses"]][,site_var]
    print((unique(sites)))
    
    if(length(unique(sites))>=2){
      out <- multi_site_genepop_basicstats(dmsx, maf, paste(genetic_group[i]), dmsx[["meta"]][["analyses"]][,site_var])
      if(isFALSE(is.null(out))){
      freq <- table(dmsx[["meta"]][["analyses"]][,site_var]) %>% data.frame()
      colnames(freq)[2]<- "n"
      out$genetic_group <- paste(genetic_group[i])
      out <- merge(out, freq, by.x=0, by.y="Var1", all.y=FALSE)
      colnames(out)[1]<- "site"
      out_list[[i]] <- out
      }
    }
    if(length(unique(sites))==1){
      out <- single_site_genepop_basicstats(dmsx, maf, paste(unique(sites)))
      if(isFALSE(is.null(out))){
        out$genetic_group <- paste(genetic_group[i])
        out$site <- rownames(out)
        rownames(out) <- NULL
        out$n <- length(dmsx[["sample_names"]])
        out_list[[i]] <- out
      }
    }
  }
  if(length(out_list)==1){
    out_df <- out_list[[1]]
    return(out_df)
  }
  if(length(out_list)>=2){
    out_df <- do.call(rbind, out_list)
    return(out_df)
  }
  if(length(out_list)==0){
    print("WARNING: no data created")
  }
}


multi_site_genepop_basicstats <- function(dms, min, group, grouping){
  # This function makes the genepop file for a dms with multiple groups with low differentiation (eg one species from multiple sites).
  # After creating the genepop file, this function runs basicstats. It removes loci where MAF is  < min specified in the input for the
  # entire set rather than per group. 
  # It is not suitable for multi-species datasets -- these should be split into single species before use. 
  # Also, loci should be filtered by missingness using `remove.poor.quality.snps(dms, min_repro=0.96,max_missing=0.3)` before usage. 
  # This method is based on Jasons dart2genepop but is suitable for single group datasets.
  
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  ds <- ds[,which(keepers>=min)] # filter it by the MAF and min MAF specified
  cat(paste(ncol(ds))," loci are being used\n") # print the final ammount of loci being used
  
  if(ncol(ds) >=50){# make into genepop format
    old <- c("0","1","2", NA) 
    new <- c("0101","0102","0202","0000")
    ds[ds %in% old] <- new[match(ds, old, nomatch = 0000)] 
    
    # populations
    pops <- unique(grouping) # get population names
    
    # write genepop file
    gf <- paste0(species, "/popgen/genepop_",group,".gen") # make genepop file path
    
    cat(paste0("genepop file: ",species, " with MAF ", paste0(min)), # first line of genepop file
        file=gf,sep="\n")
    cat(colnames(ds),file=gf,sep="\n", append=TRUE) # one loci name per line
    
    remove <- c() # vector for populations excluded from the analysis 
    
    for (i in 1:length(pops)){ #loop for making the population groups
      if (length(which(grouping %in% pops[i]))<=1){ # find if the population is n=1
        cat("Removing population ", pops[i], " due to n=1")
        remove <- c(remove, pops[i]) # add the pop name to remove vector
      }else{
        cat("pop",file=gf,sep="\n", append=TRUE) # add the data to the genepop file
        df <- ds[which(grouping %in% pops[i]),]
        for (j in 1:nrow(df)){
          cat(c(paste0(pops[i],","),df[j,], "\n"),file=gf,sep="\t", append=TRUE)
        }
      }
      
    } #end of pops loop
    
    # # do basic stats
    bs <- diveRsity::basicStats(infile = gf, outfile = NULL,
                                fis_ci = FALSE, ar_ci = TRUE,
                                ar_boots = 1000,
                                rarefaction = FALSE, ar_alpha = 0.05)
    # return(bs)
    # extract the stats
    npop <- length(pops)-length(remove)
    result <- as.data.frame(mat.or.vec(npop,11))
    measurement_names <- rownames(bs$main_tab[[1]])
    population_names  <- names(bs$main_tab) #ls() rearranges the names
    rownames(result) <- pops[!pops %in% remove]
    colnames(result) <- measurement_names
    
    for (r in 1:npop) {
      popstats <- bs$main_tab[[r]][,"overall"] ##extract from a list
      result[r,] <- popstats}
    result$loci <- ncol(ds)
    return(result)
  }else{
    print(paste("WARNING: not enough loci for", group))
    return(NULL)
  }
  
}

single_site_genepop_basicstats <- function(dms, min, group){
  # This function makes the genepop file for a dms with a single group (eg one species from one site),
  # and then runs basicstats on that genepop. It removes loci where MAF is  < min specified in the input. 
  # Also, loci should be filtered by missingness using `remove.poor.quality.snps(dms, min_repro=0.96,max_missing=0.3)` before usage. 
  # This method is based on Jasons dart2genepop but is suitable for single group datasets.
  
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  ds <- ds[,which(keepers>=min)]
  cat(paste(ncol(ds))," loci are being used\n")
  
  if (ncol(ds) >=50){
    # make into genepop format
    old <- c("0","1","2", NA)
    new <- c("0101","0102","0202","0000")
    ds[ds %in% old] <- new[match(ds, old, nomatch = 0000)]
    
    # write genepop file
    gf <- paste0(species, "/popgen/genepop_",group,".gen")
    
    cat(paste0("genepop file: ",species, " with MAF ", paste0(min)),
        file=gf,sep="\n")
    cat(colnames(ds),file=gf,sep="\n", append=TRUE)
    cat("pop",file=gf,sep="\n", append=TRUE)
    for (i in 1:nrow(ds)){
      cat(c("pop1,", ds[i,], "\n"),file=gf,sep="\t", append=TRUE)
    }
    cat("pop",file=gf,sep="\n", append=TRUE)
    for (i in 1:2){
      cat(c("pop2,", ds[i,], "\n"),file=gf,sep="\t", append=TRUE)
    }
    
    # do basic stats
    bs <- diveRsity::basicStats(infile = gf, outfile = NULL,
                                fis_ci = FALSE, ar_ci = TRUE, 
                                ar_boots = 1000, 
                                rarefaction = FALSE, ar_alpha = 0.05)
    # return(bs)
    # extract the stats
    npop <- 1
    result <- as.data.frame(mat.or.vec(npop,11))
    measurement_names <- rownames(bs$main_tab[[1]])
    population_names  <- names(bs$main_tab) #ls() rearranges the names
    rownames(result) <- {{group}}
    colnames(result) <- measurement_names
    
    for (r in 1:npop) {
      popstats <- bs$main_tab[[r]][,"overall"] ##extract from a list
      result[r,] <- popstats}
    
    result$loci <- ncol(ds)
    
    return(result)
  } else{
    print(paste("WARNING: not enough loci for", group))
    return(NULL)
  }
}



# gf <- paste0("/Users/eilishmcmaster/Documents/AcacGord/AcacGord/popgen/genepop_gordonii.gen")
# 
# bs <- diveRsity::basicStats(infile = gf, outfile = NULL,
#                             fis_ci = FALSE, ar_ci = TRUE,
#                             ar_boots = 1000,
#                             rarefaction = FALSE, ar_alpha = 0.05)
# # return(bs)
# # extract the stats
# npop <- length(pops)-length(remove)
# result <- as.data.frame(mat.or.vec(npop,11))
# measurement_names <- rownames(bs$main_tab[[1]])
# population_names  <- names(bs$main_tab) #ls() rearranges the names
# rownames(result) <- pops[!pops %in% remove]
# colnames(result) <- measurement_names



x <- species_site_stats(dms, maf=0.05, pop_var="sp", site_var="site3")
