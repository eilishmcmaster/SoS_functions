


check_relationship <- function(row, external_df){
  if (!is.na(row["relationship"])) {
    return(row["relationship"])
  }else{
    if(isTRUE(row["Var1"] == row["Var2"])) {
      return("same_sample")
    }
    else if (isTRUE(external_df$genotype[external_df$sample %in% row["Var1"]] == 
                    external_df$genotype[external_df$sample %in% row["Var2"]])) {
      # Check if the individuals in the row are siblings
      return("technical_replicate")}
    # Check if the individuals in the row are parent-offspring
    if (isTRUE(external_df$sample %in% row["Var1"]) & 
        isTRUE(external_df$mother %in% row["Var2"])) {
      return("parent-offspring")
    } else if (isTRUE(external_df$families[external_df$sample %in% row["Var1"]] == 
                      external_df$families[external_df$sample %in% row["Var2"]])) {
      # Check if the individuals in the row are siblings
      return("sibling")
    } else if (isTRUE(external_df$site[external_df$sample %in% row["Var1"]] == 
                      external_df$site[external_df$sample %in% row["Var2"]])) {
      # Check if the individuals in the row are from the same site
      return("same_site")
    } else if (isTRUE(external_df$dominant_pop[external_df$sample %in% row["Var1"]] == 
                      external_df$dominant_pop[external_df$sample %in% row["Var2"]])) {
      # Check if the individuals in the row are from the same genetic group
      return("same_genetic_group")
    }
    else{ return("unrelated") }
  }
  
}


dist_kinship_matrix <- function(gt){
  kin <- as.matrix(dist(gt, diag=TRUE))
  kin_invert <- 1- (kin/max(kin))
  return(kin_invert)
}


popkin_kinship_matrix <- function(gt){
  X <- t(gt)
  subpops <- rownames(gt)
  kin <- popkin(X, subpops)
  return(kin)
}

snprelate_kinship_matrix <- function(snpKin){
  snpKin_matrix <- snpKin$kinship
  colnames(snpKin_matrix) <- snpKin$sample.id
  rownames(snpKin_matrix) <- snpKin$sample.id
  return(snpKin_matrix)
}

kinship_by_relationship <- function(kin,parent_offspring, meta){ # kin is any pairwise kinship matrix, meta must have sample column 
  kin[upper.tri(kin, diag=FALSE)] <- NA
  kin_df <- melt(kin, na.rm=TRUE)
  
  if("tissue" %in% colnames(meta)){
    kin_df <- merge(kin_df, parent_offspring, all.x=TRUE) %>% as.data.frame()
  }
  
  kin_df$relationship <- apply(kin_df, 1, check_relationship, external_df = meta)
  return(kin_df)
}
