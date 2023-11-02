read.meta.data.new <- function (dart_data, basedir, species, dataset, 
                             nas = "-") 
{
  metafile <- paste(basedir, species, "/meta/", species, "_", 
                    dataset, "_meta.xlsx", sep = "")
  if (file.exists(metafile)) {
    cat("\n")
    cat(" Reading data file:", metafile, "\n")
  }
  else {
    cat(" Fatal Error: the metadata file", metafile, "does not exist \n")
    stop()
  }
  m <- read_excel(metafile, sheet = 1, col_names = TRUE, na = nas)
  if (any(names(m) == "sample") & any(names(m) == "lat") & 
      any(names(m) == "long")) {
    cat(" Found sample, lat and long columns in metadata \n")
  }
  else {
    cat(" Fatal Error: did not find important sample, lat or long column in metadata \n")
    stop()
  }
  mm <- match(rownames(dart_data$gt), m$sample)
  mi <- intersect(m$sample, rownames(dart_data$gt))
  mm_real <- which(m$sample %in% mi)
  num_dart_samples <- nrow(dart_data$gt)
  num_meta_samples <- nrow(m)
  num_match_samples <- length(mm_real)
  if (num_match_samples > 1) {
    cat(" Found metadata for ", num_meta_samples, " samples \n")
    cat(" This includes overlap with ", num_match_samples, 
        " samples \n")
    cat(" out of ", num_dart_samples, "in DArT genotypes \n")
  }
  else {
    cat(" Fatal Error: no matching sample information between meta-data and DArT genotypes \n")
    stop()
  }
  missing_in_meta <- rownames(dart_data$gt)[is.na(m$sample[mm])]
  meta_ordered <- m[mm_real, ]
  sample_names <- as.character(meta_ordered$sample)
  site <- as.character(meta_ordered$site)
  lat <- as.numeric(meta_ordered$lat)
  long <- as.numeric(meta_ordered$long)
  
  cat(" Adding analysis fields to meta data list \n")
  an <- meta_ordered[, 1:(ncol(meta_ordered))] %>% as.matrix()
  an <- an[,-c("sample","site","lat","long")]
  analyses <- as.matrix(an)
  
  meta_data <- list(sample_names = sample_names, site = site, 
                    lat = lat, long = long, analyses = analyses)
  return(meta_data)
}




classifier_function <- function(qdf, plateau){
  qdf$dominant_pop <- colnames(qdf[,(2:(plateau+1))])[apply(qdf[,(2:(plateau+1))],1,which.max)]
  qdf$dom_pop_proportion <- apply(qdf[,(2:(plateau+1))],1,max)
  return(qdf)
}

metadata.read <- function(species){ # read in metadata exported from rnr database
  meta <- fread(paste(species, "/meta/",species,"_meta.csv", sep=""))
  meta1 <- meta %>% select_if(~!all(is.na(.))) # remove columns of all NA
  colnames(meta1)[colnames(meta1) == 'NSWnumber'] <- 'sample'
  colnames(meta1)[colnames(meta1) == 'decimalLatitude'] <- 'lat'
  colnames(meta1)[colnames(meta1) == 'decimalLongitude'] <- 'long'
  colnames(meta1)[colnames(meta1) == 'locality'] <- 'site'
  ###sample, site, lat, long, and then 2 analysis columns with any name..
  return(meta1)
}


custom.read <- function(species, dataset){
  file <- paste(species, "/meta/",species,"_", dataset, "_meta.xlsx", sep="")
  custom <- read_excel(file,
                       col_names=TRUE, sheet=1)
  return(custom)
}

#test


#used to remove samples with high missingness
remove.by.missingness <-function(dms_o, missingness){
  dms <- dms_o
  na_per_row <- rowSums(is.na(dms[["gt"]]))/ncol(dms[["gt"]])
  high_missing <- which(na_per_row > missingness)
  
  if(length(high_missing)==0){
    print("There are no high missing samples")
    return(dms_o)
  } else{
    sample_names <- names(high_missing)
    names(high_missing) <- NULL
    
    dms$gt <- dms$gt[-high_missing, ]
    dms$sample_names <- dms$sample_names[-high_missing]
    
    if("site" %in% names(dms$meta)){
      dms$meta$site <- dms$meta$site[-high_missing]}
    
    if("lat" %in% names(dms$meta)){
      dms$meta$lat <- dms$meta$lat[-high_missing]}
    
    if("long" %in% names(dms$meta)){
      dms$meta$long <- dms$meta$long[-high_missing]}
    
    if("sample_names" %in% names(dms$meta)){
      dms$meta$sample_names <- dms$meta$sample_names[-high_missing]}
    
    if("analyses" %in% names(dms$meta)){
      dms$meta$analyses <- dms$meta$analyses[-high_missing,]}
    return(dms)
    
    if(isFALSE(unique(sample_names %in% rownames(dms$gt)))){
      print("Samples have been successfully removed from dms$gt")}
    else{stop("Huston, we have a problem (with dms$gt)")}
    
    if(isFALSE(unique(sample_names %in% dms$sample_names))){
      print("Samples have been successfully removed from dms$sample_names")  }
    else{stop("Huston, we have a problem (with dms$sample_names)")}
    
    if(isTRUE(identical(rownames(dms[["gt"]]), dms$sample_names))){
      print("Samples are in order")}
    else(stop("Huston, we have a problem (these eggs scrambled)"))
    
    paste(length(high_missing), "samples have been removed due missingness >", missingness)
    
    return(dms)
  }
}

dart.remove.samples <- remove.by.missingness

remove.by.meta <- function(dms, meta){
  missing <- which(!(rownames(dms$gt) %in% meta$sample))
  dms$gt <- dms$gt[-missing, ]
  dms$sample_names <- dms$sample_names[-missing]
  return(dms)}

remove.by.list <-function(dms, list){ # list of samples to keep
  missing <- which(!(dms$sample_names %in% list))
  if(length(missing)!=0){
    cat(length(missing))
    dms$gt <- dms$gt[-missing, ]
    dms$sample_names <- dms$sample_names[-missing]
    
    meta_names <- which(!names(dms$meta) %in% "analyses") # get the meta that re not called "analyses"
    
    for(i in  meta_names){
      main_meta <- dms$meta[[i]]
      
      dms$meta[[i]] <- main_meta[-missing]
    }
    dms$meta$analyses  <- dms$meta$analyses[-missing,]
  }
  return(dms)}

remove.by.maf <- function(dms, maf){
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  dms$gt <- ds[,which(keepers>=maf)]
  dms$locus_names <- dms$locus_names[which(keepers>=maf)]
  dms$locus_repro <- dms$locus_repro[which(keepers>=maf)]
  dms$locus_pos <- dms$locus_pos[which(keepers>=maf)]
  dms$locus_nuc <- dms$locus_nuc[which(keepers>=maf)]
  
  return(dms)
}


#continuous colour variable PCA
pca_grad_funct <- function(df, x, y, n1,n2, group, gn, colour1, colour2){
  lat_p <- ggplot(df, aes(x={{x}}, y={{y}}, color={{group}}))+geom_point()+theme_bw()+
    scale_colour_gradient(low = paste(colour1), high = paste(colour2), na.value="lightgrey")+
    labs(color=paste(gn), x=pcnames[n1], y=pcnames[n2])+
    theme(legend.position="bottom", legend.text = element_text(angle = 45, vjust = 0, hjust = 0.15))
  return(lat_p)
}

#categorical colour variable PCA
pca_cat_funct <- function(df, x, y,n1,n2, group, gn){
  plot <- ggplot({{df}}, aes(x={{x}}, y={{y}}, colour={{group}}))+
    geom_point()+theme_bw()+
    theme(legend.position="bottom")+
    labs(color=paste(gn),x=pcnames[n1], y=pcnames[n2])
  return(plot)
}

#continuous colour variable PCA
pca_grad_funct2 <- function(df, x, y, n1,n2, group, gn, colour1, colour2){
  lat_p <- ggplot(df, aes(x={{x}}, y={{y}}, color={{group}}))+geom_point()+theme_bw()+
    scale_colour_gradient(low = paste(colour1), high = paste(colour2), na.value="lightgrey")+
    theme(legend.position="bottom", legend.text = element_text(angle = 45, vjust = 0, hjust = 0.15))
  return(lat_p)
}

#categorical colour variable PCA
pca_cat_funct2 <- function(df, x, y,n1,n2, group, gn){
  plot <- ggplot({{df}}, aes(x={{x}}, y={{y}}, colour={{group}}))+
    geom_point()+theme_bw()+
    theme(legend.position="bottom")
  return(plot)
}

named_list_maker <- function(variables, palette, num){
  var <- as.character(unique(na.omit(variables)))
  getPalette2 <- colorRampPalette(brewer.pal(n={num}, paste({palette})))
  colz <- getPalette2(length(var)) #get palette for bargraphs
  names(colz) <- var
  return(colz)
}

# where bs is basic stats and old is the variable to group by (in this case,)
bstat_summary <- function(bs, old){
  Ho <- colMeans(bs$Ho) #	A table with number of populations columns and number of loci rows– of observed heterozygosities
  He <- colMeans(bs$Hs, na.rm=TRUE) # A table –with umber of populations columns and number of loci rows– of observed gene diversities (aka expected heterozygosity)
  Fis <- colMeans(bs$Fis, na.rm=TRUE)
  n1 <- as.data.frame(table(old)) # samples per site
  n <- n1[,2]
  out <- data.frame(Ho,He,Fis)
  table <- cbind(out,n)
  table$Site <- rownames(table)
  rownames(table) <- NULL
  table <- table[,c(5,1,2,3,4)]
  return(table)
}

geo_heat_function <- function(matrix){
  palette <-  colorRamp2(c(0, max(matrix)), c("white", "#80B1D3"))
  geo <- Heatmap(matrix, col=palette,
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 row_names_max_width = unit(15, "cm"),
                 border_gp = gpar(col = "black", lty = 1),
                 name="Distance (km)",
                 cluster_columns = FALSE,
                 cluster_rows=FALSE,
                 row_order=order(rownames(matrix)),
                 column_order=order(colnames(matrix))
  )
  return(geo)
}

gene_heat_function <- function(matrix, font){
  gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))
  gene <- Heatmap(matrix, col=gene_col,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8),
                  row_names_max_width = unit(15, "cm"),
                  border_gp = gpar(col = "black", lty = 1), 
                  cluster_columns = FALSE,
                  cluster_rows=FALSE,
                  row_order=order(rownames(matrix)),
                  column_order=order(colnames(matrix)),
                  name="Pairwise Fst",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", matrix[i, j]), x, y, gp = gpar(fontsize = font))})
}


geo_heat_function2 <- function(matrix, axistext){
  palette <-  colorRamp2(c(0, max(matrix)), c("white", "#80B1D3"))
  geo <- Heatmap(matrix, col=palette,
                 row_names_gp = gpar(fontsize = {{axistext}}),
                 column_names_gp = gpar(fontsize = {{axistext}}),
                 row_names_max_width = unit(15, "cm"),
                 border_gp = gpar(col = "black", lty = 1),
                 name="Distance (km)",
                 cluster_columns = TRUE,
                 cluster_rows=TRUE
  )
  return(geo)
}

gene_heat_function2 <- function(matrix,geo,anno1, anno2, insidetext, axistext){
  gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))
  gene <- Heatmap(matrix, bottom_annotation = anno1, right_annotation = anno2,
                  col=gene_col,
                  row_names_gp = gpar(fontsize = {{axistext}}),
                  column_names_gp = gpar(fontsize = {{axistext}}),
                  row_names_max_width = unit(15, "cm"),
                  border_gp = gpar(col = "black", lty = 1), 
                  # cluster_columns = FALSE,
                  # cluster_rows=FALSE,
                  # row_order=order(rownames(matrix)),
                  column_order=column_order(geo),
                  name="Pairwise Fst",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", matrix[i, j]), x, y, gp = gpar(fontsize = {{insidetext}}))})
}

admix_plotter <- function(data, splitby, textsize, textangle){
  plot <- ggplot(data, aes(x=sample_names, y=Q, fill=population))+
    geom_bar(position="stack", stat="identity")+
    theme_few()+labs(y="Admixture\ncoefficient (Q)", x=element_blank(),
                     fill="Source\npopulation")+scale_fill_manual(values=colz)+
    facet_grid(cols=vars({{splitby}}), scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size=4), strip.text.x = element_text(size = 6))+
    scale_y_continuous(limits = c(0,1.001), expand=c(0,0))+
    theme(strip.text.x = element_text(angle = textangle, size=textsize), panel.spacing = unit(0.07, "lines"))
  
  return(plot)
}


scatterpie_plot_function <- function(df, x, scale_fill, xlim, ylim, r){
  plot <- ggplot(ozmaps::abs_ste) + geom_sf(fill="#f9f9f9", colour="grey") +
    coord_sf(xlim = xlim, ylim = ylim) + labs(y=element_blank(), x=element_blank(), fill="Source\npopulation")+
    geom_scatterpie(aes(x=long.y, y=lat.y, group =x, r = {r}),data =df,
                    cols=colnames(df)[3:(plateau+2)],  alpha=1, size=0.01, colour="black")+#
    theme_few()+theme(axis.text.x = element_text(angle=90) )+ scale_fill_manual(values=scale_fill)
  return(plot)
}

map_plot_function <- function(df, colour, xlim, ylim, legend){
  plot <- ggplot(ozmaps::abs_ste) + geom_sf(fill="#f9f9f9", colour="grey") +
    coord_sf(xlim = xlim, ylim = ylim) + 
    labs(y=element_blank(), x=element_blank(), fill=legend)+
    geom_point(data = df, mapping = aes(x = long, y = lat, fill={{colour}}),
               colour="black",pch=21, size=2)+theme_few()+
    theme(legend.key.size = unit(0, 'lines'), axis.text.x = element_text(angle=90),
          legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,-5,-5,-5))+
    guides(colour = guide_legend(title.position = "top", ncol=1))
  return(plot)
}

map_plot_function2 <- function(df, colour, xlim, ylim, legend){
  plot <- ggplot(ozmaps::abs_ste) + geom_sf(fill="#f9f9f9", colour="grey") +
    coord_sf(xlim = xlim, ylim = ylim) + 
    labs(y=element_blank(), x=element_blank(), fill=legend)+
    geom_point(data = df, mapping = aes(x = long.y, y = lat.y, fill={{colour}}),
               colour="black",pch=21, size=2)+theme_few()+
    theme(legend.key.size = unit(0, 'lines'), axis.text.x = element_text(angle=90),
          legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,-5,-5,-5))+
    guides(fill = guide_legend(title.position = "top", ncol=1))
  return(plot)
}


new.read.dart.xls.onerow <- function (basedir, species, dataset, topskip, nmetavar, nas = "-", 
          altcount = TRUE, euchits = FALSE, misschar = "-", seq2fa = FALSE, 
          fnum = fnum) 
{
  require(readxl)
  datafile <- paste(basedir, species, "/dart_raw/Report-", 
                    dataset, ".xlsx", sep = "")
  if (file.exists(datafile)) {
    cat("\n")
    cat(" Reading data file:", datafile, "\n")
    x <- read_excel(datafile, sheet = 4, skip = topskip, 
                    col_names = TRUE, na = nas)
    if (any(colnames(x) == "AlleleID")) {
    }
    else {
      cat("   Dataset does not include variable AlleleID... Check names... \n")
      if (any(colnames(x) == "CloneID")) {
        cat("   CloneID column found... proceeding with CloneID as names of loci. \n")
        colnames(x)[colnames(x) == "CloneID"] <- "AlleleID"
      }
      else {
        cat("   Locus names cannot be assigned. Error. \n")
        stop()
      }
    }
    if (any(names(x) == "RepAvg")) {
      cat("  includes key variable RepAvg\n\n")
    }
    else {
      cat(" Warning: Dataset does not include variable RepAvg! Check for errors\n\n")
      stop()
    }
    NACloneID <- is.na(x$CloneID)
    if (any(NACloneID)) {
      num_NACloneID <- length(which(NACloneID))
      cat(" Found ", num_NACloneID, " missing CloneID values. Removing from data. \n\n")
      x <- x[-which(NACloneID), ]
    }
    else {
      cat(" No missing CloneID values. \n\n")
    }
    if (!euchits) {
      cat(" Ignoring information on DArT locus alignment to Eucalyptus genome \n\n")
    }
    else {
      cat(" Saving information on DArT locus alignment to Eucalyptus genome \n")
      aln_save_status <- save.eucalypt.genome.hits(x, basedir, 
                                                   species, dataset)
    }
    locus_labels <- as.character(x$CloneID)
    locus_names <- as.character(lapply(strsplit(as.matrix(locus_labels), 
                                                split = "[|]"), function(x) x[1]))
    locus_repro <- as.numeric(x$RepAvg)
    locus_calls <- as.numeric(x$CallRate)
    locus_pos <- as.integer(x$SnpPosition)
    locus_SNP <- as.character(x$SNP)
    locus_nuc <- as.character(lapply(strsplit(as.matrix(locus_SNP), 
                                              split = "[:]"), function(x) x[2]))
    num_snps <- nrow(x)
    num_loci <- length(unique(locus_names))
    num_cols <- ncol(x)
    num_samp <- num_cols - nmetavar
    cat(" Initial data scan -- \n")
    cat("   Samples: ", num_samp, " \n")
    cat("   SNPs: ", num_snps, " \n")
    cat("   Loci: ", num_loci, " \n\n")
    tgt <- x[, (nmetavar + 1):num_cols]
    gt <- as.matrix(t(tgt))
    sample_names <- colnames(x)[(nmetavar + 1):num_cols]
    colnames(gt) <- as.vector(locus_labels)
    cat(" Creating a DArT data list containing:                       \n")
    cat("             Genotypes                    --  $gt            \n")
    cat("             Sample Names                 --  $sample_names  \n")
    cat("             Locus Names                  --  $locus_names   \n")
    cat("             Locus Reproducibility Scores --  $locus_repro   \n")
    cat("             Position of SNP in locus     --  $locus_pos     \n")
    cat("             Data filtering treatments    --  $treatment     \n")
    cat("             Position of SNP in locus     --  $locus_pos     \n")
    cat("             Method of data encoding gt   --  $encoding      \n")
    cat("             Nucleotides in this SNP      --  $locus_nuc     \n\n")
    treatment <- "raw"
    encoding <- "DArT"
    dart_data <- list(gt = gt, sample_names = sample_names, 
                      locus_names = locus_names, locus_repro = locus_repro, 
                      locus_pos = locus_pos, locus_nuc = locus_nuc, encoding = encoding, 
                      treatment = treatment)
    if (altcount) {
      dart_data <- encode.dart2altcount(dart_data)
    }
    else {
      cat(" Warning: genotypes encoded in dart onerow format, 1=hom alt, 2=het \n\n")
    }
    return(dart_data)
  }
  else {
    datafile <- paste(basedir, species, "/dart_raw/Report_", 
                      dataset, "_SNP_mapping_2.csv", sep = "")
    cat("\n")
    cat(" Reading data file:", datafile, "\n")
    x <- read.delim(datafile, sep = ",", na = "-", stringsAsFactors = FALSE, 
                    header = FALSE)
    rows_before_gt <- max(which(x[, 1] == "*"))
    cols_before_gt <- max(which(x[1, ] == "*"))
    colnames(x) <- x[rows_before_gt + 1, ]
    x <- x[-(1:(rows_before_gt + 1)), ]
    if (any(colnames(x) == "AlleleID")) {
    }
    else {
      cat("   Dataset does not include variable AlleleID... Check names... \n")
      if (any(colnames(x) == "CloneID")) {
        cat("   CloneID column found... proceeding with CloneID as names of loci. \n")
        colnames(x)[colnames(x) == "CloneID"] <- "AlleleID"
      }
      else {
        cat("   Locus names cannot be assigned. Error. \n")
        stop()
      }
    }
    if (any(colnames(x) == "RepAvg")) {
    }
    else {
      cat("   Dataset does not include variable RepAvg... Check names... \n")
      stop()
    }
    NAID <- is.na(x$AlleleID)
    if (any(NAID)) {
      num_NAID <- length(which(NAID))
      cat(" Found ", num_NAID, " missing AlleleID values. Removing from data. \n\n")
      x <- x[-which(NAID), ]
    }
    if (!seq2fa) {
    }
    else {
      cat("   Writing sequences to a fasta file \n")
      seq_fname <- write_allele_seq_fasta(x, basedir, species, 
                                          dataset)
    }
    locus_labels <- as.character(x$AlleleID)
    locus_names <- as.character(lapply(strsplit(as.matrix(locus_labels), 
                                                split = "[|]"), function(x) x[1]))
    locus_repro <- as.numeric(x$RepAvg)
    locus_calls <- as.numeric(x$CallRate)
    locus_pos <- as.integer(x$SnpPosition)
    locus_SNP <- as.character(x$SNP)
    locus_nuc <- as.character(lapply(strsplit(as.matrix(locus_SNP), 
                                              split = "[:]"), function(x) x[2]))
    num_snps <- nrow(x)
    num_loci <- length(unique(locus_names))
    num_cols <- ncol(x)
    num_samp <- num_cols - cols_before_gt
    cat(" Initial data scan -- \n")
    cat("   Samples: ", num_samp, " \n")
    cat("   SNPs: ", num_snps, " \n")
    cat("   Loci: ", num_loci, " \n\n")
    tgt <- x[, (cols_before_gt + 1):num_cols]
    gt <- as.matrix(t(tgt))
    class(gt) <- "numeric"
    sample_names <- colnames(x)[(cols_before_gt + 1):num_cols]
    colnames(gt) <- as.vector(locus_labels)
    cat(" Creating a DArT data list containing:                       \n")
    cat("             Genotypes                    --  $gt            \n")
    cat("             Sample Names                 --  $sample_names  \n")
    cat("             Locus Names                  --  $locus_names   \n")
    cat("             Locus Reproducibility Scores --  $locus_repro   \n")
    cat("             Position of SNP in locus     --  $locus_pos     \n")
    cat("             Data filtering treatments    --  $treatment     \n")
    cat("             Position of SNP in locus     --  $locus_pos     \n")
    cat("             Method of data encoding gt   --  $encoding      \n")
    cat("             Nucleotides in this SNP      --  $locus_nuc     \n\n")
    treatment <- "raw"
    encoding <- "DArT"
    dart_data <- list(gt = gt, sample_names = sample_names, 
                      locus_names = locus_names, locus_repro = locus_repro, 
                      locus_pos = locus_pos, locus_nuc = locus_nuc, encoding = encoding, 
                      treatment = treatment)
    if (altcount) {
      dart_data <- encode.dart2altcount(dart_data)
    }
    else {
      cat(" Warning: genotypes encoded in dart onerow format, 1=hom alt, 2=het \n\n")
    }
    return(dart_data)
  }
}



individual_kinship_by_pop <- function(dart_data, basedir, species, dataset, pop, maf=0.1, mis=0.2, as_bigmat=TRUE) {
  require(SNPRelate)
  popvec <- unique(pop)
  kinlist <- list()
  nsamp <- nrow(dart_data$gt)
  igds_file <- dart2gds(dart_data, basedir, species, dataset)
  igds <- snpgdsOpen(igds_file)
  
  for (i in 1:length(popvec)) {
    ipop <- popvec[i]
    isamps <- which(pop == ipop)
    iout <- snpgdsIBDMoM(igds, sample.id=rownames(dart_data$gt)[isamps] , maf=maf, missing.rate=mis, num.thread=1, kinship=TRUE)
    ikout <- iout$kinship
    rownames(ikout) <- rownames(dart_data$gt)[isamps]
    colnames(ikout) <- rownames(dart_data$gt)[isamps]
    kinlist[[ ipop ]] <- ikout
    
  }
  
  snpgdsClose(igds)
  
  if (as_bigmat) {
    bigmat <- matrix( rep(0, nsamp*nsamp ), nrow=nsamp )
    rownames(bigmat) <- rep("", nrow(bigmat))
    colnames(bigmat) <- rep("", nrow(bigmat))
    
    istart <- 1
    for (i in 1:length(kinlist)) {
      
      im <- kinlist[[i]]
      iN <- nrow(im)
      
      if (istart == 1) {
        istop = iN
      } else {
        istop = istop + iN
      }
      
      cat(istart, istop, "\n")
      bigmat[istart:istop, istart:istop] <- im
      rownames(bigmat)[istart:istop] <- rownames(im)
      colnames(bigmat)[istart:istop] <- rownames(im)
      istart = istart + iN
      
    }
    return(bigmat)
  } else {
    return(kinlist)
  }
}

dart2newhy <- function(dart_data, basedir, species, dataset,meta=NULL) {
  
  if (dart_data$encoding == "altcount") {
    
    cat(" Dart data object for ", dataset, "in species", species, "\n")
    cat(" Dart data object found with altcount genotype encoding. Commencing conversion to lfmm. \n")
    nh_gt <- dart_data$gt
    
  } else {
    cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
  }
  
  if (is.null(meta)) {
    cat(" Meta data file not specified \n")
  } else {
    cat("Meta info included in samples file \n")
    meta=meta
  }
  
  nh_gt[ nh_gt == 0 ] <- 11
  nh_gt[ nh_gt == 1 ] <- 12
  nh_gt[ nh_gt == 2 ] <- 22
  
  nh_gt[ is.na(nh_gt) ] <- 0
  
  treatment <- dart_data$treatment 
  
  dir <- paste(basedir, species, "/popgen",sep="")
  if(!dir.exists(dir)) {
    cat("  Directory: ", dir, " does not exist and is being created. \n")
    dir.create(dir)
  } else {
    cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
  }
  
  dir <- paste(basedir, species, "/popgen/",treatment,sep="")
  
  if(!dir.exists(dir)) {
    cat("  Directory: ", dir, " does not exist and is being created. \n")
    dir.create(dir)
  } else {
    cat("  Directory: ", dir, " already exists...  \n")
  }
  
  nh_dir    <- paste(RandRbase,species,"/popgen/",treatment,"/newhy", sep="")
  
  if(!dir.exists(nh_dir)) {
    cat("  NewHybrids directory: ", nh_dir, " does not exist and is being created. \n")
    dir.create(nh_dir)
  } else {
    cat("  NewHybrids directory: ", nh_dir, " already exists, content will be overwritten. \n")
  }
  
  nh_gt_file   <- paste(nh_dir,"/",species,"_",dataset,".txt",sep="")
  nh_H_file    <- paste(nh_dir,"/",species,"_",dataset,"_header.txt",sep="")
  nh_S_file    <- paste(nh_dir,"/",species,"_",dataset,"_samples.txt",sep="")
  nh_L_file    <- paste(nh_dir,"/",species,"_",dataset,"_loci.txt",sep="")
  
  nS <- nrow(nh_gt); vS <- 1:nS; mS <- cbind(vS,rownames(nh_gt), meta)
  nL <- ncol(nh_gt); vL <- paste("L", 1:nL, sep=""); mL <- cbind(vL, colnames(nh_gt))
  
  write.table(mS, file=nh_S_file, sep=",",quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(mL, file=nh_L_file, sep=" ",quote=FALSE, row.names = FALSE, col.names = FALSE)
  
  sink(nh_gt_file)
  cat(c("NumIndivs ", nS, "\n"))
  cat(c("NumLoci ", nL, " \n"))
  cat(c("Digits 1\n"))
  cat(c("Format Lumped \n\n"))
  sink()
  write(c("LocusNames", vL), ncolumns=(nL+1), file=nh_gt_file, sep=" ", append=TRUE)
  sink(nh_gt_file, append = TRUE); cat(c("\n")); sink()
  write.table(cbind(vS, nh_gt), file=nh_gt_file, sep=" ",quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
  
  return(nh_dir)
  
}


common_allele_count <- function(gt, w=NULL, cthresh=1) {
  
  # estimate allele frequencies
  if (is.null(w)) {
    
    alt_counts  <- colSums(gt) 
    ref_counts  <- colSums(2-gt)
    
  } else {
    
    if(nrow(gt) == length(w)) {
      #makes a matrix the same size as gt with 0 and 1 where every column is w vector 
      wm <- matrix(rep(w,ncol(gt)),ncol=ncol(gt),byrow=FALSE)  
      
      alt_counts  <- colSums(gt*wm) #alternative allele counts for samples where w=1
      ref_counts  <- colSums((2-gt)*wm) # main allele counts where w=1
      
    } else {
      cat("   weights supplied: must have length equal to number of rows in gt \n")
    }
  }
  
  min_counts <- alt_counts #minor allele count is alt count 
  for ( i in 1:ncol(gt) ) {
    if ( alt_counts[i] > ref_counts[i]) {
      min_counts[i] <- ref_counts[i] #if ref allele count is less than alt count, alt is minor allele
    } 
    
  }
  
  minor_allele_counts <- min_counts
  number_common_alleles <- length(which( min_counts >= cthresh ))
  
  return( list(number_common_alleles=number_common_alleles, minor_allele_counts=minor_allele_counts) )
}


opt_calculator <- function(gt, group_df,N_t_vec){ # input is hirsuta+abgma gt and cluster df 
  groups <- unique(group_df$"cluster")
  out <- vector()
  cat(groups)
  for(i in 1:length(groups)){
    #get the names of the samples in the current group being analysed
    names <- as.vector(group_df[group_df$"cluster"==groups[i], "sample_names"])
    main_gt2 <- gt[which(rownames(gt) %in% names),] #cut the gt matrix down to only have samples from the group
    max_wts <- rep(2, nrow(main_gt2))

    solutions <- bigger_optimisation(N_t_vec, main_gt2, max_wts)
    out <- c(out, list(solutions))
  }
  names(out) <- groups
  return(out)
}


# opt_calculator <- function(gt, group_df){ # input is hirsuta+abgma gt and cluster df 
#   groups <- unique(group_df$"cluster")
#   out <- vector()
#   for(i in 1:length(groups)){
#     #get the names of the samples in the current group being analysed
#     names <- as.vector(group_df[group_df$"cluster"==groups[i], "sample_names"])
#     main_gt2 <- gt[which(rownames(gt) %in% names),] #cut the gt matrix down to only have samples from the group
#     max_wts <- rep(2, nrow(main_gt2))
#     out_list <- list()
#     if( nrow(main_gt2)<20){
#       N_t_vec <- c(5,10, nrow(main_gt2))
#     }else{
#       N_t_vec <- c(20, nrow(main_gt2))
#     }
#     
#     for ( i in 1:length(N_t_vec) ) {
#       
#       N_t <- N_t_vec[i]
#       cat("\n Running ", N_t, " ...\n")
#       
#       initial_weights = rep(0, nrow(main_gt2))
#       initial_weights[sample(c(1:nrow(main_gt2)))[1:N_t]] <- 1
#       out_list[[ i ]] <- run_optimization_set(main_gt2, N_t, max_wts, initial_weights)
#       
#     }
#     
#     solutions <- out_list[[1]]$nei_opt$weight[num_steps,]
#     for (i in 2:length(N_t_vec)) {
#       solutions <- cbind(solutions, out_list[[i]]$nei_opt$weight[num_steps,])
#     }
#     
#     colnames(solutions) <- as.character(N_t_vec)
#     solutions <- as.data.frame(solutions)
#     solutions$Sample <- rownames(main_gt2)
#     
#     out <- c(out, list(solutions))
#     
#   }
#   names(out) <- groups
#   return(out)
# }

dif_function <- function(x,y){
  present <- sum(x==2 & y==2)
  all <- sum(y==2)
  out <- present/all
  return(out)
}

allele_sum <- function(x, min){
  c_al_sum <- sum(x != 2, na.rm=TRUE)
  prop <- c_al_sum/length(x)
  if(prop<(1-min) & prop >min){
    out <- 2 # has both alleles
  } else{
    out <- 1 # has one allele
  }
  return(out)
}

prop_calculator <- function(major_gt, min1, minor_gt, min2, group_df){
  groups <- unique(group_df$"cluster")
  out <- vector()
  for(i in 1:length(groups)){
    #get the names of the samples in the current group being analysed
    names <- as.vector(group_df[group_df$"cluster"==groups[i], "sample_names"])
    main_gt2 <- major_gt[which(rownames(major_gt) %in% names),] #cut the gt matrix down to only have samples from the group
    minor_gt2 <- minor_gt[which(rownames(minor_gt) %in% names),]
    # proportion of the total minor alleles present in the major population that are present in the minor population
    maj_acount <- apply(main_gt2, 2, allele_sum, min1)
    min_acount <- apply(minor_gt2, 2, allele_sum, min2)
    
    prop <- dif_function(min_acount, maj_acount)
    out <- c(out, prop)
    names(out)[i] <- groups[i]
    
  }
  return(out)
}

bigger_optimisation <- function(N_t_vec, gt, max_wts){
  out_list <- list()
  gt2 <- gt

  if(nrow(gt2)==0){stop("GT is empty")}
  
  else{
  multi <- floor(max(N_t_vec)/nrow(gt2))

  if(multi>10){
    solutions2 <- cat("Requested population is more than 10x the size of source")
  }else{
    if(multi>0){
      for(i in 1:multi){
        gt2 <- rbind(gt2, gt2) # make the df bigger
      }
    }
    
    for ( i in 1:length(N_t_vec) ) {
      N_t <- N_t_vec[i]
      cat("\n Running ", N_t, " ...\n")
      initial_weights = rep(0, nrow(gt2))
      
      initial_weights = rep(0, nrow(gt2))
      initial_weights[sample(c(1:nrow(gt2)))[1:N_t]] <- 1 #randomly makes n_t number of them equal to 1
      out_list[[ i ]] <- run_optimization_set(gt2, N_t, max_wts, initial_weights)
      plot(out_list[[i]]$nei_opt$value, main = paste("max_T=",max_t))
    }
    solutions <- out_list[[1]]$nei_opt$weight[num_steps,]
    for (i in 2:length(N_t_vec)) {
      solutions <- cbind(solutions, out_list[[i]]$nei_opt$weight[num_steps,])
    }
    
    colnames(solutions) <- paste0("n",as.character(N_t_vec))
    solutions <- as.data.frame(solutions)
    solutions$Sample <- rownames(gt)
    solutions2 <- aggregate(.~Sample, data=solutions, sum) 
    }

  
  return(solutions2)
}}

#for the venn diagram alleles
get_minor_allele_frequencies <- function( gt ) {
  
  alleles   <- 2*(colSums(!is.na(gt))) # total number of alleles for a locus (all samples in gt)
  alt_freq  <- colSums(gt,na.rm=TRUE) / (alleles) # get the allele frequencies (altcount data can be summed because 0=homo1, 1=het, 2=homo2)
  ref_freq  <- 1-alt_freq # get the alternative allele frequency
  
  min_freq <- alt_freq
  for ( i in 1:ncol(gt) ) { # assign the minor allele to smaller value
    
    if ( isTRUE(alt_freq[i] > ref_freq[i] )) {
      
      min_freq[i] <- ref_freq[i]
      
    } 
    
  }
  return(min_freq) # return minor allele frequency for each locus in gt 
}

# filter <- function(x){ # find maf BAD 
#   if (length(which(!is.na(x)))==0){ # if everything is NA, MAF= 0 
#     maf <- 0
#   }else{
#     a2_freq <- sum(x, na.rm = TRUE)/(length(which(!is.na(x)))*2) # data is altcount where homo2=0 heter0=1, homo2=2, so allele2 count is sum
#     a1_freq <- 1-a2_freq #allele 1 frequency 
#     if(a1_freq<a2_freq){ #choose the smaller allele frequency
#       maf <- a1_freq
#     }else{
#       maf<- a2_freq
#     }
#   }
#   return(maf)
# }


single_site_genepop_basicstats <- function(dms, min, group){
  # This function makes the genepop file for a dms with a single group (eg one species from one site),
  # and then runs basicstats on that genepop. It removes loci where MAF is  < min specified in the input. 
  # Also, loci should be filtered by missingness using `remove.poor.quality.snps(dms, min_repro=0.96,max_missing=0.3)` before usage. 
  # This method is based on Jasons dart2genepop but is suitable for single group datasets.
  
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  ds <- ds[,which(keepers>=min)]
  cat(paste(ncol(ds))," loci are being used\n")
  
  if (ncol(ds) >=10){
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
  
  # if(ncol(ds) >=50){# make into genepop format
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
  # }else{
  #   print(paste("WARNING: not enough loci for", group))
  #   return(NULL)
  # }
  
}


multispecies_stats <- function(dms, maf, var, missing= NULL, remove_monomorphic=NULL){ # calculates whole species stats for a dms where species =sp
  species <- unique(var)
  species <- species[!is.na(species)]
  print(species)
  out_list <- list()
  
  if(is.null(missing)){
    missing <- 0.3
  }
  
  for(i in 1:length(species)){
    dmsx <- remove.by.list(dms, dms$sample_names[(var %in% paste(species[i]))]) %>% 
      remove.poor.quality.snps(., min_repro=0.96,max_missing=missing) %>%
      remove.by.maf(., maf)
    
    if(isTRUE(remove_monomorphic)){
      dmsx <- remove.fixed.snps(dmsx) 
      print("Fixed SNPs removed")
    }
    
    out <- single_site_genepop_basicstats(dmsx, maf, paste(species[i]))
    out$n <- length(dmsx$sample_names)
    
    out_list[[i]] <- out
  }
  out_df <- do.call(rbind, out_list)
  return(out_df)
}

species_site_stats <- function(dms, maf, pop_var, site_var, missing=NULL, remove_monomorphic=NULL){ 
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
  
  # remove whitespaces
  dms[["meta"]][["analyses"]][,site_var] <- gsub("\\s", "_", dms[["meta"]][["analyses"]][,site_var])
  dms[["meta"]][["analyses"]][,site_var] <- gsub(",", "", dms[["meta"]][["analyses"]][,site_var])
  
  dms[["meta"]][["analyses"]][,pop_var] <- gsub("\\s", "_", dms[["meta"]][["analyses"]][,pop_var])
  dms[["meta"]][["analyses"]][,pop_var] <- gsub(",", "", dms[["meta"]][["analyses"]][,pop_var])
  
  
  out_list <- list()
  genetic_group <- unique(dms[["meta"]][["analyses"]][,pop_var])
  
  if(is.null(missing)){
    missing <- 0.3
  }
  
  for(i in 1:length(genetic_group)){
    dmsx <- remove.by.list(dms, dms[["sample_names"]][(dms[["meta"]][["analyses"]][,pop_var] %in% paste(genetic_group[i]))]) %>% 
      remove.poor.quality.snps(., min_repro=0.96,max_missing=missing) %>%
      remove.by.maf(., maf)
    print(paste("number of samples: ",length(dmsx$sample_names)))
    print(paste("number of loci (after missingness filter): ",length(dmsx$locus_names)))
    
    
    # removes samples with only one sample per site
    tab <- table(dmsx[["meta"]][["analyses"]][,pop_var], dmsx[["meta"]][["analyses"]][,site_var]) %>% as.data.table(.)
    if(1 %in% tab[,N]){
      tab <- tab[N == 1, ]
      print(paste("removing", tab," because n=1" ))
      not_small <- dmsx[["sample_names"]][which(!(dmsx[["meta"]][["analyses"]][,pop_var] %in% tab$V1 &
                                                    dmsx[["meta"]][["analyses"]][,site_var] %in% tab$V2))]
      dmsx <- remove.by.list(dmsx, not_small)
    }
    
    if(isTRUE(remove_monomorphic)){
      dmsx <- remove.fixed.snps(dmsx) 
      print("Fixed SNPs removed")
    }
    
    sites <- dmsx[["meta"]][["analyses"]][,site_var]
    print((unique(sites)))
    
    if(length(unique(sites))>=2){
      out <- multi_site_genepop_basicstats(dmsx, 0, paste(genetic_group[i]), dmsx[["meta"]][["analyses"]][,site_var])
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
      out <- single_site_genepop_basicstats(dmsx, 0, paste(unique(sites)))
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

matcher2 <- function(df2, loci){
  df <- df2[-1]
  out <- vector()
  if(0 %in% df | 2 %in% df){
    out <- append(out, loci[loci==df2[1],2])
  }
  if(1 %in% df | 2 %in% df){
    out <- append(out,loci[loci==df2[1],3])
  }
  return(out)
}


make_allele_list <- function(dms, pops, min_af){
  # this function makes a list of vectors where each vector contains the alleles in a specified population
  # this data can then be used to calculate total alleles, private alleles, and make venn diagrams
  
  groups <- unique({{pops}})
  out <- vector()
  
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  ds <- ds[,which(keepers>=min_af)]
  cat("Found ", ncol(ds), " poly sites\n")     
  loci <- data.frame("loci"=colnames(dms$gt),
                     "allele1"=paste(dms$locus_names,substr(dms$locus_nuc, start=1, stop=1)),
                     "allele2"=paste(dms$locus_names,substr(dms$locus_nuc, start=3, stop=3)))
  for(i in 1:length(groups)){
    df <- ds[{{pops}} %in% groups[i],] # gets the gt frame of just that group
    cat("There are ", nrow(df), " samples in ", groups[i], "\n")
    df2 <- rbind("names"=colnames(df), df)
    alleles <- apply(df2, 2, matcher2, loci)
    names(alleles)<-NULL
    listed <- unlist(alleles)
    out <- c(out, list(listed))
    names(out)[i] <- groups[i]
  }
  return(out)
}

venner <- make_allele_list
# ploidy functions

calculate_private_alleles <- function(populations) {
  # takes the list of vectors produced by venner or make_allele_list 
  # calculates total alleles and private alleles 
  # produces the same results as gg_private_alleles <- poppr::private_alleles(genind_data, level="population", report="table", count.alleles=FALSE)
  
  result_table <- data.frame(population = character(0), private_allele_count = numeric(0), total_allele_count = numeric(0))
  total_alleles <- unique(unlist(populations))
  
  for (i in 1:length(populations)) {
    current_pop <- populations[[i]]
    private_alleles <- unique(current_pop)
    
    for (j in 1:length(populations)) {
      if (i != j) {
        private_alleles <- setdiff(private_alleles, populations[[j]])
      }
    }
    
    total_allele_count <- length(current_pop)
    private_allele_count <- length(private_alleles)
    
    result_table <- rbind(result_table, data.frame(population = names(populations)[i], 
                                                   private_allele_count = private_allele_count, 
                                                   total_allele_count = total_allele_count))
  }
  
  result_table <- rbind(result_table, data.frame(population = "Total", private_allele_count = length(total_alleles),
                                                 total_allele_count = length(total_alleles)))
  
  return(result_table)
}


count_subsetter <- function(dms, count, min){
  ds <- dms$gt
  keepers <- get_minor_allele_frequencies(ds)
  ds <- ds[,which(keepers>=min)]
  cat("Are there any NAs in the altcount data? ", any(is.na(ds)),"\n")
  cat("Loci with NAs:")
  print(table(apply(ds, 2, function(x) any(is.na(x)))))
  
  samples_tk <- dms$sample_names
  loci_tk <- colnames(ds)
  
  s_tk_location <- which(count$sample_names %in% samples_tk)
  l_tk_location <- which(count$locus_labels %in% loci_tk)
  
  count$c1 <- count$c1[l_tk_location,s_tk_location]
  count$c2 <- count$c2[l_tk_location,s_tk_location]
  count$sample_names <- colnames(count$c1)
  count$locus_labels <- count$locus_labels[l_tk_location]
  count$locus_names <- count$locus_names[l_tk_location]
  
  rownames(count$c1) <- count$locus_labels
  rownames(count$c2) <- count$locus_labels
  
  
  count$meta <- count$meta[,s_tk_location]
  count$sample_qc <- count$sample_qc[,s_tk_location]
  
  return(count)
}

mergeLists_internal <- function(o_element, n_element){
  if (is.list(n_element)){
    # Fill in non-existant element with NA elements
    if (length(n_element) != length(o_element)){
      n_unique <- names(n_element)[! names(n_element) %in% names(o_element)]
      if (length(n_unique) > 0){
        for (n in n_unique){
          if (is.matrix(n_element[[n]])){
            o_element[[n]] <- matrix(NA, 
                                     nrow=nrow(n_element[[n]]), 
                                     ncol=ncol(n_element[[n]]))
          }else{
            o_element[[n]] <- rep(NA, 
                                  times=length(n_element[[n]]))
          }
        }
      }
      
      o_unique <- names(o_element)[! names(o_element) %in% names(n_element)]
      if (length(o_unique) > 0){
        for (n in o_unique){
          if (is.matrix(n_element[[n]])){
            n_element[[n]] <- matrix(NA, 
                                     nrow=nrow(o_element[[n]]), 
                                     ncol=ncol(o_element[[n]]))
          }else{
            n_element[[n]] <- rep(NA, 
                                  times=length(o_element[[n]]))
          }
        }
      }
    }  
    
    # Now merge the two lists
    return(mergeLists(o_element, 
                      n_element))
    
  }
  if(length(n_element)>1){
    new_cols <- ifelse(is.matrix(n_element), ncol(n_element), length(n_element))
    old_cols <- ifelse(is.matrix(o_element), ncol(o_element), length(o_element))
    if (new_cols != old_cols)
      stop("Your length doesn't match on the elements,",
           " new element (", new_cols , ") !=",
           " old element (", old_cols , ")")
  }
  return(rbind(o_element, 
               n_element, 
               deparse.level=0))
  return(c(o_element, 
           n_element))
}
mergeLists <- function(old, new){
  if (is.null(old))
    return (new)
  
  m <- mapply(mergeLists_internal, old, new, SIMPLIFY=FALSE)
  return(m)
}


doitall <- function(dms, counts, min, name){
  test2 <- count_subsetter(dms, counts, min)
  
  tr <-  t(test2$c1)
  LMAo <- lapply(split(tr,rownames(tr)), as.list)
  
  tr2 <-  t(test2$c2)
  LMAo2 <- lapply(split(tr2,rownames(tr2)), as.list)
  
  nn <- mergeLists_internal(LMAo, LMAo2)
  
  minor <- lapply(nn, sapply, function(x) min(x)/sum(x))
  a <- do.call(rbind, minor) #make matrix
  major <- lapply(nn, sapply, function(x) max(x)/sum(x))
  b <- do.call(rbind, major) #make matrix
  c <- cbind(a,b)
  
  return(hist(c,
              breaks=50,
              main=paste(name),
              xlab="(mapped reads/ total reads) per allele"))
}
  
#

percent_polymorphic <- function(dms, missingness, maf, var){ # calculates the proportion of loci that are polymorphic vs not 
  species <- unique(var)
  print(species)
  out_df <-  as.data.frame(mat.or.vec(length(species),5))
  colnames(out_df) <- c("species", "all_loci", "poly_loci", "fixed_loci", "ppl")
  
  for(i in 1:length(species)){
    dmsx <- remove.by.list(dms, dms$sample_names[(var %in% paste(species[i]))]) %>% 
      remove.poor.quality.snps(., min_repro=0.96,max_missing=missingness)
    
    dmsx2 <- remove.by.maf(dmsx, maf)
    
    all_loci <- length(dmsx$locus_names)
    poly_loci <- length(dmsx2$locus_names)
    fixed_loci <- all_loci - poly_loci
    ppl <- 100 * (poly_loci/all_loci)
    
    out_df[i, ] <- c(paste0(species[i]), all_loci, poly_loci, fixed_loci, round(ppl, 2))
  }
  return(out_df)
}


remove_missing_loci_by_pop <- function(dms, pop_var, missingness){
  # remove loci with missingness higher than the specified value in all populations provided 
  out_list <- list()
  genetic_group <- unique(dms[["meta"]][["analyses"]][,pop_var])
  
  if(length(genetic_group)<=1){
    stop("ERROR: <=1 population found in pop_var, use a different method")
  }
  
  for(i in 1:length(genetic_group)){
    samples <- dms[["sample_names"]][(dms[["meta"]][["analyses"]][,pop_var] %in% paste(genetic_group[i]))]
    if(length(samples)>=2){ #only for groups with >2 individuals
      dmsx <- remove.by.list(dms, samples) 
      bad_loci <- find.loci.missing.gte.value(dmsx, value=missingness, return_as_names = FALSE) #get loci that fail missingness threshold
      out_list[[i]] <- bad_loci[[1]] # put them in a vector
    }
  }
  bad_indices <- Reduce(intersect, out_list) # get indices that are duplicated i.e. fail the threshold for every group
  dms_filtered <- remove.snps.from.dart.data(dms, bad_indices, input_as_names = FALSE)
  print(paste("dms had", length(dms$locus_names), "loci, now it has", length(dms_filtered$locus_names), "loci"))
  
  return(dms_filtered)
}


remove.sample.by.pop.missingness <-function(dms, pop_var, loci_missingness, sample_missingness){
  # remove samples with missingness higher than the specified value in all populations provided 
  # pop var is found in dms$meta$analyses 
  # loci_missingness is the upper limit of acceptable missingness per locus e.g. 0.3 = all loci with more than 30% NA are removed
  # sample_missingness is the upper limit of acceptable missingness for a sample e.g. 0.3 = all samples with more than 30% NA loci are removed
  
  # example: dms <- remove.sample.by.pop.missingness(dms2, "genetic_group2", loci_missingness = 0.3, sample_missingness = 0.3)

  out_list <- c() #store the good samples
  genetic_group <- unique(dms[["meta"]][["analyses"]][,pop_var]) # get genetic groups
  
  if(length(genetic_group)<=1){
    stop("ERROR: <=1 population found in pop_var, use a different method")
  }
  
  for(i in 1:length(genetic_group)){
    # get the samples in the genetic group
    samples <- dms[["sample_names"]][(dms[["meta"]][["analyses"]][,pop_var] %in% paste(genetic_group[i]))]
    
    if(length(samples)>=2){ #only for groups with >2 individuals
      # make a dms for that genetic group
      dmsx <- remove.by.list(dms, samples) 
      # remove high missingness samples
      dmsx <- remove.poor.quality.snps(dmsx, min_repro=0.96, max_missing=loci_missingness)
      # get missingness per sample
      na_per_row <- rowSums(is.na(dmsx[["gt"]]))/ncol(dmsx[["gt"]])
      # get individuals with missingness less than sample_missingness
      low_missing <- which(na_per_row <= sample_missingness)
      hist(na_per_row, main=genetic_group[i])
      print(paste("samples with low missingness:",length(low_missing), "for genetic group", genetic_group[i] ))
      # add good samples to the list to keep
      good_samples <- names(low_missing) 
      out_list <- append(out_list, good_samples) # put them in a vector
    }
  }
  dms_filtered <- remove.by.list(dms, out_list) # remove the bad samples from the final dms
  print(paste("dms had", length(dms$sample_names), "samples, now it has", length(dms_filtered$sample_names), "sample"))
  
  return(dms_filtered)
}

remove.loci.randomly <- function(dms, number_to_keep){
  # remove loci randomly 
  num_to_remove <- length(dms$locus_names)-number_to_keep
  if(num_to_remove<=0){
    print("Not enough loci to remove")
    stop()
  }
  
  loci_to_remove <- sample(1:length(dms$locus_names), num_to_remove, replace = FALSE) #get loci randomly to remove
  
  dms_filtered <- remove.snps.from.dart.data(dms, loci_to_remove, input_as_names = FALSE)
  
  print(paste("dms had", length(dms$locus_names), "loci, now it has", length(dms_filtered$locus_names), "loci"))
  
  return(dms_filtered)
}


