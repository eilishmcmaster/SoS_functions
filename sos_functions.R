classifier_function <- function(qdf, plateau){
  qdf$dominant_pop <- colnames(qdf[,(2:(plateau+1))])[apply(qdf[,(2:(plateau+1))],1,which.max)]
  qdf$dom_pop_proportion <- apply(qdf[,(2:(plateau+1))],1,max)
  return(qdf)
}

metadata.read <- function(species){
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
dart.remove.samples <- function(dms, missingness){
  na_per_row <- rowSums(is.na(dms[["gt"]]))/ncol(dms[["gt"]])
  high_missing <- which(na_per_row > missingness)
  sample_names <- names(high_missing)
  names(high_missing) <- NULL
  
  dms$gt <- dms$gt[-high_missing, ]
  dms$sample_names <- dms$sample_names[-high_missing]
  
  if(unique(sample_names %in% rownames(dms$gt))==FALSE){
    print("Samples have been successfully removed from dms$gt")}
  else{stop("Huston, we have a problem (with dms$gt)")}
  
  if(unique(sample_names %in% dms$sample_names)==FALSE){
    print("Samples have been successfully removed from dms$sample_names")  }
  else{stop("Huston, we have a problem (with dms$sample_names)")}
  
  if(identical(rownames(dms[["gt"]]), dms$sample_names)==TRUE){
    print("Samples are in order")}
  else(stop("Huston, we have a problem (these eggs scrambled)"))
  
  paste(length(high_missing), "samples have been removed due missingness >", missingness)
  
  return(dms)
}

remove.by.meta <- function(dms, meta){
  missing <- which(!(rownames(dms$gt) %in% meta$sample))
  dms$gt <- dms$gt[-missing, ]
  dms$sample_names <- dms$sample_names[-missing]
  return(dms)}

remove.by.list<- function(dms, list){
  missing <- which(!(rownames(dms$gt) %in% list))
  cat(length(missing))
  dms$gt <- dms$gt[-missing, ]
  dms$sample_names <- dms$sample_names[-missing]
  
  for(i in  1:4){
    main_meta <- dms$meta[[i]]
    
    dms$meta[[i]] <- main_meta[-missing]
  }
  dms$meta$analyses  <- dms$meta$analyses[-missing,]

  return(dms)}




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
filter <- function(x, min){
  a1_sum <- sum(x== 0 | x== 2, na.rm=TRUE)
  a1_freq <- a1_sum/(length(x)*2)
  a2_sum <- sum(x== 1 | x== 2, na.rm=TRUE)
  a2_freq <- a2_sum/(length(x)*2)
  
  if(a1_freq<(1-min) & a1_freq >min &&
     a2_freq<(1-min) & a2_freq >min
  ){ #safe zone
    return("keep")
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


venner <- function(dms, pops, min_af){
  groups <- unique({{pops}})
  out <- vector()
  
  ds <- dms$gt
  ds <- ds[,which(apply(ds,2,filter, min_af)=="keep")] #remove the low frequency snp sites
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
