# Columns RepAvg and CloneID no longer in file as DArT changed format so removed checks from script -Stephanie 20220718

#' Import DArT (onerow) genotype data and stores them in a list 
#'
#' It reads the genotype calls in dart format (encoding=dart) and
#' can covert these to counts of the alternate allele (encoding=altcount)
#'
#' @param basedir  -- base directory for analyses [required]
#' @param misschar -- missing data character [default "-"]
#' @param mincount -- minimum count of an allele [default "-"]
#' @return A list containing counts of loci formatted for use with 
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)

read_dart_counts_csv_faster <- function(datafile, misschar="-", minAlleleCount, minGenotypeCount) {

   cat("\n")
   cat(" Reading data file:", datafile,"\n")

   ### Open and parse dart datafile
   #x <- read.table(datafile,sep=",",na="-", stringsAsFactors=FALSE)
   x <- read.delim(datafile,sep=",",na=c("-", ""), stringsAsFactors=FALSE, header=FALSE)
     
   rows_before_gt <- max(which(x[,1]=="*"))
   cols_before_gt <- max(which(x[1,]=="*"))

   file_meta_info <- x[1:(rows_before_gt+1),(cols_before_gt+1):ncol(x)] 
   colnames(file_meta_info) <- NULL
   rownames(file_meta_info) <- NULL

   colnames(x) <- x[rows_before_gt+1,]
   x <- x[-(1:(rows_before_gt+1)),]

   ### Error checks
   if(any(colnames(x) == "AlleleID")) {
   
   } else {
      
      cat("   Dataset does not include variable AlleleID... Check names... \n");

#      if(any(colnames(x) == "CloneID")) {
#         cat("   CloneID column found... proceeding with CloneID as names of loci. \n");
#         colnames(x)[colnames(x)=="CloneID"] <- "AlleleID"
#      } else {
#         cat("   Locus names cannot be assigned. Error. \n"); stop()
#      }
   }
#   if(any(colnames(x) == "RepAvg")) {
#     
#   } else {
#     cat("   Dataset does not include variable RepAvg... Check names... \n"); stop()
#   }

   NAID <- is.na(x$AlleleID)
   if (any(NAID) ) {
     num_NAID <- length( which( NAID ) )
     cat(" Found ", num_NAID," missing AlleleID values. Removing from data. \n\n")
     x <- x[ -which( NAID ), ]   
   } 
 
  
   ind_onerow <- seq(2,nrow(x),2)
 

 
   # read some important data columns
   locus_labels <- as.character(x$AlleleID[ind_onerow] )
   locus_names  <- as.character(lapply(strsplit(as.matrix(locus_labels),split="[|]"), function(x) x[1]))
   # locus_repro  <- as.numeric(x$RepAvg[ind_onerow])
   locus_calls  <- as.numeric(x$CallRate[ind_onerow])
   locus_pos    <- as.integer(x$SnpPosition[ind_onerow])
   locus_SNP    <- as.character(x$SNP[ind_onerow])
   locus_nuc    <- as.character(lapply(strsplit(as.matrix(locus_SNP),split="[:]"), function(x) x[2]))

   num_cols <- ncol(x)
   num_samp <- num_cols - cols_before_gt
   xg <- x[,(cols_before_gt+1):num_cols]

   x1 <- mat.or.vec(nrow(xg)/2, ncol(xg))
   colnames(x1) <- colnames(xg)
   rownames(x1) <- locus_labels

   x2 <- mat.or.vec(nrow(xg)/2, ncol(xg))
   colnames(x2) <- colnames(xg)
   rownames(x2) <- locus_labels

   ind_row_x2 <- seq(2,nrow(xg),2)
   ind_row_x1 <- ind_row_x2 - 1

   p_x1 = xg[ind_row_x1, ]
   p_x2 = xg[ind_row_x2, ]

   x1 <- as.matrix(p_x1)
   x2 <- as.matrix(p_x2)

   ixg <- ( as.numeric(as.matrix(p_x1)) + as.numeric(as.matrix(p_x2)) < minGenotypeCount)
   ixa <- (p_x1 < minAlleleCount | p_x2 < minAlleleCount)
   x1[ ixg ] <- NA 
   x1[ ixa  ] <- NA 
   x2[ ixg ] <- NA 
   x2[ ixa  ] <- NA 

   num_snps <- nrow(x1)
   num_loci <- length(unique(locus_names))
 
   # print some stats
   cat(" Initial data scan -- \n")
   cat("   Samples: ",num_samp," \n")
   cat("   SNPs: ", num_snps," \n")
   cat("   Loci: ", num_loci," \n\n")

   # place samples in rows, loci in columns

   c1 <- as.matrix(x1)
   class(c1) <- "numeric"

   c2 <- as.matrix(x2)
   class(c2) <- "numeric"
 
   xgn <- as.matrix(xg)
   class(xgn) <- "numeric"

   sample_qc <- rbind( colSums(xgn)/(nrow(xgn)/2), colMeans(c1+c2,na.rm=TRUE), colSums(is.na(c1))/nrow(c1) )
   rownames(sample_qc) <- c("meanDepthCov", "meanDepthCovHetSites", "fractionNA")

   count_data <- list(c1=c1, c2=c2, locus_names=locus_names, sample_names=colnames(c1),
                      locus_labels=locus_labels,
                      meta=file_meta_info, sample_qc=sample_qc)

   return(count_data)
}


