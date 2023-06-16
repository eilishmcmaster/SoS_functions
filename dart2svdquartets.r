#' Returns a dart data object in appropriate files for analysis 
#' with the Peter and Slatkin method of detecting a range expansion 
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gi <- dart2svdq(dart_data, meta.csv)
#' }

dart2svdquartets <- function(dart_data, basedir, species, dataset, add_pop=FALSE, pop=NULL) {

   require(ape)

   # Step 1, get the genotypes ready
   treatment <- dart_data$treatment 
   if (dart_data$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   # genotypes -- change char for missing
   gt <- dart_data$gt
   gt[ is.na(gt) ] <- "?"

   svdq_gt     <- matrix(rep("", nrow(gt)*ncol(gt) ), nrow=nrow(gt))

   for (i in 1:ncol(gt) ) {

      get_het <- function(ref, alt) {
         if ( ref == "A" & alt == "T" | ref == "T" & alt == "A") { het = "W" }
         if ( ref == "C" & alt == "G" | ref == "G" & alt == "C") { het = "S" }
         if ( ref == "A" & alt == "C" | ref == "C" & alt == "A") { het = "M" }
         if ( ref == "G" & alt == "T" | ref == "T" & alt == "G") { het = "K" }       
         if ( ref == "A" & alt == "G" | ref == "G" & alt == "A") { het = "R" }
         if ( ref == "C" & alt == "T" | ref == "T" & alt == "C") { het = "Y" }       

         return(het)

      }

      ref_allele <- substr(dart_data$locus_nuc[i],1,1)
      alt_allele <- substr(dart_data$locus_nuc[i],3,3)
      het_allele <- get_het(ref_allele, alt_allele)

      for ( j in 1:nrow(gt) ) {

           g <- gt[j,i]
   
           if (g == 0)   { svdq_gt[ j,i ] <- ref_allele }
           if (g == 1)   { svdq_gt[ j,i ] <- het_allele }
           if (g == 2)   { svdq_gt[ j,i ] <- alt_allele }
           if (g == "?") { svdq_gt[ j,i ] <- "N"  }

      }

   }

   rownames(svdq_gt) <- rownames(gt)


   # make directory, write files 
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

   svdq_dir    <- paste(basedir,species,"/popgen/",treatment,"/svdq", sep="")
   
   if(!dir.exists(svdq_dir)) {
      cat("  RE directory: ", svdq_dir, " does not exist and is being created. \n")
      dir.create(svdq_dir)
   } else {
      cat("  RE directory: ", svdq_dir, " already exists, content will be overwritten. \n")
   }

   svdq_nexus_file   <- paste(svdq_dir,"/",species,"_",dataset,".nex",sep="")

   if(!add_pop) {

      svdq_gt_out <- svdq_gt

   } else {

     poptext <- "begin sets; \n      taxpartition Populations = \n "

     #svdq_gt_out        <- matrix(rep("", nrow(gt)*ncol(gt) ), nrow=nrow(gt))
     sample_order       <- NULL
     up <- unique(pop)

     c <- 1
     for (i in 1:length(up)) {

         p  <- up[i]
         is <- which(pop == p)
         cs <- c
         
         for (s in is) {
            #svdq_gt_out[c,] <- svdq_gt[s,]
            sample_order <- c(sample_order, s)
            c <- c + 1
         }
         ce <- c-1
         poptext <- paste0(poptext, "    ", p, " : " , cs, "-", ce )
         if (i < length(up)) { poptext <- paste0(poptext,", \n") }
         if (i == length(up)) {poptext <- paste0(poptext,"; \nend;") }
     }
      cat(sample_order)
      svdq_gt_out <- svdq_gt[sample_order,]
      #rownames(svdq_gt_out) <- rownames(gt)[sample_order,]
   } 

   
   write.nexus.data(svdq_gt_out, svdq_nexus_file, interleaved = FALSE, missing="N")
   cat(poptext)

   write(poptext,file=svdq_nexus_file,append=TRUE)
   #cat("   nexus file written... changing headers for beautii... \n")
   #snf      <- scan(svdq_nexus_file, what="character", sep="\n")
   #idt      <- grep("DATATYPE", snf)
   #snf[idt] <- gsub("DNA", "STANDARD", snf[idt]) 
   #snf[idt] <- gsub("INTERLEAVE", "SYMBOLS=\"012\" INTERLEAVE", snf[idt]) 

   #cat("   overwriting nexus file with adjusted FORMAT line... \n")
   #write(snf, svdq_nexus_file)

   return(svdq_nexus_file)

}
