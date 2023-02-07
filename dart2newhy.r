#' Writes a dart data object to a newhybrids input file
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_lfmm <- dart2lfmm(dart_data, meta=FALSE)
#' }

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

   nh_gt_file   <- paste(nh_dir,"/",species,"_",dataset,"_",analysis,".txt",sep="")
   nh_H_file    <- paste(nh_dir,"/",species,"_",dataset,"_",analysis,".header",sep="")
   nh_S_file    <- paste(nh_dir,"/",species,"_",dataset,"_",analysis,".samples",sep="")
   nh_L_file    <- paste(nh_dir,"/",species,"_",dataset,"_",analysis,".loci",sep="")

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
