##' For repeats that have read pairs flanking them, the implied expansion and length of that repeat is determined. 
##'
##' Summary stats for the implied repeat exansions/lengths are also reported. 
##' @title calculateIRLs_norm
##' @param flanking_file _flanking read pair file
##' @param imd_dist .insert_len_dist file (i.e. the inner mate distance distribution for the sample)
##' @param median_imd the median inner mate distance of your sample; can be determined empirically if imd_dist is supplied
##' @param upper_limit the value used in mapATR % processLastAlignedBase [default= 0.95]
##' @param margin the value used in mapATR % processLastAlignedBase [default= 37.5]
##' @param min_coverage minimum number of flanking read paired that map to a given repeat required for IRL calculation [default= 10]
##' @param output_dir path to output directory
##' @return _implied_repeat_lengths file
##' @author Adam Struck - Intern
calculateIRLs_norm <- function(flanking_file, imd_dist, median_imd=NULL, upper_limit=0.95, margin=37.5, min_coverage=10, output_dir) {
  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }

  flanking <- read.table(flanking_file, header=TRUE, sep="\t")
  names(flanking) <- c("repeat.id", "inner.mate.distance", "n")

  imd.dist <- read.table(imd_dist, header=TRUE, sep="\t")
  names(imd.dist) <- c("inner.mate.distance", "n")
  dist <- rep(imd.dist$inner.mate.distance, imd.dist$n)
  if( upper_limit <= 1 ) {
    upper.limit <- quantile(dist, prob=upper_limit)
  } else {
    upper.limit <- upper_limit
  }
  
  irl.file.header <- c("repeat.id", "ref.repeat.length", "implied.repeat.length", "median.irl", "mean.irl", "std.err.irl", "n")
  fh <- openOutputFileHandle(flanking_file, suffix="implied_repeat_lengths_norm", header=irl.file.header, overwrite=TRUE)
    
  repeats <- unique(flanking$repeat.id)
  if(length(repeats) > 0) {
    for(i.rep in 1:length(repeats))  {
      if ((i.rep %% 50 == 0) == TRUE) { message("Working on repeat: ", i.rep , "/", length(repeats), " ", Sys.time()) }
      repeat.id <- repeats[i.rep]
      rep <- unlist(strsplit(as.character(repeat.id), "[:-]"))
      ref.repeat.length <- as.numeric(rep[3]) - as.numeric(rep[2]) + 1

### modify the imd_dist to reflect the possible imds flanking the repeat ###
      min.imd <- ref.repeat.length - 2 * margin
      max.imd <- min.imd + 2 * upper.limit
      limits  <- (imd.dist$inner.mate.distance >= min.imd) & (imd.dist$inner.mate.distance <= max.imd)
    
      rep.imd.dist <- imd.dist[limits,]
      rep.dist     <- rep(rep.imd.dist$inner.mate.distance, rep.imd.dist$n)
      median.imd   <- median(rep.dist)
    
      implied.repeat.lengths <- getIRLs_norm(flanking, repeat.id, ref.repeat.length, median.imd, min_coverage)
    
      if(!is.null(implied.repeat.lengths)) {

        IRLs        <- paste(sort(implied.repeat.lengths), collapse=",")
        mean.irl    <- mean(implied.repeat.lengths)
        median.irl  <- median(implied.repeat.lengths)
        n           <- length(implied.repeat.lengths)
        std.err.irl <- sd(implied.repeat.lengths)/sqrt(n)
      
        results <- data.frame(repeat.id, ref.repeat.length, IRLs, median.irl, mean.irl, std.err.irl, n)
        writeResults(results, fh)
      }
    }
  }
  close(fh)
  
  file_split <- unlist(strsplit(flanking_file, "/"))
  filename <- file_split[length(file_split)]
  base_split <- unlist(strsplit(filename, "[.]"))
  base_name <- base_split[1]
  output_file <- sprintf("%s.%s", base_name, "implied_repeat_lengths_norm")
  
  sample.irls <- read.table(output_file, header=TRUE, sep="\t")
  sample.irls$mean.irl <- round(sample.irls$mean.irl, 2)
  sample.irls$std.err.irl <- round(sample.irls$std.err.irl, 2)
  sample.irls$corrected.irl <- sample.irls$mean.irl + round(mean(sample.irls$ref.repeat.length - sample.irls$mean.irl, na.rm=TRUE), 2)

  write.table(sample.irls,
              file=output_file,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=FALSE)
  
}

getIRLs_norm <- function(flanking, repeat.id, ref.repeat.length, median.imd, min_coverage) {
  repeat.subset <- flanking[flanking$repeat.id == repeat.id,]
  inner.mate.distances <- rep(repeat.subset$inner.mate.distance, repeat.subset$n)
  if(length(inner.mate.distances) >= min_coverage) {
    implied.repeat.lengths <- (median.imd - inner.mate.distances) + ref.repeat.length
    return(implied.repeat.lengths)
  } else {
    implied.repeat.lengths <- NULL
    return(implied.repeat.lengths)
  }
}
