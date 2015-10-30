##' Sorts read pairs into two categories: facing, and fully repetitive.
##'
##' Description of output files:
##' facing:            one read flanks or is anchored outside the repeat (minimum = margin) while its pair aligned inside the repeat (or didn't meet minimum anchor requirements)
##' fully_repetitive:  both reads in pair aligned inside the repeat (or didn't meet minimum anchor requirements)
##' @title processLastAlignedBase
##' @param con database connection
##' @param table_name name of table in database
##' @param input_file last_aligned_base.tsv file
##' @param imd_dist .insert_len_dist file (i.e. the inner mate distance distribution for the sample)
##' @param upper_limit values between 0 and 1 will be used to set prob in quantile() and then used to calculate the maximum allowed distance from the repeat. Values greater than 1 will be applied as the maximum. Facing reads whose ends are farther from the repeat than the value set by this parameter will be not be used. [default= 0.95]
##' @param read_length read length [default= 75]
##' @param output_dir output directory to write results [default= "."]
##' @param overwrite overwrite output files if they exist [default= FALSE]
##' @param append_chr appends "chr" to chromosome names of inner_mate_ranges file [default= TRUE]
##' @param max_gap maximum distance between the last aligned base and the start or end of the repeat [default = 500]
##' @param margin maximum bp overlap between each read and a given repeat repeat [default= read_length/2]
##' @param min_repeat_length minimum repeat length in reference genome [default= read_length/2] 
##' @return _facing, and _fully_repetitive tsv files
##' @author Adam Struck - Intern
##' @export
processLastAlignedBase <- function(con, table_name, input_file, imd_dist=NULL, upper_limit=0.95, margin=(read_length/2), read_length=75, min_repeat_length=(read_length/2), output_dir=".", overwrite=FALSE, append_chr=TRUE) {

  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }
  
  chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

  message("Reading in ", input_file, " at ", Sys.time())

  if(upper_limit <= 1) {
    if( !is.null(imd_dist) ) {
      imd.dist <- read.table(imd_dist, header=TRUE, sep="\t")
      names(imd.dist) <- c("inner.mate.distance", "n")
      dist <- rep(imd.dist$inner.mate.distance, imd.dist$n)
      limits  <- quantile(dist, probs=upper_limit)
      upper.limit <- limits[[1]]
    } else {
      stop("You need to supply an inner mate distance distribution file for your sample")
    }
  } else {
    upper.limit <- upper_limit
  }

  passed <- 0
  for(i.chr in 1:length(chrs))  {
    chr <- chrs[i.chr]
    message("Working on chr: ", i.chr , "/", length(chrs), " ", Sys.time())
    
    if(append_chr == TRUE) {
      sample   <- try(read.table(pipe(sprintf("sed 's/^/chr/g' %s | awk '/%s[\t]/ {print $0}'", input_file, chr))))
      if (class(sample) == "try-error") { next }
    } else {
      sample   <- try(read.table(pipe(sprintf("awk '/%s[\t]/ {print $0}' %s", chr, input_file))))
      if (class(sample) == "try-error") { next }
    }
    colnames(sample) <- c("chr", "pos", "strand", "n")

    passed <- passed + 1
    facing.file.header     <- c("repeat.id", "side", "distance", "n")
    fullrepeat.file.header <- c("repeat.id", "n")
    if(passed == 1) {
      overwrite.val <- overwrite
      fh.facing     <- openOutputFileHandle(input_file, suffix="facing", header=facing.file.header, overwrite=overwrite.val)
      fh.fullrepeat <- openOutputFileHandle(input_file, suffix="fullrepeat", header=fullrepeat.file.header, overwrite=overwrite.val)
    } else {
      overwrite.val <- FALSE
      fh.facing     <- openOutputFileHandle(input_file, suffix="facing", header=facing.file.header, overwrite=overwrite.val)
      fh.fullrepeat <- openOutputFileHandle(input_file, suffix="fullrepeat", header=fullrepeat.file.header, overwrite=overwrite.val)
    }
    
    repeats    <- dbGetQuery(con, sprintf("select * from %s where chr='%s'", table_name, chr))
    repeats    <- repeats[(repeats$end-repeats$start) > min_repeat_length,]
    repeat.ids <- paste(repeats$chr, paste(repeats$start, repeats$end, sep="-"), sep=":")
    
    message("Read in repeats on ", chr, " from table ", table_name, " at ", Sys.time())
    message(paste(capture.output(gc()), collapse="\n"))

    for(i.rep in 1:length(repeat.ids))  {
      repeat_id <- repeat.ids[i.rep]
      if ((i.rep %% 50 == 0) == TRUE) {
        message("Working on repeat: ", i.rep , "/", length(repeat.ids), " ", Sys.time())
      }
      facing.results <- getFacingLAB(sample, repeats, repeat_id, margin, upper.limit)
      if(class(facing.results) == "data.frame") {
        writeResults(facing.results, fh.facing)}
      
      fullrepeat.results <- getFullRepeatsLAB(sample, repeats, repeat_id, margin)
      if(class(fullrepeat.results) == "data.frame") {
        writeResults(fullrepeat.results, fh.fullrepeat)}
    }
    close(fh.facing)
    close(fh.fullrepeat)
  }
}

getFacingLAB<- function(sample, repeats, repeat_id, margin, upper.limit) {
  start <- repeats[repeats$repeat_id == repeat_id, 4]
  end <- repeats[repeats$repeat_id == repeat_id, 5]

  facing.from.right <- sample$pos >= (end - margin)
  right.prox <- (sample$pos - (end - margin)) <= upper.limit
  minus.strand <- sample$strand == "-"
  
  facing.from.left <- sample$pos <= (start + margin)
  left.prox <- ((start + margin) - sample$pos) <= upper.limit
  plus.strand <- sample$strand == "+"
  
  from.left.facing <- sample[(facing.from.left & left.prox & plus.strand),]
  if(dim(from.left.facing)[1] > 0) {
    from.left.facing.dist  <- from.left.facing$pos - start 
    from.left.facing.dist  <- rep(from.left.facing.dist, from.left.facing$n) 
    left.n                 <- table(from.left.facing.dist)
    from.left.dist.to.rep  <- as.integer(names(left.n))
    from.left.results      <- data.frame(repeat_id, side="left", distance_to_repeat=from.left.dist.to.rep, n=as.vector(left.n), row.names=NULL)
  } else { from.left.results <- NULL }
    
  from.right.facing <- sample[(facing.from.right & right.prox & minus.strand),]
  if(dim(from.right.facing)[1] > 0) {
    from.right.facing.dist <- from.right.facing$pos - end
    from.right.facing.dist <- rep(from.right.facing.dist, from.right.facing$n)
    right.n                <- table(from.right.facing.dist)
    from.right.dist.to.rep <- as.integer(names(right.n))
    from.right.results     <- data.frame(repeat_id, side="right", distance_to_repeat=from.right.dist.to.rep, n=as.vector(right.n), row.names=NULL)
  } else { from.right.results <- NULL }

  facing.results <- rbind(from.left.results, from.right.results)
  return(facing.results)
}

getFullRepeatsLAB <-  function(sample, repeats, repeat_id, margin) {
  start <- repeats[repeats$repeat_id == repeat_id, 4]
  end   <- repeats[repeats$repeat_id == repeat_id, 5]
  repeat_length <- end - start + 1
  
  from.left.read.in.repeat  <- (sample$pos > (start + margin)) & (sample$pos <= end) & (sample$strand == "+")
  from.right.read.in.repeat <- (sample$pos < (end - margin)) & (sample$pos >= start) & (sample$strand == "-")
  
  full_repeat <- sample[(from.left.read.in.repeat | from.right.read.in.repeat),]
  if(dim(full_repeat)[1] > 0) {
    n <- sum(full_repeat$n)
    fullrepeats <- data.frame(repeat_id, n, row.names=NULL)
  } else { fullrepeats <- NULL }
  return(fullrepeats)
}

