##' Sorts read pairs into three categories: flanking, facing, and fully repetitive.
##'
##' Description of output files:
##' flanking:          read pair completely flank the repeat or are anchored outside the repeat (minimum = margin)
##' facing:            one read flanks or is anchored outside the repeat (minimum = margin) while its pair aligned inside the repeat (or didn't meet minimum anchor requirements)
##' fully_repetitive:  both reads in pair aligned inside the repeat (or didn't meet minimum anchor requirements)
##' @title mapATR [mapAlignmentsToRepeats]
##' @param con database connection
##' @param table_name name of table in database
##' @param input_file inner_mate_range.tsv file (make sure chr names are compatible with your database)
##' @param imd_dist .insert_length_dist file (i.e. the inner mate distance distribution file for the sample)
##' @param upper_limit values between 0 and 1 will be used to set prob in quantile() and then used to calculate the maximum allowed distance from the repeat. Values greater than 1 will be applied as the maximum. Read pairs whose ends are farther from the repeat than the value set by this parameter will be not be used. [default= 0.95]
##' @param margin minimum bp a read must be anchored outside the repeat [default= read_length/4]
##' @param read_length read length
##' @param min_repeat_length minimum repeat length in reference genome [default= 12] 
##' @param output_dir output directory to write results [default= "."]
##' @param append_chr appends "chr" to chromosome names of inner_mate_ranges file if TRUE [default= TRUE]
##' @param overwrite should output files be overwritten? [default= TRUE]
##' @param margin maximum bp overlap between each read and a given repeat repeat [default= read_length/2]
##' @return _flanking, _facing, and _fully_repetitive .tsv files
##' @author Adam Struck
##' @export
mapATR <- function(con, table_name, input_file, imd_dist=NULL, upper_limit=0.95, margin=(read_length/4), read_length, min_repeat_length=12, output_dir=".", append_chr=FALSE, overwrite=TRUE) {

  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }
  
  chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

  message("Reading in ", input_file, " at ", Sys.time())

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
    colnames(sample) <- c("chr", "plus", "minus", "strand", "n")

    passed <- passed + 1
    flanking.file.header   <- c("repeat.id", "inner.mate.distance", "n")
    facing.file.header     <- c("repeat.id", "side", "distance", "n")
    fullrepeat.file.header <- c("repeat.id", "n")
    if(passed == 1) {
      overwrite.val <- overwrite
      fh.flanking   <- openOutputFileHandle(input_file, suffix="flanking", header=flanking.file.header, overwrite=overwrite.val)
      fh.facing     <- openOutputFileHandle(input_file, suffix="facing", header=facing.file.header, overwrite=overwrite.val)
      fh.fullrepeat <- openOutputFileHandle(input_file, suffix="fullrepeat", header=fullrepeat.file.header, overwrite=overwrite.val)
    } else {
      overwrite.val <- FALSE
      fh.flanking   <- openOutputFileHandle(input_file, suffix="flanking", header=flanking.file.header, overwrite=overwrite.val)
      fh.facing     <- openOutputFileHandle(input_file, suffix="facing", header=facing.file.header, overwrite=overwrite.val)
      fh.fullrepeat <- openOutputFileHandle(input_file, suffix="fullrepeat", header=fullrepeat.file.header, overwrite=overwrite.val)
    }

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
    
    sample.pairs.overlap    <- sample[sample$minus < sample$plus,]
    sample.pairs.overlap.gr <- GRanges(seqnames=Rle(sample.pairs.overlap$chr), ranges=IRanges(sample.pairs.overlap$minus, sample.pairs.overlap$plus), mcols=sample.pairs.overlap$n)

    sample    <- sample[sample$minus >= sample$plus,]
    sample.gr <- GRanges(seqnames=Rle(sample$chr), ranges=IRanges(sample$plus, sample$minus), mcols=sample$n)

    repeats      <- dbGetQuery(con, sprintf("select * from %s where chr='%s'", table_name, chr))
    repeats      <- repeats[(repeats$end-repeats$start) > min_repeat_length,]
    repeats.gr   <- GRanges(seqnames=Rle(repeats$chr), ranges=IRanges(repeats$start, repeats$end))

    message("Read in repeats on ", chr, " from table ", table_name, " at ", Sys.time())
    message(paste(capture.output(gc()), collapse="\n"))

    overlaps     <- BiocGenerics::as.data.frame(findOverlaps(sample.gr, repeats.gr, type="any", maxgap=0L, minoverlap=1L, ignore.strand=TRUE))
    uniq.repeats <- unique(overlaps[,2])
    if(length(uniq.repeats) > 0) {
      message("Found read pairs that overlap repeats on ", chr, " ", Sys.time())
      for(i.rep in 1:length(uniq.repeats))  {
        rep <- uniq.repeats[i.rep]
        if ((i.rep %% 50 == 0) == TRUE) {
          message("Working on repeat: ", i.rep , "/", length(uniq.repeats), " ", Sys.time())
        }
        overlap.subset <- overlaps[overlaps$subjectHits == rep,]
        sample.subset  <- BiocGenerics::as.data.frame(sample.gr)[overlap.subset$queryHits,]
        repeat.subset  <- BiocGenerics::as.data.frame(repeats.gr)[rep,]
        repeat_id      <- paste(repeat.subset$seqname, paste(repeat.subset$start, repeat.subset$end, sep="-"), sep=":")
      
        flanking.results <- getFlanking(sample.subset, repeat.subset, repeat_id, margin, upper.limit)
        if(class(flanking.results) == "data.frame") {
          writeResults(flanking.results, fh.flanking)
        }
      
        facing.results <- getFacing(sample.subset, repeat.subset, repeat_id, margin, upper.limit)
        if(class(facing.results) == "data.frame") {
          writeResults(facing.results, fh.facing)
        }
      
        fullrepeat.results <- getFullRepeats(sample.subset, repeat.subset, repeat_id, margin)
        if(class(fullrepeat.results) == "data.frame") {
          writeResults(fullrepeat.results, fh.fullrepeat)
        }
      }
    }
    ### For the cases where minus < plus
    message("Working on overlaping reads... ", Sys.time())
    overlaps.po  <- BiocGenerics::as.data.frame(findOverlaps(sample.pairs.overlap.gr, repeats.gr, type="any", maxgap=0L, minoverlap=1L, ignore.strand=TRUE))    
    uniq.repeats <- unique(overlaps.po[,2])
    if(length(uniq.repeats) > 0) {
      for(i.rep in 1:length(uniq.repeats))  {
        if ((i.rep %% 50 == 0) == TRUE) {
          message("Working on repeat: ", i.rep , "/", length(uniq.repeats), " ", Sys.time())
        }
        rep <- uniq.repeats[i.rep]
    
        overlap.subset <- overlaps.po[overlaps.po$subjectHits == rep,]
        sample.subset  <- BiocGenerics::as.data.frame(sample.pairs.overlap.gr)[overlap.subset$queryHits,]
        names(sample.subset) <- c("seqnames", "end", "start", "width", "strand", "mcols")
        sample.subset$width <- -sample.subset$width
        repeat.subset  <- BiocGenerics::as.data.frame(repeats.gr)[rep,]
        repeat_id      <- paste(repeat.subset$seqname, paste(repeat.subset$start, repeat.subset$end, sep="-"), sep=":")
      
        flanking.results <- getFlanking(sample.subset, repeat.subset, repeat_id, margin, upper.limit)
        if(class(flanking.results) == "data.frame") {
          writeResults(flanking.results, fh.flanking)
        }
        
        facing.results <- getFacing(sample.subset, repeat.subset, repeat_id, margin, upper.limit)
        if(class(facing.results) == "data.frame") {
          writeResults(facing.results, fh.facing)
        }
        
        fullrepeat.results <- getFullRepeats(sample.subset, repeat.subset, repeat_id, margin)
        if(class(fullrepeat.results) == "data.frame") {
          writeResults(fullrepeat.results, fh.fullrepeat)
        }
      }
    }
    close(fh.flanking)
    close(fh.facing)
    close(fh.fullrepeat)
  }
}

getFlanking <- function(sample.subset, repeat.subset, repeat_id, margin, upper.limit) {
  plus.read.anchored  <- sample.subset$start <= (repeat.subset$start + margin)
  minus.read.anchored <- sample.subset$end >= (repeat.subset$end - margin)

  plus.prox  <- ((repeat.subset$start + margin) - sample.subset$start) <= upper.limit
  minus.prox <- (sample.subset$end - (repeat.subset$end - margin)) <= upper.limit

  plus.conditions  <- plus.read.anchored & plus.prox
  minus.conditions <- minus.read.anchored & minus.prox
  
  flanking <- sample.subset[(plus.conditions & minus.conditions),]
  if(dim(flanking)[1] > 0) {
    imds                <- rep(flanking$width, flanking$mcols)
    n                   <- table(imds)
    inner.mate.distance <- as.integer(names(n))
    flanking.results <- data.frame(repeat_id, inner.mate.distance, n=as.vector(n), row.names=NULL)
  } else { flanking.results <- NULL }
  return(flanking.results)
}

getFacing <- function(sample.subset, repeat.subset, repeat_id, margin, upper.limit) {
  plus.read.in.repeat <- sample.subset$start > (repeat.subset$start + margin)
  minus.read.outside  <- sample.subset$end >= (repeat.subset$end - margin)
  minus.prox          <- (sample.subset$end - (repeat.subset$end - margin)) <= upper.limit
  facing.from.right   <- plus.read.in.repeat & minus.read.outside & minus.prox
  
  plus.read.outside    <- sample.subset$start <= (repeat.subset$start + margin)
  minus.read.in.repeat <- sample.subset$end < (repeat.subset$end - margin)
  plus.prox            <- ((repeat.subset$start + margin) - sample.subset$start) <= upper.limit
  facing.from.left     <- plus.read.outside & minus.read.in.repeat & plus.prox

  facing <- sample.subset[(facing.from.left | facing.from.right),]
  if(dim(facing)[1] > 0) {
    from.left.facing         <- facing[facing$end < (repeat.subset$end - margin),]
    if(dim(from.left.facing)[1] > 0) {
      from.left.facing.dist  <- from.left.facing$start - repeat.subset$start 
      from.left.facing.dist  <- rep(from.left.facing.dist, from.left.facing$mcols) 
      left.n                 <- table(from.left.facing.dist)
      from.left.dist.to.rep  <- as.integer(names(left.n))
      from.left.results      <- data.frame(repeat_id, side="left", distance_to_repeat=from.left.dist.to.rep, n=as.vector(left.n), row.names=NULL)
    } else { from.left.results <- NULL }

    from.right.facing        <- facing[facing$start > (repeat.subset$start + margin),]
    if(dim(from.right.facing)[1] > 0) {
      from.right.facing.dist <- from.right.facing$end - repeat.subset$end
      from.right.facing.dist <- rep(from.right.facing.dist, from.right.facing$mcols)
      right.n                <- table(from.right.facing.dist)
      from.right.dist.to.rep <- as.integer(names(right.n))
      from.right.results     <- data.frame(repeat_id, side="right", distance_to_repeat=from.right.dist.to.rep, n=as.vector(right.n), row.names=NULL)
    } else { from.right.results <- NULL}
  
    facing.results <- rbind(from.left.results, from.right.results)
  } else { facing.results <- NULL }
  return(facing.results)
}

getFullRepeats <-  function(sample.subset, repeat.subset, repeat_id, margin) {
  plus.read.in.repeat  <- sample.subset$start > (repeat.subset$start + margin)
  minus.read.in.repeat <- sample.subset$end < (repeat.subset$end - margin)
  
  full_repeat <- sample.subset[(plus.read.in.repeat & minus.read.in.repeat),]
  if(dim(full_repeat)[1] > 0) {
    n <- sum(sample.subset$mcols)
    fullrepeats <- data.frame(repeat_id, n, row.names=NULL)
  } else { fullrepeats <- NULL }
  return(fullrepeats)
}


