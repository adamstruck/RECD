##' Calculate the ratio of facing to flanking read pairs for each repeat. 
##'
##' ratio = [ facing / ( flanking + facing ) ] * 100
##' @title calculatePercentFacing
##' @param facing_file .facing file
##' @param flanking_file .flanking file
##' @param output_dir  path to output results
##' @return TSV file
##' @author Adam Struck - Intern
##' @export
calculatePercentFacing <- function(facing_file, flanking_file, output_dir) {

  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }

  pf.file.header <- c("repeat.id", "ref.repeat.length", "percent.facing", "n.total")
  fh <- openOutputFileHandle(facing_file, suffix="percent_facing", header=pf.file.header, overwrite=TRUE)
  
  facing   <- read.table(facing_file, sep="\t", header=TRUE, as.is=TRUE)
  flanking <- read.table(flanking_file, sep="\t", header=TRUE, as.is=TRUE)
  
  repeats <- unique(facing$repeat.id)
  results <- NULL
  if(length(repeats) > 0) {
    for (i.rep in 1:length(repeats)) {
      if ((i.rep %% 50 == 0) == TRUE) {
        message("Working on repeat: ", i.rep, "/", length(repeats), " ", Sys.time())
      }
      repeat.id <- repeats[i.rep]
      rep <- unlist(strsplit(as.character(repeat.id), "[:-]"))
      repeat.length <- as.numeric(rep[3]) - as.numeric(rep[2]) + 1
    
      facing.subset   <- facing[(facing$repeat.id == repeat.id),]
      flanking.subset <- flanking[(flanking$repeat.id == repeat.id),]
      
      n.facing   <- sum(facing.subset$n)
      n.flanking <- sum(flanking.subset$n)
      total <- n.facing + n.flanking
    
      if( n.facing >= 1 ) {
        percent.facing <- ((n.facing / (n.facing + n.flanking)) * 100)
        repeat.results <- data.frame(repeat.id, repeat.length, percent.facing, total)
        results <- rbind(results, repeat.results)
      }
    }
  }
  if(!is.null(results)) {
    ordered.results <- results[order(results$percent.facing, decreasing=TRUE, na.last=TRUE),]
    write.table(ordered.results, fh,
                quote=FALSE,
                sep="\t",
                col.names=FALSE,
                row.names=FALSE)
  }
}
