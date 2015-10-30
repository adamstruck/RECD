##' Makes the necessary input counts table and fdata table for running DEXSeq
##' @title makeDEXSeqTables
##' @param facing_file_list list of all .facing sample files
##' @param flanking_file_list list of all .flanking sample files
##' @param output_dir path to and name of file to output results to. 
##' @return TSV file
##' @author Adam Struck - Intern
##' @export
makeDEXSeqTables <- function(facing_file_list, flanking_file_list, output_dir="~", prefix) {

  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }
  
  facing.file.list   <- as.vector(facing_file_list)
  flanking.file.list <- as.vector(flanking_file_list)
  
  countsTable <- NULL
  fdata <- NULL
  for (i.file in 1:length(facing.file.list)) {
    message("Working on file: ", i.file , "/", length(facing.file.list), " ", Sys.time())
    facing_file <- facing.file.list[i.file]
    file_split  <- unlist(strsplit(facing_file, "/"))
    filename    <- file_split[length(file_split)]
    base_split  <- unlist(strsplit(filename, "[.]"))
    facing.basename <- base_split[1]
    facing <- read.table(facing_file, sep="\t", header=TRUE, as.is=TRUE)
    
    flanking_file <- flanking.file.list[i.file]
    file_split <- unlist(strsplit(flanking_file, "/"))
    filename   <- file_split[length(file_split)]
    base_split <- unlist(strsplit(filename, "[.]"))
    flanking.basename <- base_split[1]
    flanking <- read.table(flanking_file, sep="\t", header=TRUE, as.is=TRUE)

    if (facing.basename == flanking.basename) {
      message("Working on flanking file: ", flanking.basename, " & facing file: ", facing.basename)

      repeats <- unique(facing$repeat.id)
      
      sample.fdata  <- NULL
      sample.counts <- NULL
      for (i.rep in 1:length(repeats)) {
        if ((i.rep %% 50 == 0) == TRUE) {
          message("Working on repeat: ", i.rep, "/", length(repeats), " ", Sys.time())
        }
        repeat.id    <- repeats[i.rep]
        repeat.chr   <- unlist(strsplit(repeat.id, "[:-]"))[1]
        repeat.start <- unlist(strsplit(repeat.id, "[:-]"))[2]
        repeat.end   <- unlist(strsplit(repeat.id, "[:-]"))[3]
        
        facing.subset   <- facing[(facing$repeat.id == repeat.id), , drop=FALSE]
        flanking.subset <- flanking[(flanking$repeat.id == repeat.id), , drop=FALSE]
        
        n.facing   <- sum(facing.subset$n)
        n.flanking <- sum(flanking.subset$n)
        
        if (n.facing >= 1) {
          repeat.facing   <- sprintf("%s:facing", repeat.id)
          repeat.flanking <- sprintf("%s:flanking", repeat.id)

          facing.fdata   <- data.frame(repeat.facing, repeat.id, "facing" , repeat.chr, repeat.start, repeat.end, "*")
          names(facing.fdata) <- c("repeat.id", "geneID", "ExonID", "chr", "start", "end", "strand")
          flanking.fdata <- data.frame(repeat.flanking, repeat.id, "flanking" , repeat.chr, repeat.start, repeat.end, "*")
          names(flanking.fdata) <- c("repeat.id", "geneID", "ExonID", "chr", "start", "end", "strand")
          repeat.fdata <- rbind(facing.fdata, flanking.fdata)

          facing.counts   <- data.frame(repeat.facing, n.facing)
          names(facing.counts) <- c("repeat.id", facing.basename)
          flanking.counts <- data.frame(repeat.flanking, n.flanking)
          names(flanking.counts) <- c("repeat.id", flanking.basename)
          repeat.counts <- rbind(facing.counts, flanking.counts)
          
          if (is.null(sample.fdata) == TRUE ) {
            sample.fdata <- repeat.fdata
          } else {
            sample.fdata <- rbind(sample.fdata, repeat.fdata)
          }
          if (is.null(sample.counts) == TRUE ) {
            sample.counts <- repeat.counts
          } else {
            sample.counts <- rbind(sample.counts, repeat.counts)
          }
        }
      }
    }
    if (is.null(fdata) == TRUE ) {
      fdata <- sample.fdata
    } else {
      fdata <- merge(fdata, sample.fdata, all=TRUE)
    }
    if (is.null(countsTable) == TRUE ) {
      countsTable <- sample.counts
    } else {
      countsTable <- merge(countsTable, sample.counts, all=TRUE, by="repeat.id")
    }
  }

  rownames(countsTable) <- countsTable[,1]
  countsTable[,1] <- NULL
  countsTable[is.na(countsTable)] <- 0

  countsFile <- sprintf("%s.countsTable", prefix)
  write.table(countsTable,
              file=countsFile,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=TRUE)
  
  message("Created counts table at ", Sys.time())

  rownames(fdata) <- fdata[,1]
  fdata[,1] <- NULL
  fdataFile <- sprintf("%s.fdata", prefix)

  write.table(fdata,
              file=fdataFile,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=TRUE)

  message("Created fdata table at ", Sys.time())
}

