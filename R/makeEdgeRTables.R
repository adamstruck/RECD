##' Creates the necessary counts and design tables for running edgeR.
##' @title makeEdgeRTables
##' @param facing_file_list list of all .facing sample files
##' @param flanking_file_list list of all .flanking sample files
##' @param pheno_table ExpressionPlot type pheno table 
##' @param min_total_counts miniumum total number of read pairs associated with a given repeat across all samples
##' @param output_dir path to and name of file to output results to. 
##' @return TSV file
##' @author Adam Struck - Intern
##' @export
makeEdgeRTables <- function(facing_file_list, flanking_file_list, pheno_table, min_total_counts=50, output_dir="~", prefix) {

  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }
  
  facing.file.list   <- as.vector(facing_file_list)
  flanking.file.list <- as.vector(flanking_file_list)
  
  countsTable <- NULL
  
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
      message("All files were paired properly")

      repeats <- unique(facing$repeat.id)
      
      sample.results <- NULL
      for (i.rep in 1:length(repeats)) {
        if ((i.rep %% 50 == 0) == TRUE) {
          message("Working on repeat: ", i.rep, "/", length(repeats), " ", Sys.time())
        }
        repeat.id    <- repeats[i.rep]
        repeat.start <- unlist(strsplit(repeat.id, "[:-]"))[2]
        repeat.end   <- unlist(strsplit(repeat.id, "[:-]"))[3]
        
        facing.subset   <- facing[(facing$repeat.id == repeat.id), , drop=FALSE]
        flanking.subset <- flanking[(flanking$repeat.id == repeat.id), , drop=FALSE]
        
        n.facing   <- sum(facing.subset$n)
        n.flanking <- sum(flanking.subset$n)

        if (n.facing >= 1) {
          repeat.results <- data.frame(repeat.id, n.facing, n.flanking)
          facing.name    <- sprintf("%s.facing", facing.basename)
          flanking.name  <-  sprintf("%s.flanking", flanking.basename)
          names(repeat.results) <- c("repeat.id", facing.name, flanking.name)

          if (is.null(sample.results) == TRUE ) {
            sample.results <- repeat.results
          } else {
            sample.results <- rbind(sample.results, repeat.results)
          }
        }
      }
    }
    if (is.null(countsTable) == TRUE ) {
      countsTable <- sample.results
    } else {
      countsTable <- merge(countsTable, sample.results, all=TRUE, by="repeat.id")
    }
  }

  rownames(countsTable) <- countsTable[,1]
  countsTable[,1] <- NULL
  countsTable <- countsTable[which(rowSums(countsTable, na.rm=TRUE) >= min_total_counts),]
  countsTable[is.na(countsTable)] <- 0

  countsFile <- sprintf("%s.countsTable", prefix)
  write.table(countsTable,
              file=countsFile,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=TRUE)
  
  message("Created counts table at ", Sys.time())
  
  pdata <- read.table(pheno_table, header=TRUE, sep="\t", as.is=TRUE)
  
  sample.names <- names(countsTable)
  tmp <- data.frame(do.call('rbind', strsplit(sample.names, "[.]")))
  names(tmp) <- c("SAMID", "Read.Class")
  
  designTable <- merge(tmp, pdata, all=TRUE, by="SAMID")
  rownames(designTable) <- sample.names

  designFile <- sprintf("%s.designTable", prefix)
  write.table(designTable,
              file=designFile,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=TRUE)
  
  message("Created design table at ", Sys.time())

}

