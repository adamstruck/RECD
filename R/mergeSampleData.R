##' mergeFiles take in a list of files as its input and merges them into one file.
##'
##' Description of output file format:
##' rows = repeat.id
##' columns = attribute (such as mean.irl) from each input file; column names are derived from the input filename (all text preceding the first period in the filename). For CGP data, this should be the SAM identification. 
##' @title mergeSampleData
##' @param file_list character vector of sample files to merge
##' @param attribute feature to merge from files
##' @param output_file path to and name of desired output file
##' @return Summary of implied repeat lengths for all samples
##' @author Adam Struck - Intern
##' @export
mergeSampleData <- function(input_files, attribute, output_file) {
  output.table <- NULL
  files <- as.vector(input_files)
  
  for(i.file in 1:length(files))  {
      file <- files[i.file]
      message("Working on file: ", i.file , "/", length(files), " ", Sys.time())
      file_split <- unlist(strsplit(file, "/"))
      filename <- file_split[length(file_split)]
      base_split <- unlist(strsplit(filename, "[.]"))
      basename <- base_split[1]

      data <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
      data <- data.frame(data[["repeat.id"]], data[[attribute]])
      names(data) <- c("repeat.id", sprintf("%s", basename))
      
      if(is.null(output.table) == TRUE) {
        output.table <- data
        message("Finished processing ", file, " ", Sys.time())
      } else {
        output.table <- merge(output.table, data, all=TRUE, by="repeat.id")
        writeMergedResults(output.table, output_file)
        message("Finished processing ", file, " ", Sys.time())
      }
    }
}

writeMergedResults <- function(results.table, output_file) {
  write.table(results.table,
              col.names = TRUE,
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              file=output_file)
}
