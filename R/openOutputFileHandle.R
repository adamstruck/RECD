##' Open file handle for writing output. 
##'
##' @title openOutputFileHandle
##' @param input_file 
##' @param suffix 
##' @param header 
##' @param overwrite 
##' @return open file handle
##' @author Adam Struck - Intern
openOutputFileHandle <- function(input_file, suffix, header, overwrite) {
  file_split <- unlist(strsplit(input_file, "/"))
  filename <- file_split[length(file_split)]
  base_split <- unlist(strsplit(filename, "[.]"))
  base_name <- base_split[1]

  output_file <- sprintf("%s.%s", base_name, suffix)

  if(overwrite == TRUE) {
    if(file.exists(output_file) == TRUE) {
      message("MESSAGE: ", output_file, " was overwritten.")
    } else {
      message("MESSAGE: ", output_file, " does not exist. File will be created.")
    }
    fh <- file(output_file, open="w")
    header <- as.character(header)
    df <- data.frame(matrix(matrix(rep(1,length(header)),1),1)) 
    colnames(df) <- header
    df <- df[-1,]
    write.table(df, fh,
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t",
                quote = FALSE)
  } else {
    if(overwrite == FALSE) {
      if(file.exists(output_file) == TRUE) {
        message("MESSAGE: Found ", output_file, ", results will be appended to this file.")
        fh <- file(output_file, open="a+")
      } else {
        message("MESSAGE: ", output_file, " does not exist. File will be created.")
        fh <- file(output_file, open="w")
        header <- as.character(header)
        df <- data.frame(matrix(matrix(rep(1,length(header)),1),1)) 
        colnames(df) <- header
        df <- df[-1,]
        write.table(df, fh,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE)
      }
    }
  }
  return(fh)
}

##' Write a data frame to a file handle
##'
##' .. content for \details{} ..
##' @title writeResults
##' @param results.table 
##' @param output.fh 
##' @return tsv file
##' @author Adam Struck - Intern
writeResults <- function(results.table, output.fh) {
  write.table(results.table, output.fh,
              col.names = FALSE,
              row.names = FALSE,
              sep = "\t",
              quote = FALSE)
}
