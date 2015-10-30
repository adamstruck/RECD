##' Connect to a MySQL or SQLite database containing SSRs and selelct entries within a regionon a chromosome. 
##'
##' .. content for \details{} ..
##' @title RepeatsInRegion
##' @param con database connection
##' @param chr chromosome
##' @param start starting genomic coordinates
##' @param end  ending genomic coordinates
##' @return data.frame
##' @author Adam Struck
##' @export
RepeatsInRegion <- function(con, chr, start, end) {
  results <- dbGetQuery(con, sprintf("SELECT * from repeats WHERE chr = '%s' AND start >= %d AND end <= %d", CHR, START, END))
  return(results)
  dbDisconnect(con)
}



