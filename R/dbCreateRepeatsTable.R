##' Creates a table in your SQLite database.
##'
##' .. content for \details{} ..
##' @title dbCreateRepeatsTable
##' @param con database connection
##' @param file repeats data to be stored (e.g. output from parseTRFoutput.py)
##' @param table_name name of table to create
##' @return SQL/SQLite database table
##' @author Adam Struck
##' @export
dbCreateRepeatsTable  <- function(con, file, table_name) {
  data <- read.table(file, header = TRUE, sep = "\t")
  dbWriteTable(con, sprintf("%s", table_name), data)
  dbDisconnect(con)
}
