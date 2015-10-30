##' This package is used to download all tables from Satellog and store them in a local MySQL database
##'
##' .. content for \details{} ..
##' @title downloadSatellog
##' @param URL url for database download
##' @param DIR desired output path
##' @param DBNAME desired name of database that will be created to store the data
##' @param USR user name for database access
##' @return untarred Satellog table files and MySQL database
##' @author Adam Struck
downloadSatellog <- function(URL="http://satellog.bcgsc.ca/db_files/satellog_db_tables.tar.gz", DIR, DBNAME="satellog_db", USR="root") {
  setwd(DIR)
  file <- getBinaryURL(URL)
  out.path <- file.path(DIR, "satellog_db_tables.tar.gz")
  writeBin(file, out.path)
  db.path <- file.path(DIR, "satellog_db_tables")
  system(sprintf("mkdir %s", db.path))
  untar(out.path, compressed=T, verbose=T, tar='tar')
  db.files <- list.files(path=db.path, full.names=T)

  system(sprintf("mysqladmin -u %s create %s", USR, DBNAME))
  
  for(file in db.files) {
    message('Im working on', file)    
    system(sprintf("sed 's/type=MyISAM/engine=MyISAM/g' %s | mysql -u %s %s", file, USR, DBNAME))
  }
}

