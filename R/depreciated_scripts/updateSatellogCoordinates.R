##' Updates all of the genomic coordinates in the repeats table of your Satellog database to the current hg19 build.
##'
##' Since Satellog was developed using hg16 v34 of the human genome, this package updates the coordinates to the most recent hg19 build using the liftOver tool from UCSC. Updated coordinates are inserted into SQLite database
##' @title updateSatellogCoordinates
##' @param con MySQL database connection
##' @param con2 SQLite database connection
##' @param exe.path location of the liftOver executable and chain files
##' @return SQLite database
##' @author Adam Struck
updateSatellogCoordinates <- function(con, con2, exe.path="./") { 
  chrs <- dbGetQuery(con, "SELECT DISTINCT chr FROM repeats")$chr
  for(chr in chrs)  {
    message("I'm working on chr", chr)
    chr.data <- dbGetQuery(con, sprintf("SELECT * FROM repeats WHERE chr = '%s'", chr))
    chr.bed <- data.frame(chr.data[c("chr","start","end")], strand=".", chr.data["rep_id"])
    chr.bed$chr <- sprintf("chr%s", chr)
    write.table(chr.bed, sprintf("chr%s_hg16.bed", chr), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    system(sprintf("%sliftOver chr%s_hg16.bed %shg16ToHg19.over.chain chr%s_hg19.bed chr%s_hg19_unmapped.txt", exe.path, exe.path, chr, chr, chr))

### Rejoin updated coordinates with main dataframe
    hg19_chr.bed <- read.table(sprintf("chr%s_hg19.bed", chr), sep="\t")
    header <- colnames(chr.data)
    hg19_chr.data <- merge(hg19_chr.bed, chr.data, by.x = "V5", by.y = "rep_id")
    hg19_chr.data <- data.frame(hg19_chr.data[c("V5", "V1", "V2", "V3", "period", "unit", "class", "seq", "length", "pvalue")])
    colnames(hg19_chr.data) <- header
    if(dbExistsTable(con2, "repeats")){
      dbWriteTable(con2, "repeats", hg19_chr.data, append=TRUE)
    } else {
      dbWriteTable(con2, "repeats", hg19_chr.data)
    }
    system(sprintf("rm chr%s*", chr))
  }  
  dbSendQuery(con2, "CREATE INDEX repeats_rep_id ON repeats(rep_id)")
  dbDisconnect(con2)
}

### Update MySQL database
##     primary <-  'rep_id'
##     vars <- header[2:10]
##     pastedvars <- paste("'", apply(chr.data.hg19[, vars], 1, paste, collapse="', '"), "'", sep="")
##     varlist <- paste("repeats", "(", paste(c(primary, vars), collapse=", "), ")", sep="") 
##     datastring <- paste("(", paste(paste(chr.data.hg19[, primary], pastedvars, sep=", "), collapse="), ("), ")", sep="")
##     toupdate <- paste(paste(vars, "=VALUES(", vars, ")", sep=""), collapse=", ")
##     sprintf("INSERT INTO %s VALUES %s ON DUPLICATE KEY UPDATE %s",varlist, datastring, toupdate)
##     dbSendQuery(con, sprintf("INSERT INTO %s VALUES %s ON DUPLICATE KEY UPDATE %s",varlist, datastring, toupdate))
##     dbDisconnect(con)

# con <- dbConnect(MySQL(), user='root', dbname='satellog_db')

