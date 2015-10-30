##' Try to detect repeats that have gotten longer than the typical inner mate distance
##'
##' Uses the DEXSeq package to try and detect differential read class frequencies (facing versus flanking) between tumor and normal samples. 
##' @title DEXSeqAnalysis
##' @param countsTable 
##' @param fdata 
##' @param pheno_table 
##' @param condition condition to base the differential analysis on.
##' @param grouping list specifying a specific variable name and level to subset the data by. Only samples that have this trait (level) will be examined. [ Example: grouping = list(name="Tissue.Group.Oncology", level="Colon") ]
##' @param output_file path to and name of file to output results to. 
##' @return TSV file
##' @author Adam Struck - Intern
##' @export
DEXSeqAnalysis <- function(counts_table, fdata_table, pheno_table, condition, grouping, output_file) {

  countsTable <- read.table(counts_table, header=TRUE, sep="\t")
  fdata <- read.table(fdata_table, header=TRUE, sep="\t")

  pdata <- read.table(pheno_table, header=TRUE, row.names=1, sep="\t")
  pdata[, names(pdata)] <- lapply(pdata[, names(pdata)], factor)
  names(pdata)[names(pdata) == condition] <- "test_condition"

  message("pdata: ", paste(dim(pdata), collapse = " "))
  message("countsTable: ", paste(dim(countsTable), collapse = " "))
  
  if(is.list(grouping) == TRUE) {
    varname <- grouping$name
    if ( !varname %in% names(pdata) ) {
      stop(varname, " not in names(pdata) = [", paste(collapse = " ", names(pdata)), "]")
    }
    if ( "level" %in% names(grouping) ) {
      level <- grouping$level
      if ( !level %in% levels(pdata[[varname]]) ) {
        stop(level, sprintf(" not in levels(pdata$%s)", varname), " = [", paste(collapse = " ", levels(pdata[[level]])), "]")
      }
      limit <- pdata[[varname]] == level
      pdata <- pdata[limit, , drop = FALSE]
    }
  }
  
  countsTable <- countsTable[,rownames(pdata)]
  countsTable <- countsTable[which(rowSums(countsTable) > 0),]

  fdata <- fdata[rownames(countsTable),]
  
  message("pdata: ", paste(dim(pdata), collpase = " "))
  message("countsTable: ", paste(dim(countsTable), collpase = " "))
  message("Loaded input files...")
  
  ecs <- newExonCountSet(countData = countsTable,
                         design = pdata,
                         geneIDs = fdata$geneID,
                         exonIDs = fdata$ExonID)

  message("Created new exon count set...")
  
  ecs <- estimateSizeFactors(ecs)

  dispersionformula <- count ~ sample + exon * test_condition
  ecs <- estimateDispersions(ecs, formula=dispersionformula, minCount=10)
  ecs <- fitDispersionFunction(ecs)

  message("Estimated dispersions...")
  
  formula0 <- count ~ sample + exon + test_condition
  formula1 <- count ~ sample + exon + test_condition * I(exon == exonID)
  
  ecs <- testForDEU(ecs, formula0 = formula0, formula1 = formula1)
  ecs <- estimatelog2FoldChanges(ecs, fitExpToVar = "test_condition")
  res <- DEUresultTable(ecs)
  res <- res[order(res$pvalue), ]

  numVars <- sapply(res, is.numeric)
  res[numVars] <- lapply(res[numVars], round, digits=4)
  
  message("Writing results to: ", output_file)
  
  write.table(res,
              file=output_file,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=FALSE)
}

