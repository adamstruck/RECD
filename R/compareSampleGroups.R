##' Perform a Student's t-test or a Spearmans's correlation test on different groups of your data.
##'
##' See below:
##' @title compareSampleGroups
##' @param merged_sample_data A data frame of numerical data that you wish to group and compare using a student's t-test. 
##' @param pheno_table A data frame giving the group labels (phenotypic data) for each point in ‘merged_sample_data’. It should have the same number of rows as ‘merged_sample_data’ is long. If ‘merged_sample_data’ is named and ‘pdata’ has rownames then they are checked for equality.
##' @param grouping A list specifying the levels of grouping your data based on the ‘pdata’. Each entry of the list has several componenets:
##'                 $name The name of a field from ‘pdata’
##'                 $level a specific level of the specified $name to filter the data by. Only samples with this level will be considered.
##'                 ONLY ONE $NAME CAN BE SPECIFIED WITHOUT A CORRESPONDING LEVEL. 
##' @param test "ttest" or "cor"; perform either a student's t test on the supplied data or calculate the Spearmans' correlation between patient tumor/normal samples. 
##' @param matched_patient_tissues logical, if TRUE only samples with matched tumor and normal samples will be considered. 
##' @param output_file Desired name of output file
##' @return TSV file containing t-test statistics for each repeat
##' @author Adam Struck - Intern
##' @export
compareSampleGroups <- function(merged_sample_data, pheno_table, grouping=list(list(name="plot.metaclass")), test="ttest", matched_patient_tissues=FALSE, output_file) {

  samples <- read.table(merged_sample_data, header=TRUE, row.names=1, sep="\t")
  pdata <- read.table(pheno_table, header=TRUE, row.names=1, sep="\t")

  if ( test == "cor" ) {
    matched_patient_tissues = TRUE
  }
  
  grouping <- lapply(grouping, function(gp) {
    varname <- gp$name
    if ( !varname %in% names(pdata) ) {
      stop(varname, " not in names(pdata) = [", paste(collapse = " ", names(pdata)), "]") }
    if ( "level" %in% names(gp) ) {
      level <- gp$level
      if ( !level %in% levels(pdata[[varname]]) ) {
        stop(level, sprintf(" not in levels(pdata$%s)", varname), " = [", paste(collapse = " ", levels(pdata[[level]])), "]") }
      limit <- pdata[[varname]] == level
      pdata <<- pdata[limit, , drop = FALSE] }
    gp
  })

  if( matched_patient_tissues == TRUE ) {
    pdata  <- pdata[which(pdata$Patient.ID %in% names(which(table(pdata$Patient.ID) == 2))),]
    pdata  <- pdata[order(pdata$Patient.ID),]
  }
  
  subsets <- NULL
  for ( group in grouping ) {
    if ( !"level" %in% names(group) ) {
      lvls <- levels(pdata[[group$name]])
      if ( !length(lvls) == 2 ) {
        stop("group ", group$name, " does not have 2 levels, it has: ", length(lvls)) }
      else {
        for ( i.lvl in 1:length(lvls) ) {
          lvl <- lvls[i.lvl]
          subsets <- c(subsets, paste("group.",lvl,sep=""))
          assign(paste("group.",lvl,sep=""), samples[,rownames(pdata[pdata[[group$name]] == lvl,])])
        }
      }
    }
  }

  if ( !length(subsets) == 2 ) {
    stop("Missing level on at least one of your grouping variables. /n A level must be specified on all grouping variables EXCEPT the one you want to make comparisons between.")
  }

  if( matched_patient_tissues == TRUE ) {
    if( dim(get(subsets[1]))[2] != dim(get(subsets[2]))[2] ) {
      tmp1 <- get(subsets[1])[,which(pdata[names(get(subsets[1])),]$Patient.ID %in% pdata[names(get(subsets[2])),]$Patient.ID)]
      tmp2 <- get(subsets[2])[,which(pdata[names(get(subsets[2])),]$Patient.ID %in% pdata[names(get(subsets[1])),]$Patient.ID)]
      tmp <- cbind(tmp1, tmp2)
      first.subset <- 1:ncol(tmp1)
      second.subset <- (1+ncol(tmp1)):(ncol(tmp1)+ncol(tmp2))
    } else {
      tmp <- cbind(get(subsets[1]), get(subsets[2]))
      first.subset <- 1:ncol(get(subsets[1]))
      second.subset <- (1+ncol(get(subsets[1]))):(ncol(get(subsets[1]))+ncol(get(subsets[2])))
    }
  } else {
    tmp <- cbind(get(subsets[1]), get(subsets[2]))
    first.subset <- 1:ncol(get(subsets[1]))
    second.subset <- (1+ncol(get(subsets[1]))):(ncol(get(subsets[1]))+ncol(get(subsets[2])))
  }

  ttest <- function(x) {
    out <- try( t.test(x[first.subset], x[second.subset]) )
    if( is( out, "try-error" )) {
      data.frame(statistic = NA,
                 df   = NA,
                 pval = NA)
    } else {
      data.frame(statistic = round(out$statistic, 2),
                 df   = round(out$parameter, 2),
                 pval = out$p.value)
    }
  }

  spear.cor <- function(x) {
    cor(x[first.subset], x[second.subset],
        method="spearman",
        use = "na.or.complete")
  }

  if( test == "ttest" ) {
    results <- apply(tmp, 1, ttest)
  
    comps <- do.call(rbind.data.frame, results)
    comps$adj.pval <- p.adjust(comps$pval, method="BH")
    ranked.comps <- comps[order(comps$adj.pval), , drop=FALSE]
  
    ## middle <- function(x) round(median(x, na.rm=TRUE), 2)
    ## assign(paste(unlist(strsplit(subsets[1], "[.]"))[2], ".median", sep=""), apply(get(subsets[1])[rownames(ranked.comps),], 1, middle))
    ## assign(paste(unlist(strsplit(subsets[2], "[.]"))[2], ".median", sep=""), apply(get(subsets[2])[rownames(ranked.comps),], 1, middle))

    avg <- function(x) round(mean(x, na.rm=TRUE), 2)
    assign(paste(unlist(strsplit(subsets[1], "[.]"))[2], ".mean", sep=""), apply(get(subsets[1])[rownames(ranked.comps),], 1, avg))
    assign(paste(unlist(strsplit(subsets[2], "[.]"))[2], ".mean", sep=""), apply(get(subsets[2])[rownames(ranked.comps),], 1, avg))

    ranked.comps$pval <- formatC(ranked.comps$pval, format="e", digits=3)
    ranked.comps$adj.pval <- formatC(ranked.comps$adj.pval, format="e", digits=3)

    repeat.lengths <- data.frame(do.call('rbind', strsplit(as.character(rownames(ranked.comps)), "[:-]")))
    names(repeat.lengths) <- c("chr", "start", "end")
    repeat.lengths$start  <-  as.numeric(as.character(repeat.lengths$start))
    repeat.lengths$end    <-  as.numeric(as.character(repeat.lengths$end))
    lengths <- repeat.lengths$end - repeat.lengths$start + 1
  
    name.subset1 = paste(unlist(strsplit(subsets[1], "[.]"))[2], ".mean", sep="")
    name.subset2 = paste(unlist(strsplit(subsets[2], "[.]"))[2], ".mean", sep="")
    
    ranked.comparisons <- data.frame(ranked.comps, get(name.subset1), get(name.subset2), lengths)
    names(ranked.comparisons) <- c("t.statistic", "df", "pval", "adj.pval", name.subset1, name.subset2, "repeat.length")

  } else {
    if( test == "cor" ) {
      results <- apply(tmp, 1, spear.cor)

      countPairs <- function(x) {
        x <- data.frame(x[first.subset], x[second.subset])
        x <- x[which(rowSums(is.na(x)) == 0),]
        n <- dim(x)[1]
      }

      n <- apply(tmp, 1, countPairs)
      
      ranked.comparisons <- data.frame(results, n)
      names(ranked.comparisons) <- c("correlation", "n")
      ranked.comparisons$correlation <- round(ranked.comparisons$correlation, 3)
      ranked.comparisons <- ranked.comparisons[ranked.comparisons$n >= 2, ]
      ranked.comparisons <- ranked.comparisons[order(ranked.comparisons$n, decreasing=TRUE, na.last=TRUE, ranked.comparisons$correlation),]

    }
  }
  
  write.table(ranked.comparisons,
              file=output_file,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=TRUE)

}
