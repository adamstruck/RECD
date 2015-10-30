##' Try to detect repeats that have gotten longer than the typical inner mate distance
##'
##' Uses the edgeR package to try and detect differential read class frequencies (facing versus flanking) between tumor and normal samples. 
##' @title edgeRAnalysis
##' @param countsTable 
##' @param designTable 
##' @param output_file path to and name of file to output results to. 
##' @return TSV file
##' @author Adam Struck - Intern
##' @export
edgeRAnalysis <- function(countsTable, designTable, output_file) {

  countsTable <- read.table(countsTable, header=TRUE, sep="\t")
  designTable <- read.table(designTable, header=TRUE, sep="\t")
  
  designTable$Tissue.Group.Oncology <- as.character(designTable$Tissue.Group.Oncology)
  designTable$Tissue.Group.Oncology[designTable$Tissue.Group.Oncology == "Urinary Bladder"] <- "Urinary_Bladder"
  designTable$Tissue.Group.Oncology <- as.factor(designTable$Tissue.Group.Oncology)
  
  ## edgeR
  cds  <- DGEList(counts=countsTable, genes=rownames(countsTable))
  keep <- rowSums(cpm(cds) > 1) >= 66
  cds  <- cds[keep,]
  cds  <- calcNormFactors(cds)

  Group <- factor( paste( designTable$Read.Class, designTable$Diagnosis.Group.Oncology, sep="."))
  #Group <- factor( paste( designTable$Read.Class, designTable$Tissue.Group.Oncology, designTable$Diagnosis.Group.Oncology, sep="."))
  designTable <- cbind(designTable, Group)
  design <- model.matrix(~0+Group, data=designTable)
  colnames(design) <- levels(Group)
  rownames(design) <- colnames(cds)

  cds <- estimateGLMCommonDisp(cds, design, verbose=TRUE)
  cds <- estimateGLMTrendedDisp(cds, design)
  cds <- estimateGLMTagwiseDisp(cds, design)
  fit <- glmFit(cds, design)

  CGP.contrast <- makeContrasts((facing.Cancer-flanking.Cancer)-(facing.Normal-flanking.Normal), levels=design)  
  lrt <- glmLRT(fit, contrast=CGP.contrast)

  Colon.contrast   <- makeContrasts((facing.Colon.Cancer-flanking.Colon.Cancer)-(facing.Colon.Normal-flanking.Colon.Normal), levels=design)
  lrt <- glmLRT(fit, contrast=Colon.contrast)

  Stomach.contrast <- makeContrasts((facing.Stomach.Cancer-flanking.Stomach.Cancer)-(facing.Stomach.Normal-flanking.Stomach.Normal), levels=design)
  lrt <- glmLRT(fit, contrast=Stomach.contrast)

  Liver.contrast   <- makeContrasts((facing.Liver.Cancer-flanking.Liver.Cancer)-(facing.Liver.Normal-flanking.Liver.Normal), levels=design)
  lrt <- glmLRT(fit, contrast=Liver.contrast)

  Lung.contrast    <- makeContrasts((facing.Lung.Cancer-flanking.Lung.Cancer)-(facing.Lung.Normal-flanking.Lung.Normal), levels=design)  
  lrt <- glmLRT(fit, contrast=Lung.contrast)

  ## PD.contrast <- makeContrasts((facing.Parkinsons_Disease-flanking.Parkinsons_Disease)-(facing.Normal-flanking.Normal), levels=design)  
  ## lrt <- glmLRT(fit, contrast=PD.contrast)
  
  edgeR.results <- BiocGenerics::as.data.frame(topTags(lrt, n=dim(lrt$table)[1]))
  names(edgeR.results) <- c("repeat.id", "logFC", "logCPM", "LR", "PValue", "FDR")
  
  write.table(edgeR.results,
              file=output_file,
              quote=FALSE,
              sep="\t",
              col.names=TRUE,
              row.names=FALSE)
  
}

## DESeq
##   cdsFull <- newCountDataSet(countsTable, designTable)
##   cdsFull <- estimateSizeFactors(cdsFull)
##   cdsFull <- estimateDispersions(cdsFull, sharingMode="gene-est-only", method="pooled-CR", modelFormula = count ~ Read.Class + Diagnosis.Group.Oncology) 
##   fit1 = fitNbinomGLMs(cdsFull, count ~ Read.Class + Diagnosis.Group.Oncology )
##   fit0 = fitNbinomGLMs(cdsFull, count ~ Read.Class )

##   pvalsGLM = nbinomGLMTest( fit1, fit0 )
##   padjGLM = p.adjust( pvalsGLM, method="BH" )

##   DEseq_results <- data.frame(fit1, pval=pvalsGLM, padj=padjGLM)
##   DEseq_results <- DEseq_results[with(DEseq_results, order(DEseq_results$padj)),]
 
