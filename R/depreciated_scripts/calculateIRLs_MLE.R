##' For repeats that have read pairs flanking them, the implied expansion and length of that repeat is determined. 
##'
##' Summary stats for the implied repeat exansions/lengths are also reported. 
##' @title calculateIRLs_MLE
##' @param flanking_file _flanking read pair file
##' @param imd_distribution .insert_len_dist file (i.e. the inner mate distance distribution for the sample)
##' @param upper_limit the value used in mapATR % processLastAlignedBase [default= 0.95]
##' @param margin the value used in mapATR % processLastAlignedBase [default= 37.5]
##' @param min_coverage minimum number of flanking read pairs that map to a given repeat required for IRL calculation [default= 10]
##' @param output_dir path to output directory
##' @return _implied_repeat_lengths file
##' @author Adam Struck - Intern
calculateIRLs_MLE <- function(flanking_file, imd_distribution, upper_limit=0.95, margin=37.5, min_coverage=10, output_dir) {
  if( file.exists(output_dir) ) {
    setwd(output_dir)
  } else {
    dir.create(output_dir)
    setwd(output_dir)
  }
 
  irl.file.header <- c("repeat.id", "ref.repeat.length", "implied.repeat.length", "n")
  fh <- openOutputFileHandle(flanking_file, suffix="implied_repeat_lengths_MLE", header=irl.file.header, overwrite=TRUE)
    
  imd.distribution <- read.table(imd_distribution, header=TRUE, sep="\t")
  distribution     <- rep(imd.distribution$inner.mate.distance, imd.distribution$n)

  if( upper_limit <= 1 ) {
    d95 <- quantile(distribution, prob=upper_limit)
  } else {
    d95 <- upper_limit
  }

  imd.distribution          <- imd.distribution[imd.distribution$inner.mate.distance <= 10000,]
  imd.prob.distribution     <- calculate.IMD.prob.distribution(imd.distribution)
  log.imd.prob.distribution <- log(imd.prob.distribution)
  
  possible.repeat.lengths <- 0 : 1000
  
  mass.sums     <- IMD.mass.sums(possible.repeat.lengths, imd.prob.distribution, margin, d95)
  log.mass.sums <- log(mass.sums)

  flanking <- read.table(flanking_file, header=TRUE, sep="\t")  
  repeats  <- unique(flanking$repeat.id)

  if(length(repeats) > 0) {
    for(i.rep in 1:length(repeats))  {
      if ((i.rep %% 50 == 0) == TRUE) { message("Working on repeat: ", i.rep , "/", length(repeats), " ", Sys.time()) }

      repeat.id         <- repeats[i.rep]
      rep.info          <- unlist(strsplit(as.character(repeat.id), "[:-]"))
      ref.repeat.length <- as.numeric(rep.info[3]) - as.numeric(rep.info[2]) + 1

      flanking.subset     <- flanking[flanking$repeat.id == repeat.id,]
      flanking.imd.counts <- setNames(flanking.subset$n, flanking.subset$inner.mate.distance)
      total.n.flanking    <- sum(flanking.imd.counts)

      if(total.n.flanking >= min_coverage) {
        
        likelihoods          <- IMD.logLikelihoodOffset(flanking.imd.counts, log.imd.prob.distribution, log.mass.sums, ref.repeat.length)
        ## implied.repeat.length <- as.numeric(names(which(likelihoods == max(likelihoods, na.rm=TRUE))))
        
        likelihoods.loess    <- loess(likelihoods ~ possible.repeat.lengths, span=0.20)
        obs.repeat.lengths   <- likelihoods.loess$x
        smoothed.likelihoods <- predict(likelihoods.loess, obs.repeat.lengths)
        smoothed.likelihoods <- setNames(smoothed.likelihoods, obs.repeat.lengths)

        ## smoothed.likelihoods  <- predict(likelihoods.loess, possible.repeat.lengths)
        ## smoothed.likelihoods  <- setNames(smoothed.likelihoods, possible.repeat.lengths)
        ## implied.repeat.length <- as.numeric(names(which(smoothed.likelihoods == max(smoothed.likelihoods, na.rm=TRUE))))
        
        if( names(which(smoothed.likelihoods == max(smoothed.likelihoods, na.rm=TRUE))) == max(obs.repeat.lengths, na.rm=TRUE) ) {
          local.minima <- localMinima(smoothed.likelihoods)
          implied.repeat.length <- as.numeric(names(which(smoothed.likelihoods[1:local.minima[length(local.minima)]] == max(smoothed.likelihoods[1:local.minima[length(local.minima)]], na.rm=TRUE))))
        } else {
          implied.repeat.length <- as.numeric(names(which(smoothed.likelihoods == max(smoothed.likelihoods, na.rm=TRUE))))
        }

        results <- data.frame(repeat.id, ref.repeat.length, implied.repeat.length, total.n.flanking)
        writeResults(results, fh)
      }
    }
  }
  close(fh)
}

localMinima <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}

calculate.IMD.prob.distribution <- function(imd.distribution) {
  total        <- sum(imd.distribution$n)
  probs        <- imd.distribution$n / total
  imd.probs    <- setNames(probs, imd.distribution$inner.mate.distance)
  return(imd.probs)
}  

IMD.mass.sum <- function(repeat.length, imd.prob.distribution, margin, d95) {
  candidate.imds <- (repeat.length - 2 * margin) : (repeat.length - 2 * margin + 2 * d95)
  prob.obs.imd   <- imd.prob.distribution[as.character(candidate.imds)]
  prob.obs.imd[is.na(prob.obs.imd)] <- 0
  n.flanking.pos <- candidate.imds - repeat.length + 2 * margin + 1
  n.flanking.pos[n.flanking.pos < 0] <- 0
  n.pos.within.d95 <- 2 * d95 - candidate.imds + repeat.length - 2 * margin + 1
  n.pos.within.d95[n.pos.within.d95 < 0] <- 0
  n.positions <- pmin(n.flanking.pos, n.pos.within.d95)
  masses <- prob.obs.imd * n.positions
  C.replen <- sum(masses)
  return(C.replen)
}

IMD.mass.sums <- function(possible.repeat.lengths, imd.prob.distribution, margin, d95)  {
  names(possible.repeat.lengths) <- possible.repeat.lengths
  all.mass.sums <- sapply(possible.repeat.lengths, IMD.mass.sum,
                      imd.prob.distribution = imd.prob.distribution,
                      margin = margin,
                      d95 = d95)
  return(all.mass.sums)
}

IMD.logLikelihoodOffset <- function(flanking.imd.counts, log.imd.prob.distribution, log.mass.sums, ref.repeat.length)  {
  possible.repeat.lengths <- as.numeric(names(log.mass.sums))
  flanking.imds           <- as.numeric(names(flanking.imd.counts))

  logLikelihoodOffset <- sapply(possible.repeat.lengths, function(repeat.length)  {
    replen.string     <- as.character(repeat.length)
    implied.true.imds <- flanking.imds + repeat.length - ref.repeat.length

    first.term  <- -flanking.imd.counts * log.mass.sums[replen.string]
    second.term <- flanking.imd.counts * log.imd.prob.distribution[as.character(implied.true.imds)]

    results <- first.term + second.term
    sum(results)
  })
  logLikelihoodOffset <- setNames(logLikelihoodOffset, possible.repeat.lengths)
  return(logLikelihoodOffset)
}
