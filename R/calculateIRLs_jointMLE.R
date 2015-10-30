##' For repeats that have read pairs flanking them, the implied expansion and length of that repeat is determined. 
##'
##' Summary stats for the implied repeat exansions/lengths are also reported. 
##' @title calculateIRLs_MLE
##' @param flanking_file _flanking read pair file
##' @param facing_file  _facing reads file
##' @param imd_distribution .insert_len_dist file (i.e. the inner mate distance distribution for the sample)
##' @param upper_limit the value used in mapATR % processLastAlignedBase [default= 0.95]
##' @param margin the value used in mapATR % processLastAlignedBase [default= 37.5]
##' @param output_dir path to output directory
##' @return _implied_repeat_lengths_jointMLE file
##' @author Adam Struck - Intern
##' @export
calculateIRLs_jointMLE <- function(flanking_file, facing_file, imd_distribution, upper_limit=0.95, margin=37.5, output_dir) {
    if( file.exists(output_dir) ) {
        setwd(output_dir)
    } else {
        dir.create(output_dir)
        setwd(output_dir)
    }
 
    irl.file.header <- c("repeat.id", "ref.repeat.length", "mean.imd", "implied.repeat.length", "n.flanking", "n.facing")
    fh <- openOutputFileHandle(flanking_file, suffix="implied_repeat_lengths_jointMLE", header=irl.file.header, overwrite=TRUE)
    
    imd.distribution <- read.table(imd_distribution, header=TRUE, sep="\t")
    names(imd.distribution) <- c("inner.mate.distance", "n")
    distribution <- rep(imd.distribution$inner.mate.distance, imd.distribution$n)

    if( upper_limit <= 1 ) {
        d95 <- quantile(distribution, prob=upper_limit)
    } else {
        d95 <- upper_limit
    }
    
    imd.distribution          <- imd.distribution[imd.distribution$inner.mate.distance <= 10000,]
    imd.prob.distribution     <- calculate.IMD.prob.distribution(imd.distribution)
    log.imd.prob.distribution <- log(imd.prob.distribution)

    possible.repeat.lengths <- 0 : 1000
    flanking.mass.sums <- get.flanking.mass.sums(possible.repeat.lengths, imd.prob.distribution, margin, d95)
    log.flanking.mass.sums <- log(flanking.mass.sums)

    facing.probabilities <- get.facing.probabilities(possible.repeat.lengths, imd.prob.distribution, margin, d95)
    log.facing.probabilities <- log(facing.probabilities)

    flanking.probabilities <- get.flanking.probabilities(facing.probabilities)
    log.flanking.probabilities <- log(flanking.probabilities)

    flanking <- read.table(flanking_file, header=TRUE, sep="\t")  
    facing <- read.table(facing_file, header=TRUE, sep="\t")

    repeats  <- unique(c(as.character(flanking$repeat.id), as.character(facing$repeat.id)))
  
    if(length(repeats) > 0) {
        for(i.rep in 1:length(repeats))  {
            if ((i.rep %% 50 == 0) == TRUE) { message("Working on repeat: ", i.rep , "/", length(repeats), " ", Sys.time()) }

            repeat.id         <- repeats[i.rep]
            rep.info          <- unlist(strsplit(as.character(repeat.id), "[:-]"))
            ref.repeat.length <- as.numeric(rep.info[3]) - as.numeric(rep.info[2]) + 1

            flanking.subset <- flanking[flanking$repeat.id == repeat.id,]
            flanking.counts <- setNames(flanking.subset$n, flanking.subset$inner.mate.distance)
            mean.imd <- mean(flanking.subset$inner.mate.distance, na.rm=TRUE)
            total.n.flanking <- sum(flanking.counts)

            facing.subset <- facing[facing$repeat.id == repeat.id,]
            total.n.facing <- sum(facing.subset$n)

            likelihoods <- logLikelihoodOffset(flanking.counts, total.n.facing, log.imd.prob.distribution, log.flanking.mass.sums, log.flanking.probabilities, log.facing.probabilities, ref.repeat.length) 
            implied.repeat.length <- as.numeric(names(which(likelihoods == max(likelihoods, na.rm=TRUE))))
            
            results <- data.frame(repeat.id, ref.repeat.length, mean.imd, implied.repeat.length, total.n.flanking, total.n.facing)
            writeResults(results, fh)
        }
    }
    close(fh)
}


calculate.IMD.prob.distribution <- function(imd.distribution) {
    total        <- sum(imd.distribution$n)
    probs        <- imd.distribution$n / total
    imd.probs    <- setNames(probs, imd.distribution$inner.mate.distance)
    return(imd.probs)
}  

n.facing <- function(repeat.length, margin, d95) {
    candidate.imds <- -margin : (repeat.length - 2 * margin + d95)
    candidate.imds <- round(candidate.imds, 0)
    less.or.equal.end <- pmax( 1 , (candidate.imds - d95 + margin) )
    greater.or.equal.end <- pmin( (repeat.length - margin), (candidate.imds + margin - 1) )
    n.positions <- 2 * pmax( 0, (greater.or.equal.end - less.or.equal.end + 1) )
    n.positions <- setNames(n.positions, candidate.imds)
    return(n.positions)
}

n.flanking <- function(repeat.length, margin, d95) {
    candidate.imds <- (repeat.length - 2 * margin) : (repeat.length - 2 * margin + 2 * d95)
    n.flanking.pos <- candidate.imds - repeat.length + 2 * margin + 1
    n.pos.within.d95 <- 2 * d95 - candidate.imds + repeat.length - 2 * margin + 1
    n.positions <- pmin(n.flanking.pos, n.pos.within.d95)
    n.positions <- setNames(n.positions, candidate.imds)
    return(n.positions)
}

get.flanking.mass.sum <- function(repeat.length, imd.prob.distribution, margin, d95) {
    candidate.imds <- (repeat.length - 2 * margin) : (repeat.length - 2 * margin + 2 * d95)
    prob.obs.imd   <- imd.prob.distribution[as.character(candidate.imds)]
    prob.obs.imd[is.na(prob.obs.imd)] <- 0

    n.positions <- n.flanking(repeat.length, margin, d95)
    masses <- prob.obs.imd * n.positions
    C.replen <- sum(masses)

    return(C.replen)
}

get.flanking.mass.sums <- function(possible.repeat.lengths, imd.prob.distribution, margin, d95)  {
    names(possible.repeat.lengths) <- possible.repeat.lengths
    flanking.mass.sums <- sapply(possible.repeat.lengths, get.flanking.mass.sum,
                                 imd.prob.distribution = imd.prob.distribution,
                                 margin = margin,
                                 d95 = d95)
    return(flanking.mass.sums)
}


get.facing.probability <- function(repeat.length, imd.prob.distribution, margin, d95) {
    candidate.imds <- -margin : (repeat.length - 2 * margin + 2 * d95)
    candidate.imds <- round(candidate.imds, 0)
    prob.obs.imd   <- imd.prob.distribution[as.character(candidate.imds)]
    prob.obs.imd[is.na(prob.obs.imd)] <- 0
    
    n.flanking.positions <- n.flanking(repeat.length, margin, d95)
    n.facing.positions <- n.facing(repeat.length, margin, d95)

    facing.pos <- n.facing.positions[as.character(candidate.imds)]
    facing.pos[is.na(facing.pos)] <- 0

    flanking.pos <- n.flanking.positions[as.character(candidate.imds)]
    flanking.pos[is.na(flanking.pos)] <- 0
  
    facing.mass <- sum( facing.pos * prob.obs.imd)
    combined.mass <- sum( (facing.pos + flanking.pos) * prob.obs.imd )
    p.facing <- facing.mass / combined.mass

    return(p.facing)
}

get.facing.probabilities <- function(possible.repeat.lengths, imd.prob.distribution, margin, d95)  {
    names(possible.repeat.lengths) <- possible.repeat.lengths
    all.facing.probs <- sapply(possible.repeat.lengths, get.facing.probability,
                               imd.prob.distribution = imd.prob.distribution,
                               margin = margin,
                               d95 = d95)
  
    return(all.facing.probs)
}

get.flanking.probabilities <- function(facing.probabilities) {
    1 - facing.probabilities
}

logLikelihoodOffset <- function(flanking.counts, total.n.facing, log.imd.prob.distribution, log.flanking.mass.sums, log.flanking.probabilities, log.facing.probabilities, ref.repeat.length)  {
    possible.repeat.lengths <- as.numeric(names(log.flanking.mass.sums))
    flanking.imds           <- as.numeric(names(flanking.counts))
    
    if( length(flanking.counts) == 0) {
        flanking.counts <- 0
    }
  
    logLikelihoodOffset <- sapply(possible.repeat.lengths, function(repeat.length)  {
        replen.string     <- as.character(repeat.length)
        implied.true.imds <- flanking.imds + repeat.length - ref.repeat.length
        
        first.term  <- -flanking.counts * log.flanking.mass.sums[replen.string]
        second.term <- flanking.counts * log.imd.prob.distribution[as.character(implied.true.imds)]
        third.term  <- flanking.counts * log.flanking.probabilities[replen.string]
        fourth.term <- total.n.facing * log.facing.probabilities[replen.string]

        if( total.n.facing > 0 ) {
            results <- sum(first.term + second.term + third.term) + fourth.term
            results
        } else {
            results <- sum(first.term + second.term + third.term)
            results
        }
    })
  
    logLikelihoodOffset <- setNames(logLikelihoodOffset, possible.repeat.lengths)
    return(logLikelihoodOffset)
}
