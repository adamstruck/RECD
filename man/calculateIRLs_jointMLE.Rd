\name{calculateIRLs_jointMLE}
\alias{calculateIRLs_jointMLE}
\title{calculateIRLs_MLE}
\usage{
  calculateIRLs_jointMLE(flanking_file, facing_file,
    imd_distribution, upper_limit = 0.95, margin = 37.5,
    output_dir)
}
\arguments{
  \item{flanking_file}{_flanking read pair file}

  \item{facing_file}{_facing reads file}

  \item{imd_distribution}{.insert_len_dist file (i.e. the
  inner mate distance distribution for the sample)}

  \item{upper_limit}{the value used in mapATR %
  processLastAlignedBase [default= 0.95]}

  \item{margin}{the value used in mapATR %
  processLastAlignedBase [default= 37.5]}

  \item{output_dir}{path to output directory}
}
\value{
  _implied_repeat_lengths_jointMLE file
}
\description{
  For repeats that have read pairs flanking them, the
  implied expansion and length of that repeat is
  determined.
}
\details{
  Summary stats for the implied repeat exansions/lengths
  are also reported.
}
\author{
  Adam Struck - Intern
}

