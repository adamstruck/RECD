\name{mapATR}
\alias{mapATR}
\title{mapATR [mapAlignmentsToRepeats]}
\usage{
  mapATR(con, table_name, input_file, imd_dist = NULL,
    upper_limit = 0.95, margin = (read_length/4),
    read_length, min_repeat_length = 12, output_dir = ".",
    append_chr = FALSE, overwrite = TRUE)
}
\arguments{
  \item{con}{database connection}

  \item{table_name}{name of table in database}

  \item{input_file}{inner_mate_range.tsv file (make sure
  chr names are compatible with your database)}

  \item{imd_dist}{.insert_length_dist file (i.e. the inner
  mate distance distribution file for the sample)}

  \item{upper_limit}{values between 0 and 1 will be used to
  set prob in quantile() and then used to calculate the
  maximum allowed distance from the repeat. Values greater
  than 1 will be applied as the maximum. Read pairs whose
  ends are farther from the repeat than the value set by
  this parameter will be not be used. [default= 0.95]}

  \item{margin}{minimum bp a read must be anchored outside
  the repeat [default= read_length/4]}

  \item{read_length}{read length}

  \item{min_repeat_length}{minimum repeat length in
  reference genome [default= 12]}

  \item{output_dir}{output directory to write results
  [default= "."]}

  \item{append_chr}{appends "chr" to chromosome names of
  inner_mate_ranges file if TRUE [default= TRUE]}

  \item{overwrite}{should output files be overwritten?
  [default= TRUE]}

  \item{margin}{maximum bp overlap between each read and a
  given repeat repeat [default= read_length/2]}
}
\value{
  _flanking, _facing, and _fully_repetitive .tsv files
}
\description{
  Sorts read pairs into three categories: flanking, facing,
  and fully repetitive.
}
\details{
  Description of output files: flanking: read pair
  completely flank the repeat or are anchored outside the
  repeat (minimum = margin) facing: one read flanks or is
  anchored outside the repeat (minimum = margin) while its
  pair aligned inside the repeat (or didn't meet minimum
  anchor requirements) fully_repetitive: both reads in pair
  aligned inside the repeat (or didn't meet minimum anchor
  requirements)
}
\author{
  Adam Struck
}

