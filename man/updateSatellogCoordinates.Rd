\name{updateSatellogCoordinates}
\alias{updateSatellogCoordinates}
\title{updateSatellogCoordinates}
\usage{
  updateSatellogCoordinates(con, con2, exe.path = "./")
}
\arguments{
  \item{con}{MySQL database connection}

  \item{con2}{SQLite database connection}

  \item{exe.path}{location of the liftOver executable and
  chain files}
}
\value{
  SQLite database
}
\description{
  Updates all of the genomic coordinates in the repeats
  table of your Satellog database to the current hg19
  build.
}
\details{
  Since Satellog was developed using hg16 v34 of the human
  genome, this package updates the coordinates to the most
  recent hg19 build using the liftOver tool from UCSC.
  Updated coordinates are inserted into SQLite database
}
\author{
  Adam Struck
}

