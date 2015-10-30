\name{RepeatsInRegion}
\alias{RepeatsInRegion}
\title{RepeatsInRegion}
\usage{
  RepeatsInRegion(con, chr, start, end)
}
\arguments{
  \item{con}{database connection}

  \item{chr}{chromosome}

  \item{start}{starting genomic coordinates}

  \item{end}{ending genomic coordinates}
}
\value{
  data.frame
}
\description{
  Connect to a MySQL or SQLite database containing SSRs and
  selelct entries within a regionon a chromosome.
}
\details{
  .. content for \details{} ..
}
\author{
  Adam Struck
}

