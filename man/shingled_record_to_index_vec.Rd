% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minhash_v2.R
\name{shingled_record_to_index_vec}
\alias{shingled_record_to_index_vec}
\title{Function to convert to tell what index the shingle corresponds to in the record}
\usage{
shingled_record_to_index_vec(shingled_record, universal_set)
}
\arguments{
\item{shingled_record}{Shingled record}

\item{universal_set}{Universal set of all shingles}
}
\value{
the index regarding where the shingle falls in the record
}
\description{
Function to convert to tell what index the shingle corresponds to in the record
}
\examples{
shingles("Alexander",2)
shingles("Alexander Smith", 2)
shingled_record_to_index_vec(shingles("Alexander",2), unique(shingles("Alexander Smith", 2)))
}
