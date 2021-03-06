% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minhash_v2.R
\name{compare_buckets}
\alias{compare_buckets}
\title{Function that creates a similarity graph and divides it into communities (or blocks) for entity resolution}
\usage{
compare_buckets(hashed_signatures, max_bucket_size = 1000)
}
\arguments{
\item{hashed_signatures}{The hashed signatures}

\item{max_bucket_size}{The largest block size allowed by user}
}
\value{
max_bucket_size The largest bucket size (or block size) that one
can handle
}
\description{
Function that creates a similarity graph and divides it into communities (or blocks) for entity resolution
}
\examples{
head(data <- RLdata500[-c(2,4)])
minidata <- data[1:2,]
head(all_the_shingles <- apply(minidata,1,shingles,k=8))
head(minhash.minidata <- minhash_v2(all_the_shingles, p=10))
hashed_signature <- hash_signature(minhash.minidata, b=5)
compare_buckets(hashed_signature, max_bucket_size=200)
}
