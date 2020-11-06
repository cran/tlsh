# This is the one of the main blocking methods in Steorts, Ventura, Sadinle,
# Fienberg (2014), Privacy in Statistical Databases.
# If you use this code, please cite Steorts, R., Ventura, S., Sadinle, M., and
# Fienberg, S. (2014). "Blocking Comparisons for Record Linkage." Privacy
# in Statistical Databases (Lecture Notes in Computer Science 8744), ed. J
# Domingo-Ferrer, Springer, 252-268, doi:10.1007/978-3-319-11257-220.

# tlsh Copyright 2018 Rebecca C. Steorts (beka@stat.duke.edu)

# tlsh is free software: you can redistribute it and/or modify it
# under the terms of the Creative Commons license, either version 3 of the
# license, or (at your option) any later version.

# tlsh is distributed in the hope it will be useful, but without ANY WARRANTY;
# without even the implied warranty of merchantability or fitness for a particular
# purpose. Specifically, you may share the software in any medium or format and
# you may adapt the software. Credit must be given when either of these are
# given to indicate if and what changes were made. The software may not be
# used for noncommerical purposes. If you are interested in using the software
#  for commercial purposes, please contact the author above.
###########################################################################################################


#Begin working example
# TODO: make sure the blocks are saved.

#library(plyr)
#library(digest)
#library(RecordLinkage)
#data(RLdata500)
#minidata <- RLdata500[-c(2,4)]
#The command
#rl_data_500_b26 <- adply(1:5, .margins=1, .fun = eval.blocksetup,  dat=minidata, b=26, .expand=F,key=identity.RLdata500)
#plot(1:5, rl_data_500_b26[,2],xlab="k",ylab="Recall")
#plot(1:5, rl_data_500_b26[,3],xlab="k",ylab="Elaped Time")


# will loop through shingles 1:5 and save the recall and the runtime. We should also
# save the precision and reduction ratio as well.

#rl_data_500_b22_30 <- adply(2:8, .margins=1, .fun = eval.blocksetup, dat = #RLdata500, b=22, .expand=F,key=identity.RLdata500)
#save(rl_data_500_b22_30, file="rl_data_500_b22_10.Rdata")

# plot(2:8, rl_data_500_b22_30[,2],xlab="k",ylab="Recall",ylim=c(0,0.95),type="b")
# points(2:8, rl_data_500_b22_30[,2], xlab="k",ylab="Recall",ylim=c(0,0.1),pch=2,type="b")
# points(2:8, rl_data_500_b22_50[,2], xlab="k",ylab="Recall",ylim=c(0,0.1),pch=3,type="b")
# legend("bottomright", legend= c("10%", "30%","50%"), pch=c(1,2,3))

#End working example

# ATTN: There are additional functions below that will allow TLSH
# to be integrated into random forests with a mapping function for
# parallezation.


#' Function to evaluate the blocking step
#'
#' import blink
#' @param dat Data set
#' @param b Number of buckets
#' @param k Parameter k, which is the number of shingle, tokens, or grams to break the string into
#' @param key Unique identifier
#' @return Recall and runtime
#' @export
#' @examples
#' r.set <- RLdata500[1:50,c(-2)]
#' eval.blocksetup(r.set, k=2, b=22, key=identity.RLdata500)

eval.blocksetup <- function(dat, k=5, b=21, key){
	#runtime <- as.numeric((mapping <- block_setup_v2(dat, b=b, k=k))[3] )
  mapping <- block_setup_v2(dat, b=b, k=k)
	recall<- confusion.from.blocking (blocking=mapping,true_ids=key,recall.only=TRUE)[[1]]
	return(data.frame(recall))
}



#' Function that divides all records into bins using locality sensitive hashing and using TLSH (based upon community detection technique)
#'
#' import blink
#' @param r.set Record set (shingled records)
#' @param b Band
#' @param save_signature Flag of whether or not to save the signature
#' @param k Shingle size
#' @return List of blocks where a particular index is the record id in the original
#' data set
#' @export
#' @examples
#' r.set <- RLdata500[1:3,c(-2)]
#' block_setup_v2(r.set = RLdata500[1:3,c(-2)], b=22, save_signature=FALSE, k=2)

block_setup_v2 <- function(r.set, b=22, save_signature=FALSE,  k=5) {
	# for each record r in r.set
	  # calculate the hash function of the record r, say h
	  # store r under h in the hash map
	# return hash map from hash values to sets of records

	# Convert each record (= row of r.set) to k-token shingles
	shingled_records <- apply(r.set,1,shingles,k=k)
	# Create the matrix of minhashed signatures, using p random permutations
	# ATTN: Put this in parallel and test that it works

	minhash_time <- system.time(minhashed_records <- minhash_v2(shingled_records,p=100),gcFirst=FALSE)
	print(minhash_time)
	if(save_signature) {
		timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
		save(minhashed_records, file=paste("minhashed_signature", timestamp))
	}

	# Get rid of the shingled records as they've served their purpose
	rm(shingled_records)

	# Calculate signatures, put into buckets, make the graph, return blocks
	return(compare_buckets(hash_signature(minhashed_records,b=b)))
}
