#' Function to shingle (token or gram) a string into its k components
#'
#' @param record String or record
#' @param k Parameter k, which is the number of shingle, tokens, or grams to break the string into
#' @return Computes the shingled (tokened or grammed) version of a string
#' @export
#' @examples
#' shingles("Alexander",2)
#' shingles("Alexander Smith", 2)

shingles <- function(record,k){
	# factors only convert to characters properly elementwise,
	# not by applying as.character() to a whole row of a data
	# frame
	char_record <- lapply(record,as.character)
	string <- paste(char_record, collapse=" ")
	k_substring <- function(start){
		substring(string,start,start + k - 1)
	}
	tokens <- lapply(X=seq(1,nchar(string)-k+1), k_substring)
	return(tokens)
}

#' Function to convert to tell what index the shingle corresponds to in the record
#'
#' @param shingled_record Shingled record
#' @param universal_set Universal set of all shingles
#' @return the index regarding where the shingle falls in the record
#' @export
#' @examples
#' shingles("Alexander",2)
#' shingles("Alexander Smith", 2)
#' shingled_record_to_index_vec(shingles("Alexander",2), unique(shingles("Alexander Smith", 2)))

shingled_record_to_index_vec <- function(shingled_record, universal_set) {
	unique_tokens_in_record <- unique(shingled_record)
	token_indices <- match(unique_tokens_in_record, universal_set)
	return(token_indices)
	# TODO: Check if next line ever needed?
	#return(1:length(universal_set) %in% token_indices)
}

#' Function to create a matrix of minhashed signatures
#'
#' @import blink
#' @import plyr
#' @param shingled_records Shingled records
#' @param p Number of permutations to be applied to the hash function
#' @param do_one_hash_and_record Combination of one hash and one record
#' @return Computes an integer-valued matrix of minhash signatures with one row per permutation and one column per record
#' @export
#' @examples
#' head(data <- RLdata500[-c(2,4)])
#' minidata <- data[1:2,]
#' head(all_the_shingles <- apply(minidata,1,shingles,k=8))
#' head(minhash.minidata <- minhash_v2(all_the_shingles, p=10))

minhash_v2 <- function(shingled_records, p, do_one_hash_and_record=do_one_hash_and_record) {
	n.records <- length(shingled_records)

	# Figure out the universal set of all tokens
	print("Creating the universal set of tokens")
	print(system.time(universal_set_tokens <- unique(do.call(c, shingled_records)))[3])
	n.shingles <- length(universal_set_tokens)
	print("Number of tokens in universal set")
	print(n.shingles)

	# Generate a vector of p random hash functions
	print("Generating a vector of random hash functions")
	print(system.time(vector_of_hash_funcs <- rhash_funcs(n=p, size=n.shingles, vector.valued=FALSE))[3])

	# prepare a function to do the combination of one hash and one record
	# Presumes: rec_col is a vector saying which shingles (from the universal set) # are present in the shingled record
	do_one_hash_and_record <- function(h, rec_col) {
		# if (!timing) {
		# 	updated_v <- vector_of_hash_funcs[[h]](rec_col)
		# } else {
			applying_hash <- (updated_v <- vector_of_hash_funcs[[h]](rec_col))[3]
		# }
		# if (!timing) {
		# 	min_value <- min(updated_v)
		# } else {
			taking_min <- (min_value <- min(updated_v))[3]
		# }
		# if (timing) {
		# 	print("Applying hash function")
		# 	print(applying_hash)
		# 	print("Getting values")
		# 	print(getting_values)
		# 	print("Taking the minimum")
		# 	print(taking_min)
		#}
		return(min_value)
	}
	# Create a function to apply all the hash functions to one record
	 # This function must turn the shingled record into an indicator vector
	 # then apply all the functions
	  # then discard the indicator vector
	multiple_hash_one_record <- function(record) {
		index_vec <- shingled_record_to_index_vec(record, universal_set_tokens)
		multi_hash <- sapply(1:p, do_one_hash_and_record, rec_col=index_vec)
		return(multi_hash)
	}
	# Timing applying all the hash functions at once to one record
	print("Creating index vector and applying hash functions to first record")
	print(system.time(multiple_hash_one_record(shingled_records[[1]])))
	# Apply that multi-hash-function to all records
	signatures <- sapply(shingled_records, multiple_hash_one_record)
	# Return the matrix of minhash signatures
	return(signatures)
}

#' Function to generate all primes larger than an integer n1 (lower limit) and less than any other integer n2 (upper limit)
#'
#' @param n1 An integer taken to be 1 as the default
#' @param n2 Any integer n2
#' @return Generates all prime numbers with the above constraints
#' @export
#' @examples
#' primest(1, 5)
#' primest(1, 17)
primest <- function(n1=1, n2){
    p <- (n1+1):n2
    i <- 1
    while (p[i] <= sqrt(n2)) {
        p <-  p[p %% p[i] != 0 | p==p[i]]
        i <- i+1
    }
    p
}

#' Function to generate a vector of random hash functions (or optionally one vector-valued function)
#'
#' @param n Number of random hash functions
#' @param size Range of each size
#' @param vector.valued Flag for outputing vector of functions or vector-valued function
#' @param perfect Flag for whether a perfect permutation should be done, or just a hash function
#' @return Vector of n hash functions or a function which will take a number and return a vector of n different hashes of it
#' @export
#' @examples
#' rhash_funcs(1, 1, vector.valued=FALSE, perfect=FALSE)
#' rhash_funcs(5, 1, vector.valued=FALSE, perfect=FALSE)

# TODO: replace this with digest
rhash_funcs <- function(n, size, vector.valued, perfect=FALSE) {
	# Determine a suitable prime greater than size and =< 2*size
	candidate_primes <- primest(size,2*size)
	# Take the first suitable prime for simplicity's sake
	the_prime <- candidate_primes[1]
	# Create a single random hash function and return it
	if (!perfect) {
		# Make up a random hash function by modulo arithmetic
		make_one_hash_func <- function() {
			# Generate a function of the form ((ax+b) mod the_prime ) mod size
			# a,b < the_prime, a non-zero
			a <- as.integer64(sample(1:(the_prime-1),size=1))
			b <- as.integer64(sample(0:(the_prime-1),size=1))
			# Cast to a 64-bit integer
			the_prime <- as.integer64(the_prime)
			size <- as.integer64(size)
			hash_func <- function(x) {
				x <- as.integer64(x)
				h <- ((a*x+b) %% the_prime) %% size
				return(as.integer(h))
			}
			return(hash_func)
		}
	} else {
		# Make a perfect hash function that permutes the whole domain
		make_one_hash_func <- function() {
				perm <- sample(size)
				hash_func <- function(x) { perm[x] }
				return(hash_func)
		}
	}
	# Make a list of n random hash functions
	hash_func_list <- replicate(n, make_one_hash_func())
	if (vector.valued) {
		# Create a function which takes a number and returns a vector,
		# each component a different hash function's evaluation
  		   # TODO: replace iteration with something more vectorized
		   # want many functions with one input, not many inputs
		   # to one function
		vector_hash_func <- function(x) {
			h <- vector(length=n)
			for (i in 1:length(h)) {
				h[i] <- (hash_func_list[[i]])(x)
			}
			return(h)
		}
		# return the vector-valued hash function
		return(vector_hash_func)
	} else {
		return(hash_func_list)
	}
}


#' Function to take a signature matrix M composed of b bands and r rows and return
#' a bucket for each band for each record
#'
#' @param signature Signature matrix M composed of b bands and r rows
#' @param b Number of bands
#' @return Bucket for each band for each record
#' @export
#' @examples
#' head(data <- RLdata500[-c(2,4)])
#' minidata <- data[1:2,]
#' head(all_the_shingles <- apply(minidata,1,shingles,k=8))
#' head(minhash.minidata <- minhash_v2(all_the_shingles, p=10))
#' hash_signature(minhash.minidata, b=2)
#' hash_signature(minhash.minidata, b=5)

#assumes that signature matrix has been computed
	#take signature matrix M into b bands of r rows
	#returns the bucket for each band for each record
hash_signature <- function(signature,b){
	# need to divide signature into b bands of rows
	# so r = nrow(signature)/b, rounded down
	r = floor(nrow(signature)/b)
	extract_band <- function(i) { signature[((i-1)*r+1):(i*r),] }
	bands <- lapply(1:b,extract_band)
	# for each band, hash each portion of its columns to a hash table with k buckets
		#make k as large as possible
	band_hash <- sapply(bands,my_hash)
	# sapply builds up its output columnwise, so transpose to keep records as columns
	return(t(band_hash))
}

#' Function that applies a hash function to each column of the band from the
#' signature matrix
#' import bit64
#'
#' @import bit64
#' @param a_band Band from the signature matrix M
#' @return a 64 bit integer
#' @export
#' @examples
#' band1 <- c(2,1,2,1,2)
#' band2 <- c(4,5,2,1,9)
#' combined_band <- rbind(band1,band2)
#' my_hash(combined_band)

my_hash <- function(a_band) {
	hash64 <- function(x) {
		# if x is a vector, concatenate all its digits into one long number
		# use paste() to concatenate, then as.integer64 to turn into a big number
		y <- as.integer64(paste(x,sep="",collapse=""))
		return(hashfun(y,hashbits=64))
	}
	return(apply(a_band,2,hash64))
}


#' Function that extracts pairs of records from a band in the signature matrix M
#' import bit64
#'
#' @param a_band Band of the signature matrix M
#' @return The edgelist of record pairs that are connected
#' @export
#' @examples
#' band1 <- c(2,1,2,1,2)
#' extract_pairs_from_band(band1)
#' band2 <- c(6,7,8,9,6)
#' extract_pairs_from_band(band2)
#' band.12 <- rbind(band1, band2)
#' apply(band.12,1,extract_pairs_from_band)

extract_pairs_from_band <- function(a_band) {
   	# Each record has been mapped to some bucket within this band
   	# We now want to note down which pairs of records got mapped to the _same_
   	# bucket in this band (not caring about whether they got put in the same
   	# bucket in other bands)
   	record_pairs_in_bucket <- function(a_bucket) {
   		# print(paste("In bucket",a_bucket))
   		records_in_the_bucket <- which(a_band==a_bucket)
   		# print(paste(length(records_in_the_bucket),"records"))
   		if (length(records_in_the_bucket) > 1) {
	   		recpairs <- as.matrix(combn(records_in_the_bucket,m=2))
	   	} else {
	   		recpairs <- matrix(records_in_the_bucket,nrow=2,ncol=1)
	   	}
   		return(recpairs)
   	}
   	# Which buckets did we actually see in this band?
   	observed_buckets <- unique(a_band)
   	# Extract common pairs of records for each observed bucket
   	edgelist <- lapply(observed_buckets, record_pairs_in_bucket)
   	# We'll get back a list so bind them together columnwise
   	edgelist <- do.call(cbind,edgelist)
   	rownames(edgelist) <- c("rec1","rec2")
   	# We want the transpose
   	return(t(edgelist))
}

#' Function that creates a similarity graph and divides it into communities (or blocks) for entity resolution
#'
#' @import igraph
#' @import plyr
#' @param hashed_signatures The hashed signatures
#' @param max_bucket_size The largest block size allowed by user
#' @return max_bucket_size The largest bucket size (or block size) that one
#' can handle
#' @export
#' @examples
#' head(data <- RLdata500[-c(2,4)])
#' minidata <- data[1:2,]
#' head(all_the_shingles <- apply(minidata,1,shingles,k=8))
#' head(minhash.minidata <- minhash_v2(all_the_shingles, p=10))
#' hashed_signature <- hash_signature(minhash.minidata, b=5)
#' compare_buckets(hashed_signature, max_bucket_size=200)

# Create a similarity graph and divide it into communities, as blocks for record linkage
# Inputs: minhashed signature matrix, maximum block size
# Presumes: the signature matrix has been created by minhashing (or something like
  # it), so 2 records matching in some row indicates non-trivial similarity
compare_buckets <- function(hashed_signatures, max_bucket_size=1000) {
	# Create blocks from the buckets the bands were mapped to
	# General idea: each record gets put in multiple buckets from the multiple minhashes
	# Two records are a "candidate pair" if they get mapped to the same bucket
	# by some minhash or other
	# Form a graph, with records as nodes, and edges between candidate pairs
	# Divide the graph into dense sub-graphs (communities), subject to a maximum
	# size limit

    # Each row of hashed_signatures represents the bucket-mapping of
    # the records for a different minhash permutation
    # Apply extract_pairs_from_band to each row, and then combine the resulting
    # matrices of candidate-pair records into one big edgelist

    # TODO: Try using plyr rather than apply + do.call to see about speed
    print("Creating edgelist")
    edgelisting <- system.time(edgelist <- as.matrix(do.call(rbind,apply(hashed_signatures,1,extract_pairs_from_band))), gcFirst=FALSE)
    print(edgelisting)
    print(dim(edgelist))
	# Actually build the graph
	    # edgelist contains only edges in one direction, so we need to tell igraph
	    # that edges are directionless
	print("Building graph from edgelist")
	graphing <- system.time(candidate_pairs_graph <- graph.edgelist(edgelist, directed=FALSE), gcFirst=FALSE) # Actually build the graph
	# edgelist isn't needed any more and can be quite big, so remove it from memory
	rm(edgelist)
	print(graphing)
	# Remove multiple and self edges, if they exist
	candidate_pairs_graph <- simplify(candidate_pairs_graph)

	# Try dividing the graph into communities. Use a hierarchical community method
	# so that if the initial cut has communities which are too big, we can go further down
	# until they are small enough to work with.
	print("Dividing graph into communities initially")
	communitying <- system.time(initial_community <- fastgreedy.community(candidate_pairs_graph), gcFirst=FALSE)
	print(communitying)
	#save(candidate_pairs_graph,file = "candidate_pairs_graph.Rdata")


	# The graph has served its purpose and should go away
	rm(candidate_pairs_graph)

	# Sub-divide communities if too big
	max_comm_size <- max(sizes(initial_community))
	comm_number <- length(initial_community)
	comm_membership <- membership(initial_community)
	print("Subdividng communities")
	subdividing <- system.time(
	while(max_comm_size > max_bucket_size) {
		comm_number <- comm_number+1
		comm_membership <- cutat(initial_community, no=comm_number)
		max_comm_size <- max(table(comm_membership))
	}
	,gcFirst=FALSE)
	print(subdividing)
	blocks_members <- comm_membership
	num_blocks <- comm_number
	#save(blocks_members, file="blocks_members.Rdata")

	# Now create a list, saying which records are in which block
	records_per_block <- function(b) { which(blocks_members == b)}
	blocks <- lapply(1:num_blocks,records_per_block)
	return(blocks)
}
