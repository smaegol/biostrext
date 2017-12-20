
#' Checks if provided object is of class DNAString
#'
#' @param x - R object
#'
#' @return True if object is DNAString, false otherwise
#' @export
#'
#' @examples
#' is_DNAstring(Biostrings::DNAString("AGCTAGCA"))
#' is_DNAstring("AGCTAGCA")
is_DNAstring <- function(x) {
  return(isTRUE(class(x)[1] == 'DNAString'))
}

#' Checks if provided object is of class DNAStringSet
#'
#' @param x - R object
#'
#' @return True if object is DNAString, false otherwise
#' @export
#'
#' @examples
#' is_DNAstringSet(Biostrings::DNAStringSet(c("AGCTAGCA","CACGATC","CGACGAC")))
#' is_DNAstringSet("AGCTAGCA")
is_DNAstringSet <- function(x) {
  return(isTRUE(class(x)[1]=='DNAStringSet'))
}

#' Calculates %GC content for a single DNAString object
#'
#' @param x DNAString object
#'
#' @return %GC content of provided sequence
#' @export
#'
#' @examples
#' calc_GC_single(Biostrings::DNAString("CGAGCACACATCACAC"))
calc_GC_single <- function(x) {
  #check if Biostrings package is loaded
  sequence_length_min = 4 # set the minimum sequence length requried for %GC calculation
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #check if parameters are provided
  if (missing(x)) {
    stop("The sequence object (argument x) is missing",
         call. = FALSE)
  }
  #check if input is of DNAString class
  assertthat::assert_that(is_DNAstring(x))
  #get length of provided DNAString object (sequence)
  seq_length = length(x)
  #calculate %GC if sequence length is > sequence_length_min (default = 4)
  if(seq_length>=sequence_length_min) {
    GC <- sum(Biostrings::alphabetFrequency(x)[c(2,3)])/seq_length*100
  }
  else{
    #if sequence is too short - return NA
    GC <- NA
  }
    return(GC)
}


#' Split a DNAString object in smaller chunks
#'
#' @param x DNAString object
#' @param window size of a sequence chunk
#' @param overlap size of the overlap between chunks (have to be < window, default = 0)
#'
#' @return DNAStringSet object
#' @export
#'
#' @examples
#' split_in_windows(Biostrings::DNAString("ACGTAGCATGCTAGCTACGTAGCTA"),4,1)
split_in_windows <- function(x,window,overlap = 0) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #check if parameters are provided
  if (missing(x)) {
    stop("The sequence object (argument x) is missing",
         call. = FALSE)
  }
  if (missing(window)) {
    stop("Window argument is missing",
         call. = FALSE)
  }
  #check if DNAString object is provided as input
  assertthat::assert_that(is_DNAstring(x))
  #make sure that window and overlap parameters are numeric
  assertthat::assert_that(is.numeric(window))
  assertthat::assert_that(is.numeric(overlap))
  #assure that window size is less (or equal) than the length of sequence
  assertthat::assert_that(window<=length(x))
  #assure that window size is a positive integer
  assertthat::assert_that(window>0)
  #make sure that overlap is not less than 0
  assertthat::assert_that(overlap>=0)
  #assure that overlap is not larger (or equal to) than window
  assertthat::assert_that(overlap<window)

  #get provided sequence length
  seq_length = length(x)
  #initialize vector for return values
  return_vector <- character(0)
  #initialize start position for the first chunk
  start=1
  #initialize end position for the first chunk
  end=window
  #chunk number iterator
  chunk_no = 1
  #while the end fo the chunk is not after the sequence end, keep processing
  while (end<=seq_length) {
    #create temporary chunk, using subseq function from Biostrings
    temp_seq<-as.character(Biostrings::subseq(x,start,end))
    #store created chunk in the return vector
    return_vector[chunk_no]<-temp_seq
    #increase iterators
    chunk_no=chunk_no+1
    start=start+window-overlap
    end=end+window-overlap

  }
  #if after the while loop start position is still within the sequence, create the last chunk
  if(start<seq_length) {
    temp_seq<-as.character((Biostrings::subseq(x,start,seq_length)))
    return_vector[chunk_no]<-temp_seq
  }
  #create the DNAStringSet object, which will be returned by function
  return_set<-Biostrings::DNAStringSet(return_vector)
  return(return_set)

}

#' Wrapper for calc_GC_single
#'
#' @param x DNAString or DNAStringSet object
#'
#' @return %GC for single sequence or vector of %GC if DNAStringSet was provided as input
#' @export
#'
#' @examples
#' calc_GC(Biostrings::DNAString("ACAGCAGTCA"))
#' calc_GC(Biostrings::DNAStringSet(c("ACAGCAGTCA","ACGTAGCATCGA","CGATCGACA")))
calc_GC <- function(x) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #check if parameters are provided
  if (missing(x)) {
    stop("The sequence object (argument x) is missing",
         call. = FALSE)
  }
  #make sure that DNAString of DNAStringSet is provided as an input
  assertthat::assert_that(is_DNAstring(x) || is_DNAstringSet(x))
  #if DNAString is an input simply call calc_GC_single
  if(is_DNAstring(x)) {
    GC = calc_GC_single((x))
  }
  #if an input is a DNAStringSet object, calculate GC for each element fro mthe set
  else if (is_DNAstringSet(x)) {
    #get the number of elements in DNAStringSet
    no_elements_in_set = length(x)
    #initialize numeric vector for return values
    GC <- numeric(0)
    #calculate %GC for each element in set
    for (z in seq(1,sno_elements_in_set)) {
      GC_temp <- calc_GC_single(x[[z]])
      GC<-c(GC,GC_temp)
    }
  }
  return(GC)
}

#' Split provided DNAString in chunks and calculate %GC for each chunk
#'
#' @param x DNAString object
#' @param window size of a sequence chunk
#' @param overlap size of the overlap between chunks (have to be < window, default =0)
#'
#' @return vector of calculated %GC content for each chunk generated
#' @export
#'
#' @examples
#' calc_windowed_GC(Biostrings::DNAString("CAGTCAGTCGATC"),10,2)
calc_windowed_GC <- function(x,window,overlap = 0) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #check if parameters are provided
  if (missing(x)) {
    stop("The sequence object (argument x) is missing",
         call. = FALSE)
  }
  if (missing(window)) {
    stop("Window argument is missing",
         call. = FALSE)
  }
  #make sure that DNAString object is provided as an input
  assertthat::assert_that(is_DNAstring(x))
  #split provided input in chunks
  splitted_seq <- split_in_windows(x,window,overlap)
  #calculate %GC for created chunks
  GC <- calc_GC(splitted_seq)
  return(GC)
}

#' Title Get kmer which is most frequent in provided sequence
#'
#' @param x DNAString object
#' @param k kmer length
#' @param force_large_kmer logical, should calcualtion of large kmers be possible (but shorter than 15)
#'
#' @return named numeric, with count of given kmer and its sequence as name
#' @export
#'
#' @examples
#' get_max_kmer_single(Biostrings::DNAString("CAGCTGCACACATCGTACA"),3)
get_max_kmer_single <- function(x,k,force_large_kmer = FALSE) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #check if parameters are provided
  if (missing(x)) {
    stop("The sequence object (argument x) is missing",
         call. = FALSE)
  }
  if (missing(k)) {
    stop("Kmer size(k) argument is missing",
         call. = FALSE)
  }
  #make sure that DNAString object is provided as an input
  assertthat::assert_that(is_DNAstring(x))
  #make sure that kmer size k is a numeric variable
  assertthat::assert_that(is.numeric(k))
  #make sure that kmer size k is a positive number
  assertthat::assert_that(k>0)
  #make sure that kmer size k is not larger than sequence length
  assertthat::assert_that(k<length(x))
  #if kmer size k is larger than 10 - stop calculation unless force_large_kmer is TRUE
  if (k>10 && !isTRUE(force_large_kmer)) {
    stop("Calculation of kmers larger than 10 is computationally intensive and may take a long time.
If you need that, specify force_large_kmer=TRUE (use with caution!)",
         call. = FALSE)
  }
  #if mer size is larger than 15 - stop calculation without any option of forcing calculation
  if (k>15) {
    stop("Calculation of kmers larger than 15 is extremally computationally intensive and it is discouraged to use.",
         call. = FALSE)
  }
  #calculate counts for each possible kmer of size k
  kmersCount <- Biostrings::oligonucleotideFrequency(x,k)
  #get most frequent kmer
  max_kmer_count <- max(kmersCount)
  #get sequence of the most frequent kmer
  max_kmer_seq <- names(which.max(kmersCount))
  #create return variable
  names(max_kmer_count) <- max_kmer_seq
  return(max_kmer_count)
}


#' Calculates maximum kmer for DNA sequence object (DNAString or DNAStringSet)
#'
#' @param x DNAString or DNAStringSet object
#' @param mer kmer size
#'
#' @return named vector of max kmer counts with kmer sequences as names
#' @export
#'
#' @examples
#' get_max_kmer(Biostrings::DNAString("CAGCGATCGTACACGATC"),4)
#' get_max_kmer(Biostrings::DNAStringSet(c("CAGCGATCGTACACGATC","CAGCGTAGCTAG","CAGCATCACACACAC")),3)
get_max_kmer <- function(x,k) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #check if parameters are provided
  if (missing(x)) {
    stop("The sequence object (argument x) is missing",
         call. = FALSE)
  }
  if (missing(k)) {
    stop("Kmer size(k) argument is missing",
         call. = FALSE)
  }
  #make sure an input is provided as DNAString or DNAStringSet object
  assertthat::assert_that(is_DNAstring(x) || is_DNAstringSet(x))
  #make sure that kmer size k is a numeric variable
  assertthat::assert_that(is.numeric(k))
  #make sure k is a positive number
  assertthat::assert_that(k>0)
  #if DNAString object is provided as an input - simply call get_max_kmer_single
  if(is_DNAstring(x)) {
    max_kmer = get_max_kmer_single(x,k)
  }
  #for DNAStringSet object - calculate get_max_kmer_single for each sequence in the set
  else if (is_DNAstringSet(x)) {
    no_elements_in_set = length(x)
    max_kmer <- numeric(0)
    for (z in seq(1,no_elements_in_set)) {
      #if a sequence is longer than kmer size k, get max kmer
      if (length(x[[z]])>k) {
        max_kmer_temp <- get_max_kmer_single(x[[z]],k)
      }
      else {
        #if kmer size k is larger than sequence length - return NA
        max_kmer_temp <- NA
      }
      max_kmer<-c(max_kmer,max_kmer_temp)
    }
  }
  return(max_kmer)
}

#' Plot %GC distribution across sequence in fragments of given size
#'
#' @param x DNAString object
#' @param window size of the sequence chunk
#' @param overlap size of the overlap between chunks
#'
#' @return ggplot2 geom_smooth plot with %GC distribution across sequence
#' @export
#'
#' @examples
#' plot_GC_distribution(Biostrings::DNAString("CAGCTAGCTAGCTAGGCAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTACGTACGTAGCTACACTAGCTAGCTAGCTAGCTAATAATTATGGCGCGC"),6,2)
plot_GC_distribution <- function(x,window,overlap) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #make sure an input is provided as DNAString object
  assertthat::assert_that(is_DNAstring(x))
  #make sure that input sequence will be divided in proper number of chunks (at least 20)
  assertthat::assert_that(length(x)/(window-overlap)>20)
  #calculate step size - how sliding window is moving
  step=window-overlap
  #calculate start position for each chunk based on calculated step size
  positions<-seq(1,length(x)-window+step+1,step)
  #calculate %GC for chunks
  GC <- calc_windowed_GC(x,window,overlap)
  #name each calculated value with the start position of chunk
  GC_frame <- data.frame(position=positions,GC=GC)
  #create plot using geom_smooth from ggplot2
  plot <- ggplot2::ggplot(GC_frame,ggplot2::aes(x=position,y=GC)) + ggplot2::geom_smooth(span=0.1) + ggplot2::scale_y_continuous() + ggplot2::ggtitle("%GC content distribution accross sequence")
  return(plot)
}


#' Plots a distribution of given sequence fragment across the whole sequence
#'
#' @param x DNAString object
#' @param kmer - string
#'
#' @return plot with distribution of given kmer across sequence
#' @export
#'
#' @examples
#' plot_kmer_distribution(Biostrings::DNAString("CAGCTAGCTAGCTAGGCAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTACGTACGTAGCTACACTAGCTAGCTAGCTAGCTAATAATTATGGCGCGC"),"TAC")
plot_kmer_distribution <- function(x,kmer) {
  #check if Biostrings package is loaded
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  #define limits for linearange geom
  min_plot_value = -0.1
  max_plot_value = 0.1
  #make sure an input is provided as DNAString object
  assertthat::assert_that(is_DNAstring(x))
  #make sure that provided kmer sequence is a character variable
  assertthat::assert_that(is.character(kmer))
  #make sure that kmer is not an empty string
  assertthat::assert_that(nchar(kmer)>0)
  #make sure that kmer variable contain only DNA nucleotide letters
  if (grepl("[^AGCT]",kmer)) {
    stop("Wrong character in the kmer sequence provided. Please use only {A,G,C,T} letters.",
         call. = FALSE)
  }
  #find all occurences of kmer in sequence x using matchPattern from Biostrings
  kmer_matches<-Biostrings::matchPattern(kmer,x)
  if (length(kmer_matches)>0) {
    #if any occurence was found, get position of each occurence using start function from Biostrings
    kmer_starts=Biostrings::start(kmer_matches)
    #create data frame with values
    positions_df<-data.frame(position=kmer_starts,value_min=rep(min_plot_value,length(kmer_starts)),value_max=rep(max_plot_value,length(kmer_starts)))
    #plot distribution of given kmer using geom_linerange from ggplot2
    plot <- ggplot2::ggplot(positions_df,ggplot2::aes(x=position)) + ggplot2::geom_linerange(ggplot2::aes(ymin=value_min,ymax=value_max)) + ggplot2::scale_y_continuous(limits = c(-2,2)) + ggplot2::ggtitle(paste(kmer," occurences in sequence (",length(kmer_starts)," times)"))
    return(plot)
  }
  else {
    #if pattern was not found throw error
    stop("Can't find provided pattern in the sequence")
  }
}
