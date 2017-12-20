
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

#' #' Produce error if is_DNAstring is false
#' #'
#' #' @param call is_DNAstring invocation
#' #' @param env
#' #'
#' #' @return Error if is_DNAstring call is false
#' #'
#' on_failure(is_DNAstring) <- function(call, env) {
#'   paste0(deparse(call$x), " is not a DNAString object")
#' }

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

#' #' Produce error if is_DNAstring is false
#' #'
#' #' @param call is_DNAstring invocation
#' #' @param env
#' #'
#' #' @return Error if is_DNAstring call is false
#' #'
#' on_failure(is_DNAstringSet) <- function(call, env) {
#'   paste0(deparse(call$x), " is not a DNAStringSet object")
#' }


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
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  assertthat::assert_that(is_DNAstring(x))

  seq_length = length(x)
  if(seq_length>=4) {
    GC <- sum(Biostrings::alphabetFrequency(x)[c(2,3)])/seq_length*100
  }
  else{
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

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  assertthat::assert_that(is_DNAstring(x))
  assertthat::assert_that(window<=length(x))
  assertthat::assert_that(window>0)
  assertthat::assert_that(overlap>=0)
  assertthat::assert_that(overlap<window)

  seq_length = length(x)
  return_vector <- character(0)
  start=1
  end=window
  z = 1
  while (end<=seq_length) {
    temp_seq<-as.character(Biostrings::subseq(x,start,end))
    return_vector[z]<-temp_seq
    z=z+1
    start=start+window-overlap
    end=end+window-overlap

  }
  if(start<seq_length) {
    temp_seq<-as.character((Biostrings::subseq(x,start,seq_length)))
    return_vector[z]<-temp_seq
  }
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
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  assertthat::assert_that(is_DNAstring(x) || is_DNAstringSet(x))
  if(is_DNAstring(x)) {
    value = calc_GC_single((x))
  }
  else if (is_DNAstringSet(x)) {
    set_length = length(x)
    value <- numeric(0)
    for (z in seq(1,set_length)) {
      value_temp <- calc_GC_single(x[[z]])
      value<-c(value,value_temp)
    }
  }
  return(value)
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
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  assertthat::assert_that(is_DNAstring(x))
  assertthat::assert_that(window<=length(x))
  splitted_seq <- split_in_windows(x,window,overlap)
  GC <- calc_GC(splitted_seq)
  return(GC)
}

#' Title Get kmer which is most frequent in provided sequence
#'
#' @param x DNAString object
#' @param mer kmer length
#' @param force_large_kmer logical, should calcualtion of large kmers be possible (but shorter than 15)
#'
#' @return named numeric, with count of given kmer and its sequence as name
#' @export
#'
#' @examples
#' get_max_kmer_single(Biostrings::DNAString("CAGCTGCACACATCGTACA"),3)
get_max_kmer_single <- function(x,mer,force_large_kmer = FALSE) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  assertthat::assert_that(is_DNAstring(x))
  assertthat::assert_that(mer>0)
  assertthat::assert_that(mer<length(x))
  if (mer>10 && !isTRUE(force_large_kmer)) {
    stop("Calculation of kmers larger than 10 is computationally intensive and may take a long time.
If you need that, specify force_large_kmer=TRUE (use with caution!)",
         call. = FALSE)
  }
  if (mer>15 && isTRUE(force_large_kmer)) {
    stop("Calculation of kmers larger than 15 is extremally computationally intensive and it is discouraged to use.",
         call. = FALSE)
  }
  oligoFreq <- Biostrings::oligonucleotideFrequency(x,mer)
  max_kmer_count <- max(oligoFreq)
  max_kmer_seq <- names(which.max(oligoFreq))
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
get_max_kmer <- function(x,mer) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  assertthat::assert_that(is_DNAstring(x) || is_DNAstringSet(x))
  assertthat::assert_that(mer>0)
  if(is_DNAstring(x)) {
    value = get_max_kmer_single(x,mer)
  }
  else if (is_DNAstringSet(x)) {
    set_length = length(x)
    value <- numeric(0)
    for (z in seq(1,set_length)) {
      if (length(x[[z]])>mer) {
        value_temp <- get_max_kmer_single(x[[z]],mer)
      }
      else {
        value_temp <- 'NA'
      }
      value<-c(value,value_temp)
    }
  }
  return(value)
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
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  assertthat::assert_that(is_DNAstring(x))
  assertthat::assert_that(length(x)/(window-overlap)>20)
  splitted_seq <- split_in_windows(x,window,overlap)
  step=window-overlap
  positions<-seq(1,length(x)-window+step+1,step)
  GC <- calc_GC(splitted_seq)
  names(GC)<-positions
  GC_melt <- reshape::melt(GC)
  GC_melt[['position']]<-as.numeric(row.names(GC_melt))
  plot <- ggplot2::ggplot(GC_melt,ggplot2::aes(x=position,y=value)) + ggplot2::geom_smooth(span=0.1) + ggplot2::scale_y_continuous() + ggplot2::ggtitle("%GC content distribution accross sequence")
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
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  assertthat::assert_that(is_DNAstring(x))
  if(kmer=='') {
    stop("Empty kmer sequence provided",
         call. = FALSE)
  }
  kmer_matches<-Biostrings::matchPattern(kmer,x)
  if (length(kmer_matches)>0) {
    kmer_starts=Biostrings::start(kmer_matches)
    positions_df<-data.frame(position=kmer_starts,value_min=rep(-0.1,length(kmer_starts)),value_max=rep(0.1,length(kmer_starts)))
    plot <- ggplot2::ggplot(positions_df,ggplot2::aes(x=position)) + ggplot2::geom_linerange(ggplot2::aes(ymin=value_min,ymax=value_max)) + ggplot2::scale_y_continuous(limits = c(-2,2)) + ggplot2::ggtitle(paste(kmer," occurences in sequence (",length(kmer_starts)," times)"))
    return(plot)
  }
  else {
    stop("Can't find provided pattern in the sequence")
  }
}
