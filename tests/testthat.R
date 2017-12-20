library(testthat)
library(biostrext)

test_DNAString1_seq <- "AGTCATGCATCGATCGATGCATCGATCGATCGATCGATCGATCGATGCTAGCTAGCTACGATCGATCACGATGCATGCTAGCTAGCTAGCTAGCTCGATCGATGCATCGATCGATCGATCGTAGCTACGTAGCCATGCTACGTACGTACGTACGATCGATCGATCGATGCTAGCTACGATCGTACGTAGCTACGATCGTA"
test_DNAString1 <- Biostrings::DNAString(test_DNAString1_seq)
test_DNAString2_seq <- "CAGTCAGCATCGATCGTACGTAGCTAGCTGCACCAACACCATCGATCGATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGTATGCATACCGGGTAAT"
test_DNAString2 <- Biostrings::DNAString(test_DNAString2_seq)
test_DNAString3_seq <- "CAGCTACAC"
test_DNAString3 <- Biostrings::DNAString(test_DNAString3_seq)
test_DNAString_short_seq <- "AG"
test_DNAString_short <- Biostrings::DNAString(test_DNAString_short_seq)
test_DNAString_GC_seq <- "GCGCGCGCGCGCGCGCGCGCGCCCCGCGGC"
test_DNAString_GC <- Biostrings::DNAString(test_DNAString_GC_seq)
test_DNAString_AT_seq <- "ATATATATATATATATTATATATATATAA"
test_DNAString_AT <- Biostrings::DNAString(test_DNAString_AT_seq)
test_DNAStringSet <- Biostrings::DNAStringSet(c(test_DNAString1_seq,test_DNAString2_seq,test_DNAString3_seq))
kmer_seq = "GCAT"
kmer_seq_not_present = "CAGTC"
kmer_seq_empty = ""
kmer_seq_wrong_letters = "ARG"

test_check("biostrext")
