context("Test kmer calculation functions")
library(biostrext)
library(Biostrings)
library(assertthat)


test_that("kmers are properly found in DNAString objects",{
  expect_error(get_max_kmer_single(test_DNAStringSet,3))
  expect_is(get_max_kmer_single(test_DNAString1,3),"integer")
  expect_equal(unname(get_max_kmer_single(test_DNAString1,3)),22)
  expect_equal(names(get_max_kmer_single(test_DNAString1,3)),"CGA")
  expect_error(get_max_kmer_single(test_DNAString1))
  expect_error(get_max_kmer_single(test_DNAString1,11))
  expect_error(get_max_kmer_single(test_noDNAString1,3))
  expect_equal(unname(get_max_kmer(test_DNAString1,3)),22)
  expect_equal(names(get_max_kmer(test_DNAString1,3)),"CGA")
  expect_error(get_max_kmer(test_DNAString1))
  expect_error(get_max_kmer(test_DNAString1,11))
  expect_error(get_max_kmer(test_noDNAString1,3))
  expect_error(get_max_kmer_single(test_DNAString3,10))
  expect_error(get_max_kmer(test_DNAString3,10))
  expect_error(get_max_kmer_single(test_DNAString2,16,force_large_kmer = TRUE))
  expect_error(get_max_kmer_single(test_DNAString1,0))
  expect_error(get_max_kmer_single(test_DNAString1,-1))
})


test_that("kmers are properly found in DNAStringSet objects",{
  expect_is(get_max_kmer(test_DNAStringSet,3),"numeric")
  expect_equal(unname(get_max_kmer(test_DNAStringSet,3)),c(22,10,1))
  expect_equal(names(get_max_kmer(test_DNAStringSet,3)),c("CGA","AGC","ACA"))
  expect_error(get_max_kmer(test_DNAStringSet))
  expect_error(get_max_kmer(test_DNAStringSet,11))
  expect_error(get_max_kmer(test_DNAString1_seq,3))
  expect_error(get_max_kmer(test_DNAStringSet,0))
  expect_error(get_max_kmer(test_DNAStringSet,-2))
})


test_that("kmer distribution is properly plotted",{
  expect_is(plot_kmer_distribution(test_DNAString1,kmer_seq),"gg")
  expect_error(plot_kmer_distribution(test_DNAString1,kmer_seq_not_present))
  expect_error(plot_kmer_distribution(test_DNAString1_seq,kmer_seq))
  expect_error(plot_kmer_distribution(test_DNAStringSet,kmer_seq))
  expect_error(plot_kmer_distribution(test_DNAStringSet,""))
})

