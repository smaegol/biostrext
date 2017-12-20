context("test-biostrext.R")
library(biostrext)
library(Biostrings)
library(assertthat)

test_DNAString <- Biostrings::DNAString("AGTCATGCATCGATCGATGCATCGATCGATCGATCGATCGATCGATGCTAGCTAGCTACGATCGATCACGATGCATGCTAGCTAGCTAGCTAGCTCGATCGATGCATCGATCGATCGATCGTAGCTACGTAGCCATGCTACGTACGTACGTACGATCGATCGATCGATGCTAGCTACGATCGTACGTAGCTACGATCGTA")
test_DNAString_short <- Biostrings::DNAString("AG")
test_noDNAString <- "AGTCATGCATCGATCGATGCATCGATCGATCGATCGATCGATCGATGCTAGCTAGCTACGATCGATCACGATGCATGCTAGCTAGCTAGCTAGCTCGATCGATGCATCGATCGATCGATCGTAGCTACGTAGCCATGCTACGTACGTACGTACGATCGATCGATCGATGCTAGCTACGATCGTACGTAGCTACGATCGTA"

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("DNAString is properly recognized",{
  expect_true(is_DNAstring(test_DNAString))
  expect_false(is_DNAstring(test_noDNAString))
})

test_that("GC is properly calculated",{
  expect_is(calc_GC_single(test_DNAString),"numeric")
  expect_equal(length(test_DNAString),200)
  expect_equal(calc_GC_single(test_DNAString),50)
  expect_equal(length(test_DNAString_short),2)
  expect_equal(calc_GC_single(test_DNAString_short),NA)
})


