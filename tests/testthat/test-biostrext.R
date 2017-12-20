context("Test main package functions")
library(biostrext)
library(Biostrings)
library(assertthat)


test_that("DNAString is properly recognized",{
  expect_true(is_DNAstring(test_DNAString1))
  expect_false(is_DNAstring(test_DNAString1_seq))
})


test_that("DNAString is properly splitted in chunks",{
  expect_is(split_in_windows(test_DNAString1,15,1),"DNAStringSet")
  expect_error(split_in_windows(test_DNAString1,10,-3))
  expect_error(split_in_windows(test_DNAString1,-2,-3))
})
