context("Test %GC related functions")
library(biostrext)
library(Biostrings)
library(assertthat)


test_that("GC is properly calculated for DNAString",{
  expect_is(calc_GC_single(test_DNAString1),"numeric")
  expect_equal(length(test_DNAString1),200)
  expect_equal(calc_GC_single(test_DNAString1),50)
  expect_equal(calc_GC_single(test_DNAString2),50.49505,tolerance=1e-4)
  expect_equal(calc_GC_single(test_DNAString3),55.55556,tolerance=1e-4)
  expect_equal(length(test_DNAString_short),2)
  expect_equal(calc_GC_single(test_DNAString_short),NA)
  expect_error(calc_GC_single(test_DNAString1_seq))
  expect_equal(calc_GC_single(test_DNAString_AT),0)
  expect_equal(calc_GC_single(test_DNAString_GC),100)
  expect_is(calc_GC(test_DNAString1),"numeric")
  expect_equal(calc_GC(test_DNAString1),50)
})

test_that("GC is properly calculated for chunks",{
  expect_is(calc_GC(test_DNAStringSet),"numeric")
  expect_equal(calc_GC(test_DNAStringSet),c(50,50.49505,55.55556),tolerance=1e-4)
  expect_equal(calc_windowed_GC(test_DNAString1,22,2),c(45.45455,45.45455,45.45455,50,54.54545,50,54.54545,54.54545,54.54545,50),tolerance=1e-4)
  expect_error(calc_windowed_GC(test_DNAStringSet,22,2))
  expect_error(calc_windowed_GC(test_DNAString1,22,22))
  expect_error(calc_windowed_GC(test_DNAString3,22,2))
  expect_error(calc_windowed_GC(test_DNAString1,22,-1))
  expect_error(calc_windowed_GC(test_DNAString1,0,2))
  expect_error(calc_windowed_GC(test_DNAString1,-2,-3))
})

test_that("GC is properly plotted",{
  expect_is(plot_GC_distribution(test_DNAString1,5,2),"gg")
  expect_error(plot_GC_distribution(test_DNAStringSet,5,2))
  expect_error(plot_GC_distribution(test_DNAString1,12,0))
  expect_error(plot_GC_distribution(test_DNAString1,10,-1))
  expect_error(plot_GC_distribution(test_DNAString1,-2,-3))
})
