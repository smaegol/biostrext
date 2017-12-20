# biostrext
Simple extension to [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) R package

[![Build Status](https://travis-ci.org/smaegol/biostrext.svg?branch=master)](https://travis-ci.org/smaegol/biostrext)
[![codecov](https://codecov.io/gh/smaegol/biostrext/branch/master/graph/badge.svg)](https://codecov.io/gh/smaegol/biostrext)

## Description

Provide simple functions easing %GC calculation and plotting. %GC can be calculated for single DNAString or DNAStringSet objects. For the DNAString objects it is possible to calculate %GC using sliding window approach and plot %GC distribution with ggplot.

Another functionality allows to find the most frequent pattern of given length (k) in the provided sequence. It also allows plotting the positions of given pattern in the sequence.

## Functions

#### split_in_windows(x,window,overlap)

Splits a DNAString object using sliding window approach. Window parameter defines the size of the window whereas overlap defines the overlap between consecuting windows.


#### calc_GC_single(x)

Returns %GC for provided DNAString object (x)

#### calc_GC(x)

Returns %GC for provided DNAString or DNAStringSet object (x). In the case of DNAStringSet the vector containing %GC for each sequence in the set is returned

#### calc_windowed_GC(x,window,overlap)

Calculates %GC of provided DNAString object (x) using sliding window approach. Window parameter defines the size of the window whereas overlap defines the overlap between consecuting windows.
Return numeric vector containing %GC for each created chunk

#### plot_GC_distribution(x,window,overlap)

Calculates %GC of provided DNAString object (x) using sliding window approach and creates the plot. Window parameter defines the size of the window whereas overlap defines the overlap between consecuting windows.

#### get_max_kmer(x,k)

Finds the most frequently occuring pattern of size k in the DNAString object (x)
