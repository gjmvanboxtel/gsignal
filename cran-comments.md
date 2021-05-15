## Resubmission

This is a resubmission. In this version an issue resulting from a valgrind check was fixed:
==1100340== Invalid read of size 8
==1100340==    at 0x48C55C9: rfilter(Rcpp::Vector<14, Rcpp::PreserveStorage>, Rcpp::Vector<14, Rcpp::PreserveStorage>, Rcpp::Vector<14, Rcpp::PreserveStorage>, Rcpp::Vector<14, Rcpp::PreserveStorage>) (/tmp/gsignal.Rcheck/00_pkg_src/gsignal/src/filter.cpp:40) 

Thank you.

## Submission of gsignal v0.3-2

* corrected CRAN WARNINGs on ATLAS, MKL, valgrind, fedora, solaris:
- corrected import NOTE for grDevices
- adapted code in vignette "gsignal"
- use explicit tolerance in testthat tests
- added badges and logo (just for the fun of it)
- sort zeros and poles on output in sos2zp(), tf2zp()
- use matrix() instead of as.matrix() in functions
    dct(), idct(), czt(), dst(), idst(), fht(), ifht()
- minor bugfix in arburg()

## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)
- R-hub linux-x86_64-rocker-gcc-san (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  
  Maintainer: 'Geert van Boxtel <G.J.M.vanBoxtel@gmail.com>'
  Version contains large components (0.3-2.9000)

> On ubuntu-gcc-release (r-release)
  checking installed package size ... NOTE
    installed size is  5.5Mb
    sub-directories of 1Mb or more:
      doc    1.0Mb
      libs   3.3Mb

It is a large package with 195 objects, compiled C++ code and a vignette.

0 errors ✓ | 0 warnings ✓ | 2 notes x

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
