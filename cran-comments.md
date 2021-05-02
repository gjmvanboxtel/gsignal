## Resubmission
This is a resubmission. The following changes were made:
* More details about the package functionality and implemented methods
  were added in Description text
* Package names, software names and API names were written in single quotes
* It was checked that all .Rd files contain a \value{} paragraph, even if
  a function does not return a value
* commented source code in examples was removed or uncommented
* It was ensured that the global environment was not modified.

I then ran devtools::check() and rhub::check_for_cran() again.
Thank you.

## Test environments
* local R installation, R 4.0.5
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.

This is my first submission to CRAN

* checking installed package size ... NOTE
  installed size is  5.4Mb
  sub-directories of 1Mb or more:
    doc    1.0Mb
    libs   3.2Mb

It is a big package consisting of 334 objects, including compiled C++ code and
an extensive vignette.

## Downstream dependencies
There are currently no downstream dependencies for this package

