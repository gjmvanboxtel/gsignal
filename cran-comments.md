## Submission of gsignal v0.3-3

Changes with respect to v0.3-2 (on CRAN):
- Fixed Github Issue #3: Problems with fftfilt when FFT length is provided by user
- copy attributes of input object x to output in functions filter, filtfilt, sosfilt, fftfilt
- copy dimnames of input object x to output in functions upfirdn, resample, upsample, upsamplefill,
       downsample, decimate, detrend, fht, sgolayfilt, 
- added ultrwin() function
- adapted filter() to allow data and filter coefficients to be of type complex
- adapted sosfilt() to allow data and filter coefficients to be of type complex
- bugfix in pwelch() for multivariate input
- Fixed Github Issue #4: Problem with hilbert for small amplitude signals
- added isConjSymm() function to gsignal-internal
- adapted ifft() to use isConjSymm instead of ZapIm
- reduced default tolerance for isWhole() and zapIm()
- bugfix in detrend(): function now returns a vector if input was a vector
- bugfix in filtfilt(): corrected bug in computing filter ends (default and Sos methods)

## Test environments
- Windows Server 2022, R-devel, 64 bit
- Fedora Linux, R-devel, clang, gfortran
- Ubuntu Linux 20.04.1 LTS, R-release, GCC
- Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results

All tests return 'Status: success'

There were three NOTES:

>* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Geert van Boxtel <G.J.M.vanBoxtel@gmail.com>’

This note always appears when I submit a package for testing

>Version contains large components (0.3-3.9000)
* checking installed package size ... NOTE
  installed size is  5.6Mb
  sub-directories of 1Mb or more:
    doc    1.0Mb
    libs   3.3Mb

It is a large package with 196 objects, compiled C++ code and a vignette.

>* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

In https://github.com/r-hub/rhub/issues/503, gaborcsardi noted:
Seems like a bug/crash in miktex, so you can ignore this.

## revdepcheck results

revdepcheck::revdep_check() reports: *Wow, no problems at all. :)*

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
 