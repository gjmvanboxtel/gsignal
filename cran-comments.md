## Resubmission of gsignal v0.3-5

> Thanks, we see:

> Found the following (possibly) invalid URLs:
    URL: https://codecov.io/gh/gjmvanboxtel/gsignal (moved to https://app.codecov.io/gh/gjmvanboxtel/gsignal)
      From: README.md
      Status: 301
      Message: Moved Permanently

> Please change http --> https, add trailing slashes, or follow moved content as appropriate.

> Please fix and resubmit. 

Sorry for missing that one.
Now fixed.
Thank you,
Geert van Boxtel

## Submission of gsignal v0.3-5

Changes with respect to v0.3-4:
- Fixed Github Discussion #6: remove padding to nearest power of 2 in pwelch()
- Fixed Github Issue #5: returning matrix when input is matrix in pwelch()
- Adapted e-mail addresses in mexihat, morlet, nutallwin, pburg, pyulear
- use inherits() instead of direct comparison of class name in ar_psd, findpeaks,
    pwelch, sgolayfilt, upfirdn
- defined plot methods for ar_psd, pwelch, specgram classes
- Fixed Github Issue #7: decimate with a matrix;, added "fir" argument to ftype

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

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
 