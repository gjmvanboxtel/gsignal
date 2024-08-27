## Resubmission
This is a resubmission. In this version I have:

* Fixed issue #9: differences between `freqz()` and `signal.scipy.freqz`
  - updated `freqz()` to match current Octave version
  - implemented proper `freqz.Sos()` instead of converting to Arma
  - `freqz.Zpg()` converts to Sos instead of Arma
* Implemented changes in 'Octave' 'signal' 1.4.5 with respect to 1.4.1
  - added tests for `cplxreal()` (gsignal did not have Octave bug #60606)
  - added test for `cheb2ap()` (gsignal did not have Octave bug #58218)
  - all other changes concerned code style or documentation
* correct typos in stft.R (H. Dieter Wilhelm - pull request #10)
* Negate the signal reversed in time at both ends in `filtfilt()`
  (Rafael Laboissière - pull request #12)
* add check on NULL value in `zapIm()`
* Fixed issue #15 (loeriver): incorrect results and error in `residuez()`
  - removed `Cong()` on calculation of k
  - made calculation of `r` and `p` conditional
* Changed coefficients of the `hamming()` window function as in Matlab/Octave
* Fixed issue #17 (dipterix): `decimate()` not compatible with Matlab
  - replaced call to `fftfilt()` with `filtfilt()` if FIR filter is requested
* Bugfix in `fir1()`: adapted calculation of w_o
* Fixed issue #19 (jefferis): changed control flow logic in `findpeaks()`
  (line 174)
* Adapt `freqs()` to match `freqz()` - delete freqs_plot.R and add function
  `freqs_plot()` to freqs.R.


## R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 note ✖
Version contains large components (0.3-6.9010)

 * This is a big package with nearly 200 functions and C++ code
   with large object files. It also include example data. These components
   will rarely be updated.

## revdepcheck results

We checked 9 reverse dependencies (5 from CRAN + 4 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages




