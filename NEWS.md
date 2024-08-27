# gsignal 0.3.6

* date: 20240827
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

---

# gsignal 0.3-5

* date: 20220514
* Fixed  Discussion #6: remove padding to nearest power of 2 in `pwelch()`
* Fixed Github Issue #5: returning matrix when input is matrix in `pwelch()`
* Adapted e-mail addresses in mexihat, morlet, nutallwin, pburg, pyulear
* use `inherits()` instead of direct comparison of class name in ar_psd,
    findpeaks, pwelch, sgolayfilt, upfirdn
* defined plot methods for ar_psd, pwelch, specgram classes
* Fixed Github Issue #7: decimate with a matrix;, added "fir" argument to ftype

---

# gsignal 0.3-4

* date: 20220404
* Fixed test failure in tests/testthat/test_miscellaneous_Functions.R

---

# gsignal 0.3-3

* date: 20220330
* Fixed Github Issue #3: Problems with `fftfilt()` when FFT length is 
    provided by user
* copy attributes of input object x to output in functions filter, filtfilt,
    sosfilt, fftfilt
* copy dimnames of input object x to output in functions upfirdn, resample,
    upsample, upsamplefill, downsample, decimate, detrend, fht, sgolayfilt 
* added `ultrwin()` function
* adapted `filter()` to allow data and filter coefficients to be of type
    complex
* adapted `sosfilt()` to allow data and filter coefficients to be of type
    complex
* bugfix in `pwelch()` for multivariate input
* Fixed Github Issue #4: Problem with `hilbert()` for small amplitude signals
* added `isConjSymm()` function to gsignal-internal
* adapted `ifft()` to use `isConjSymm()` instead of `ZapIm()`
* reduced default tolerance for `isWhole()` and `zapIm()`
* bugfix in `detrend()`: function now returns a vector if input was a vector
* bugfix in `filtfilt()`: corrected bug in computing filter ends
    (default and Sos methods)

---

# gsignal 0.3-2

* date: 20210518
* corrected CRAN WARNINGs on ATLAS, MKL, valgrind, fedora, solaris
* corrected import NOTE for grDevices
* adapted code in vignette "gsignal"
* use explicit tolerance in testthat tests
* added badges and logo (just for the fun of it)
* sort zeros and poles on output in sos2zp(), tf2zp()
* use `matrix()` instead of `as.matrix()` in functions
    `dct()`, `idct()`, `czt()`, `dst()`, `idst()`, `fht()`, `ifht()`
* minor bugfix in `arburg()`
* bugfix in filter.cpp: resize `a` and `b` vectors, length of `zi`
* adapted some examples

---

# gsignal 0.3-1

- date: 20210502
- added 'signals' data frame
- cleaned code
- minor bugfixes
- updated tests
- added vignette

---

# gsignal 0.3-0

- date: 20210411
- Bugfixes in pwelch(), butter(), as.Arma.Sos(), as.Sos.Zpg(), sos2tf()
- Added 'output' parameter to butter(), cheby1(), cheby2(), ellip()
- Redesigned filter() to correct problems with initial conditions (now direct-form II)
- Return objects of respective class tf2zp(), tf2sos(), zp2sos(), zp2tf(), sos2tf(), sos2zp()
- Added 'order' argument to zp2sos(); changed default ordering
- Added filter_zi() function
- Redesigned sosfilt(); added support for initial conditions
- Adapted filtfilt(): use padding and Gustafsson method for initial conditions
- Function sgolayfilt(): if input x is a matrix, filter its columns
- Replace 'dim' argument by 'MARGIN' in functions cplxpair(), cplxpair(), medfilt1()
- Solved numerical precision error in chebwin() and cheb() functions
- Bugfix in cplxpair(); added additional tests
- Return vector if input is a vector in function czt()
- Added range parameter ('half' or 'whole') in function pwelch()
- Bugfix in calculating time points, function stft()
- Changed plotting to S3 functions in function specgram()

---

# gsignal 0.2-0

- date: 20201218
- completed documentation
- adapted examples and links
- license changed to GPL-3 in accordance with Octave signal package license

---

# gsignal 0.1-0

- date: 20201213
- initial setup. Functions:

- **Signals**:
  buffer, chirp, cmorwavf, diric, gauspuls, gmonopuls, mexihat, meyeraux, morlet, pulstran, rectpuls,
  sawtooth, shanwavf, shiftdata, sigmoid_train, sinetone, sinewave, specgram, square, tripuls,
  udecode, uencode, unshiftdata
- **Signal Measurement**:
  findpeaks, peak2peak, peak2rms, rms, rssq
- **Correlation and Convolution**:
  cconv, convmtx, wconv, xcorr, xcorr2, xcov
- **Filtering**:
  fftfilt, filter, filter2, filtfilt, filtic, medfilt1, movingrms, sgolayfilt, sosfilt
- **Filter Analysis**:
  freqs, freqs_plot, freqz, freqz_plot, fwhm, grpdelay, impz, zplane
- **Filter Conversion**:
  polystab, residued, residuez, sos2tf, sos2zp, tf2sos, tf2zp, zp2sos, zp2tf
- **IIR Filter Design**:
  besselap, besself, bilinear, buttap, butter, buttord, cheb, cheb1ap, cheb1ord, cheb2ap, cheb2ord, cheby1, cheby2
  ellip, ellipap, ellipord, iirlp2mb, impinvar, invimpinvar, ncauer, pei_tseng_notch, sftrans
- **FIR Filter Design**:
  cl2bp, fir1, fir2, firls, kaiserord, qp_kaiser, remez, sgolay
- **Transforms**:
  bitrevorder, cceps, cplxreal, czt, dct, dct2, dctmtx, dftmtx, digitrevorder, dst, dwt, fftshift, fht, fwht,
  hilbert, idct, idct2, idst, ifft, ifftshift, ifht, ifwht, imvfft, rceps, stft
- **Power Spectrum Analysis**:
  ar_psd, cohere, cpsd, csd, db2pow, mscohere, pburg, pow2db, pwelch, pyulear, tfe, tfestimate
- **Window Functions**:
  barthannwin, bartlett, blackman, blackmanharris, blackmannuttall, bohmanwin, boxcar, chebwin, flattopwin,
  gaussian, gausswin, hamming, hanning, hann, kaiser, nuttallwin, parzenwin, rectwin, triang, tukeywin, ultrwin, welchwin
- **System Identification**:
  arburg, aryule, invfreq, invfreqs, invfreqz, levinson
- **Sample Rate Change**:
  decimate, downsample, interp, resample, upfirdn, upsample
- **Utility**:
  clustersegment, fracshift, marcumq, primitive, sampled2continuous, schtrig, upsamplefill, wkeep, zerocrossing
- **Standard Functions**:
  detrend, pad, postpad, prepad
- **Miscellaneous**:
  Arma, Ma, Zpg, Sos, FilterSpecs,
  fftconv, unwrap, cplxpair, poly, conv, conv2, residue, polyreduce, mpoles
