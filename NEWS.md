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
- licence changed to GPL3 in accordance with Octave signal package licence

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
