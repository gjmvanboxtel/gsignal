//     mmfir.c: MiniMax FIR filter design per McClellan, Parks & Rabiner
//      Copyright (c) 2014/15/21 Rob Sykes <robs@users.sourceforge.net>
//
// mmfir is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// mmfir is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

/*
 * GvB 20261111:
 * Adapted diagnostic output for use in an R package.
 * Replaced direct stderr output with R's console interface.
 */

#include "mmfir.h"
#include <R_ext/Print.h>


#if defined QUAD_PREC
  #include <quadmath.h>
#endif
#if defined _OPENMP
  #include <omp.h>
#endif

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if INT_MAX == 2147483647
  typedef int rand_t;
#elif LONG_MAX == 2147483647L
  typedef long rand_t;
#else
  #error no suitable int32 type
#endif

#if defined _MSC_VER
  #define inline __inline
#endif



// Real-number precision, quad, double or long double (the default).

#if defined QUAD_PREC // Very slow without native support.
  typedef __float128 real;
  #define cosr cosq
  #define sinr sinq
  #define PIr M_PIq
#elif defined DOUBLE_PREC
  typedef double real;
  #define cosr cos
  #define sinr sin
  #define PIr PI
#else
  typedef long double real;
  #define cosr cosl
  #define sinr sinl
  #define PIr 3.1415926535897932384626433832795029L
#endif



#define bandL           (&bandSpecs[0])
#define bandR           (&bandSpecs[numBands-1])
#define clamp(x,lo,hi)  (((x)>(hi))?(hi):(((x)<(lo))?(lo):(x)))
#define debugLevel      (flags & 3)
#define dir(x)          ((x)<0? -1: (x)>0? 1 : 0)
#define epsilon         1e-11 // In some circumstances, treat < epsilon as 0
#define E(x)            peaks[x].e
#define e(x)            space[x].e
#define F(x)            peaks[x].f
#define f(x)            space[x].f
#define maxExtras       (numBands+2)
#define maxPeaks        (R+1+maxExtras)
#define PI              3.1415926535897932384626433832795029
#define rand(r,x)       (ranqd1(r) * ((x) / (65536. * 32768.))) // in [-x,x)
#define ranqd1(r)       ((r) = 1664525 * (r) + 1013904223) // int32_t r
#define returnError(x)  do {report.result=x; goto END;} while(0)
#define iround(x)       (int)floor((x)+.5)
#define setF(F)         for (f=(F); f>b->freqR; ++b)
#define userEpsilon     1e-6
#define hasConverged(len) (((minimality = minimality_(peaks, len, \
 minimalityThreshold))>0 && (density<density2 || stability>2)) || stability>2.7)



static inline double maxD(double a, double b) {return a>=b? a:b;}
static inline int    maxI(int    a, int    b) {return a>=b? a:b;}

#if defined _OPENMP
static int ompThreadCount(size_t iterations, size_t work)
{
  size_t const grain = 8192;
  int const max = omp_get_max_threads();
  size_t n = (work + grain - 1) / grain;

  if (n > iterations)
    n = iterations;
  if (max < 2 || n < 2)
    return 1;
  return n < (size_t)max? (int)n : max;
}
#endif



static char const * const errorText[]={
  "out-of-memory",
  "number of bands must be at least 1",
  "filter-order must be at least 1",
  "band-edge frequencies must increase in [0,1]",
  "weighting function out-of-range",
  "amplitude function out-of-range",
  "type III/IV DC amplitude response must be 0",
  "type II/III Nyquist amplitude response must be 0",
};



// Type definitions:

typedef enum {TypeI, TypeII, TypeIII, TypeIV} FirType;

typedef struct {real x, beta, gamma;} LagrangeCoef;

typedef MmfirPoint Point;



// Determine coefs of a 2nd-order polynomial from 3 of its points:

static void poly2(
    double x0, double x1, double x2,
    double y0, double y1, double y2,
    double * a, double * b, double * c)
{
  double a0 = y0/((x0-x1)*(x0-x2));
  double a1 = y1/((x1-x2)*(x1-x0));
  double a2 = y2/((x2-x0)*(x2-x1));
  *a = a0+a1+a2;
  *b = -(a0*(x1+x2) + a1*(x0+x2) + a2*(x0+x1));
  *c = a0*x1*x2 + x0*a1*x2 + x0*x1*a2;
}



// Evaluate a 2nd-order polynomial at given point:

static double poly2Val(double a, double b, double c, double x)
{
  return fma(fma(a, x, b), x, c);
}



/* For filter-response analysis, use a discrete frequency space that is
 * 'warped' around transition-bands.  The warping function is derived from the
 * observation by Shpak/Antoniou (1990) that an 'initial guess' for the peak
 * frequencies based on adjacent-peak spacing proportional to s(n) =
 * exp(c1+c2*n), improves upon one where the adjacent-peak spacing is constant.
 *
 * E.g., for a band-pass filter with band-edges 0, s1, p1, p2, s2, 1, the
 * analysis points within warped space compared to those in a linear space
 * might look like the following (second stop-band not shown):
 *          0                       s1    p1                                  p2
 * Warped:  o      o     o    o  o oo     oo o  o    o      o      o    o  o oo
 * Linear:  o o o o o o o o o o o o o     o o o o o o o o o o o o o o o o o o o
 *
 * Warp function derivation is as follows: approx. peak position continuous
 * function, warp(x) =
 *   x>=1: f(x)-f(0)
 *   x<1:  2nd-order approximation (since exp. approx. breaks down here)
 * where f(x) = Indef. integral of s(x) = exp(c1+c2*x)/c2
 *
 * Constants c1, c2 are determined empirically.  Dependent constants c3=-f(0),
 * ca,cb,cc are calculated as follows: c3 = -exp(c1)/c2;
 *   poly2(0,1,2, 0,warp(1),warp(2), &ca,&cb,&cc);
 */

static double const c1 = +.1235;
static double const c2 = -.8;

static double const c3 = 1.41431;
static double const ca = .214437;
static double const cb = .006742;
static double const cc = 0;



static double warp(double x) // Calculate warp in arbitrary units:
{
  double t = c1 + c2*x;
  return x<0? 0 : x<1? poly2Val(ca,cb,cc,x) : x-(t>-38? exp(t)/c2 : 0)-c3;
}



static double iwarp(double y) // Inverse of warp().
{
  double t, x1=y+c3, x;

  if (y<0)
    return 0;

  if (y<warp(1))
    return (-cb + sqrt(cb*cb-4*ca*(cc-y)))/(2*ca);

  // Newton iteration:
  do x=x1, t=c1+c2*x, x1=x-(warp(x)-y)/(1-(t>-38? exp(t) : 0));
  while (fabs(x1-x)>1e-5);
  return x1;
}



// Generate a set of warped frequency-space points.
// N.B. rand helps to avoid Lagrange interpolation instability:

static int warpSpace(Point * p, double m, double f1, double f2,
    int density, int step, rand_t * random)
{
  int i, n = (int)ceil(m*density);
  double scale, df = f2 - f1;

  p->f = f1, p+=step; // End points must be exact (not calculated).
  if (n > 0) {
    scale = df/warp(1.*n/density);
    for (i = 1; i < n; ++i)
      p->f = f1+scale*warp((i+rand(*random, .001))/density), p+=step;
    p->f = f2; // End points must be exact (not calculated).
  }
  return 1+n;
}



static int warpSpaceR( // Right of transition-band.
    Point * p, double m, double f1, double f2, int density, rand_t * random)
{
  return warpSpace(p, m, f1, f2, density, +1, random);
}



static int warpSpaceL ( // Left of transition-band.
    Point * p, double m, double f1, double f2, int density, rand_t * random)
{
  return warpSpace(p+(int)ceil(m*density), m, f2, f1, density, -1, random);
}



static int warpSpaceRL( // Between transition-bands.
    Point * p, double m, double f1, double f2, int density, rand_t * random)
{
  double df1 = (f2 - f1) * .5, m1 = m * .5;
  int n = warpSpaceR(p, m1, f1, f1+df1, density, random);
  warpSpaceL(p + n-1, m1, f1+df1, f2, density, random);
  return n*2-1;
}



// Apportion space to each band, taking account of warping:

static void apportionSpace(MmfirBandSpec * b, int numBands, int n)
{
  int i, changed;
  double sum, t, av;

  do { // Use an iterative method.
    for (av=changed=0, i=0; i<numBands; ++i)
      av += warp(b[i].portion*n/b[i].ends)/(b[i].portion0/b[i].ends);
    av/=numBands;

    for (sum=i=0; i<numBands; ++i) {
      t=iwarp(av*b[i].portion0/b[i].ends)/n*b[i].ends;
      changed += fabs(t-b[i].portion) > 1e-6;
      sum+= b[i].portion=t;
    }
    for (i=0; i<numBands; b[i++].portion/=sum); // Normalise.
  } while (changed);
}



static inline real // Modifier at given f in [0,1], for filter-type T
typeMod(real f, FirType T)
{
  switch (T) {
    case TypeI  : default: return 1;
    case TypeII : return cosr(f *(PIr*.5)); // Has zero at f=1.
    case TypeIII: return sinr(f * PIr);     // Has zero at f=0, f=1.
    case TypeIV : return sinr(f *(PIr*.5)); // Has zero at f=0.
  }
}



// Call the response fn.; optionally check its output for validity;
// apply typeMod if applicable:

static MmfirResult respFnMod(MmfirRespFn  respFn, FirType T,
    MmfirBandSpec const * b, int len, Point * p, int check, va_list args)
{
  va_list tmpArgs;
  MmfirResult err;
  va_copy(tmpArgs, args), err=respFn(b, len, p, tmpArgs), va_end(tmpArgs);
  if (err) return err;

  int i;

  if (check) {
    double   maxA=0;
    for (i=0; i<len; ++i) {
      double const d=fabs(p[i].a), w=p[i].w;
      maxA=maxD(maxA, d);
      if (d > 1e+10)
        return MmfirInvalidAmplitude;
      if (w < 1e-10 || w > 1e+10)
        return MmfirInvalidWeight;
    }
    if (maxA<1e-20)
      return MmfirInvalidAmplitude;
  }

  if (T) for (i=0; i<len; ++i) {
    double f;
    setF(p[i].f);
    real const t = typeMod(f,T);
    p[i].a=(double)(p[i].a/t);
    p[i].w=(double)(p[i].w*t/(b->weightF?f:1));
  }
  return MmfirSuccess;
}



// Calculate the estimated amplitude response at a single, given
// frequency using precalculated Lagrange interpolation coefficients.

static real A(LagrangeCoef const * coefs, int R, double f)
{
  int i = 0;
  real n = 0, d = 0, x = cosr(PIr*f), t;
  for (; i < R; d += t = coefs[i].beta / t, n += coefs[i].gamma * t, ++i)
    if (fabs((double)(t = x - coefs[i].x)) < epsilon)
      return coefs[i].gamma;
  return n/d;
}



static MmfirResult builtInRespFn( // For simple response, given in bandSpecs:
    MmfirBandSpec const * b, int len, Point * p, va_list args)
{
  (void)args;
  int i;
  double f;
  for (i=0; i<len; ++i) {
    setF(p[i].f); // '...R' values here contain pre-calculated gradient:
    p[i].w= (b->weightL + b->weightR*(f-b->freqL));
    p[i].a= (b->ampL    + b->ampR   *(f-b->freqL));
  }
  return MmfirSuccess;
}



static double minimality_(Point const * p, int len, double threshold)
{                      // Determine the minimality of the maximum error:
  double max=1+1e-9;
  while (len--) max=maxD(max, fabs(p[len].w));
  return -log10(max-1)-threshold;
}



MmfirReport mmfir(
    double            h[],
    MmfirFilterClass  filterClass,
    int               N0, // = filter order + 1
    int               numBands,
    MmfirBandSpec     bandSpecs[/* numBands */], // N.B. destroyed.
    double            accuracy,
    double            persistence,
    double            robustness,
    double            target,
    int               flags,
    MmfirState        * state,
    MmfirRespFn       respFn, ...)
{
  int      N = maxI(N0,2);
  FirType  T = (FirType)((filterClass!=MmfirSymmetric)*2 + !(N&1));
  int      R = N/2 + !T; // R = # of cosine functions; num peaks = R+1

  MmfirReport report = {MmfirSuccess, 0, 0, 0, 0, {0,0,0}};

  // For forward-compatibility of the use of the control
  // variables, (almost) silently ignore out-of-range values:
  int      density2 = iround((5+(numBands>2))*(1+(
             report.controls[0]=clamp(accuracy,0,7)     )));
  int      maxIterations = iround(128*pow(1.587,
             report.controls[1]=clamp(persistence,-3,3) ));
  int      density = 3+(numBands>2)+iround(
             report.controls[2]=clamp(robustness,0,3)   );

  int      doneInit=0, doLG=1; // For the lower/higher density passes.
  int      spaceLength=0, i, j, _1; // _1 holds +/- 1
  int      prevNumPeaks=0;
  double   f;
  double   delta_1 = 1e-30; // 1/delta
  double   totalWidth=0; // Of given bands; in which peaks will be distributed.
  MmfirBandSpec * b;

  // Allocate the working arrays:
  LagrangeCoef * coefs = calloc((size_t)(R+1), sizeof(*coefs));
  Point        * peaks = calloc((size_t)(maxPeaks), sizeof(*peaks));
  Point    * prevPeaks = calloc((size_t)(maxPeaks), sizeof(*prevPeaks));
  Point       * space0 = 0, * space = 0; // Need to determine its size first.

  MmfirResult err;
  va_list respFnArgs, tmpArgs;
  va_start(respFnArgs, respFn);
  if (!respFn)
    respFn = builtInRespFn;

  // Initial checks:
  if (numBands<=0)
    returnError(MmfirInvalidNumBands);
  if (N!=N0)
    returnError(MmfirInvalidOrder);
  if (bandL->freqL<0 || bandR->freqR>1)
    returnError(MmfirInvalidFrequency);

  // Further checks, and pre-calc. linear interp. gradients, etc.:
  for (b=bandL; b<=bandR; ++b) {
    double df = b->freqR - b->freqL;

    if ((df <= userEpsilon) || (b<bandR && b->freqR>(b+1)->freqL))
      returnError(MmfirInvalidFrequency);

    if (respFn==builtInRespFn) {
      if (filterClass==MmfirDifferentiator)
        b->ampL *= -1, b->ampR *= -1;

      b->weightF = filterClass==MmfirDifferentiator &&
        (fabs(b->ampL) > userEpsilon || fabs(b->ampR) > userEpsilon);
      if (b->weightF)
        b->weightL *= 2, b->weightR *= 2; // For 2/f weighting.

      b->ampR    = (b->ampR    - b->ampL   ) / df;
      b->weightR = (b->weightR - b->weightL) / df;

    }
    else b->weightF=0;
    totalWidth += df;
  }

  // Avoid dividing by typeMod(0,T)=0:
  if ((T == TypeIII || T == TypeIV) && (b=bandL)->freqL < (f=epsilon)) {
    Point p = {b->freqL,0,0,0};
    va_copy(tmpArgs, respFnArgs), err=respFn(b,1,&p,tmpArgs), va_end(tmpArgs);
    if (err) returnError(err);
    if (fabs(p.a) > userEpsilon)
      returnError(MmfirInvalidDcAmplitude);
    if (respFn==builtInRespFn) {
      b->ampL += b->ampR*(f-b->freqL);
      b->weightL += b->weightR*(f-b->freqL);
    }
    totalWidth -= f-b->freqL, b->freqL = f;
  }

  // Avoid dividing by typeMod(1,T)=0:
  if ((T==TypeII || T==TypeIII) && (b=bandR)->freqR > (f=1-epsilon)) {
    Point p = {b->freqR,0,0,0};
    va_copy(tmpArgs, respFnArgs), err=respFn(b,1,&p,tmpArgs), va_end(tmpArgs);
    if (err) returnError(err);
    if (fabs(p.a) > userEpsilon)
      returnError(MmfirInvalidNyquistAmplitude);

    totalWidth -= b->freqR-f, b->freqR = f;
  }

  // Separate any contiguous bands:
  for (b=bandL; b<bandR; ++b)
    if (b->freqR > (f = b[1].freqL-epsilon))
      totalWidth -= b->freqR-f, b->freqR = f;

  // Usually, start with a lower density, then switch to the final value:
  for (density2=maxI(density,density2);
      report.result==MmfirSuccess && density<=density2;
      doLG=0, density+=maxI(density2-density,1)) {

    // Map current density to a convergence minimality threshold:
    double minimalityThreshold = density<density2? density*.35-.45 :
      (density+35)*.05; // Though in this case, stability usually prevails.
    if (debugLevel>1)
      REprintf("minimalityThreshold: %g\n", minimalityThreshold);

    { // Allocate & populate the analysis frequency space:
      rand_t random = 0;

      for (b=bandL; b<=bandR; ++b) {
        double df = b->freqR - b->freqL;
        b->portion = b->portion0 = df/totalWidth;
        b->ends=2-(b->freqL <= epsilon || b->freqR >= 1-epsilon);
      }
      apportionSpace(bandL, numBands, R+1);

      // Determine # of points in analysis space:
      for (spaceLength=0, b=bandL; b<=bandR; b->endP=spaceLength, ++b) {
        double m = (R+1)*b->portion;
        if (b->freqL <= epsilon || b->freqR >= 1-epsilon)
          spaceLength += 1 + (int)ceil(m*density);
        else
          spaceLength += 2 * (1 + (int)ceil(m*.5*density)) - 1;
      }

      // Allocate space and calculate frequencies:
      free(space0);
      if ((space0 = (Point *)calloc((size_t)(spaceLength+2), sizeof(*space)))) {
        space = space0+1; // Simplifies peak detection at edges of space.
        for (j=0, b = bandL; b<=bandR; ++b) {
          double m = (R+1)*b->portion;
          if (b->freqL <= epsilon)
            j+= warpSpaceL(space+j, m, b->freqL, b->freqR, density, &random);
          else if (b->freqR >= 1-epsilon)
            j+= warpSpaceR(space+j, m, b->freqL, b->freqR, density, &random);
          else
            j+= warpSpaceRL(space+j, m, b->freqL, b->freqR, density, &random);
        }
        err = respFnMod(
            respFn, T, bandL, spaceLength, space, !doneInit, respFnArgs);
        if (err) returnError(err);
      }
    }

    if (!space0 || !coefs || !peaks || !prevPeaks)
      returnError(MmfirOutOfMemory);

    // Initial 'guess' distributes peaks evenly through warped space.  The
    // offset of half a step at each end facilitates longer filters:
    if (!doneInit++) for (i=0; i<=R; ++i)
      peaks[i].f=space[iround((i+.5)/(R+1)*(spaceLength-1))].f;

    report.result=MmfirOngoing;
    while (1) { // Perform the Remez exchange (until break, below):
      double max, minimality, stability=0;
      int numPeaks, ok, extras=0;

      if (doLG++) { // Calculate Lagrange coefficients for response estimate:
        real denom = 0, numer = 0, t1, t, delta;

        err=respFnMod(respFn, T, bandL, R+1, peaks, 0, respFnArgs);
        if (err) returnError(err);

        for (i=0; i<=R; ++i)
          coefs[i].x = cosr(PIr * peaks[i].f);

        for (_1 = -1, b=bandL, i=0; i <= R; ++i) {
          for (t=1, j=0; j<=R; t1=t, t*=2*(coefs[i].x-coefs[j].x) +(i==j), ++j);
          if (!t1 || !t)
            returnError(MmfirNumericalError);
          coefs[i].beta = 1/t1; // Value at i==R will not be used.
          numer += peaks[i].a/t, denom += (_1 = -_1)/(t*peaks[i].w);
        }
        if (!numer || !denom)
          returnError(MmfirNumericalError);
        delta_1 = 1/(double)(delta = numer/denom);

        for (_1 = -1, b=bandL, i=0; i < R; ++i)
          coefs[i].gamma = peaks[i].a-(_1 = -_1)*delta/peaks[i].w;
      }

      // Stop Remez here if converged or no iterations remain:
      if (report.result == MmfirSuccess)
        break;
      else if (report.iterations == maxIterations) {
        report.result = density<density2? MmfirGaveUp1 : MmfirGaveUp;
        break;
      }
      else ++report.iterations;

      // Evaluate the normalised (i.e. converges to [-1,1]) error function:
      for (i=0, b=bandL; b<=bandR; i=b->endP, ++b) {
        #ifdef _OPENMP
        size_t const iterations = (size_t)(b->endP - i);
        int const threads =
          ompThreadCount(iterations, iterations * (size_t)R);
        #pragma omp parallel for if(threads > 1) num_threads(threads)
        #endif
        for (j=i; j<b->endP; ++j) { // omp needs separate counter j here.
          Point * const p = space+j;
          p->e = delta_1 * (double)(p->w * (p->a - A(coefs, R, p->f)));
        }
      }
      report.FEs += spaceLength;

      // Find and store all local peaks in error magnitude:
      b=bandL, numPeaks=0;
      for (i=0; i<spaceLength; ++i) if (dir(e(i+1)-e(i))!=dir(e(i)-e(i-1))) {
        double e, A, B, C, fp;

        // Approx. continuous-space peak (by interpolation in discrete-space):
        setF(f(i));
        j = i - (f(i) == b->freqR) + (f(i) == b->freqL); // Band-edge back-off
        poly2(f(j-1), f(j), f(j+1), e(j-1), e(j), e(j+1), &A, &B, &C);
        fp = B/A*-.5; // Freq. at which 1st derivative of poly is 0.
        ok = j==i || (dir(e(j+1)-e(j))*dir(e(j)-e(j-1)) > 0 &&
            fp >= f(i-(j<=i)) && fp <= f(i+(j>=i))); // Band-edge needs care.
        f = ok? fp : f, e = ok? poly2Val(A,B,C,f) : e(i);

        // Store the peak:                             But avoid twin-    \_
        if (!numPeaks || f>peaks[numPeaks-1].f) {   /* peaks in this case:  \ */
          if (numPeaks == maxPeaks) // Likely due to numerical error.
            returnError(MmfirTooManyPeaks);
          peaks[numPeaks].w=peaks[numPeaks].e=e, peaks[numPeaks++].f=f; // (1)
        }
      }

      if (debugLevel>2) { // Dump arrays to files for debugging:
        char name[9];
        FILE * f = fopen(name,"w");
        for (i=0; i<spaceLength;++i) fprintf(f, "%.16f %.9g\n", f(i), e(i));
        fclose(f);
        name[0]='p';
	f = fopen(name,"w");
        for (i=0; i<numPeaks;++i) fprintf(f, "%.16f % .9g\n", F(i), E(i));
        fclose(f);
      }

      // Check that there are at least R+1 peaks:
      if ((extras = numPeaks - (R+1)) < 0)
        returnError(MmfirTooFewPeaks); // Perhaps due to numerical error.

      // Before discarding any peaks, check freq. stability of the entire set:
      if (numPeaks == prevNumPeaks) {
        for (max=1e-9, i=0; i<numPeaks-1; ++i)
          max=maxD(max,fabs(1-(prevPeaks[i+1].f-prevPeaks[i].f)/(F(i+1)-F(i))));
        stability = -log10(max);
      }
      else stability = -9.99;

      // If converged, don't change peak-set:
      if (extras && density==density2 && hasConverged(numPeaks))
        for (i=0; i<numPeaks; peaks[i].e=prevPeaks[i].w, ++i); // Per (1) above
      memcpy(prevPeaks, peaks, (unsigned)numPeaks*sizeof(*peaks));
      prevNumPeaks = numPeaks;

      // If >R+1 peaks, reduce to R+1 by discarding lesser peaks/pairs:
      while (numPeaks>R+1) {
        int n=1, try2, d_na, d; // d = index of peak to discard.

        // Find a lesser peak, either overall or at one end:
        try2=numPeaks-(R+1)>1, i=try2? 1 : numPeaks-1;
        for (d_na=d=0; !d_na && i<numPeaks; ++i) {
          d_na = E(i)*E(i-1)>0? i : d_na; // Check alternating.
          d = fabs(E(i)) < fabs(E(d))? i : d; // Check magnitude.
        }
        d=d_na? d_na: d; // Prefer to discard non-alternating.
        if (try2 && d && d != numPeaks-1) // Discard pair only if not at end.
          n=2, d-=fabs(E(d-1))<fabs(E(d+1)); // Choose lesser adjacent.
        memmove(peaks+d, peaks+d+n, (size_t)((numPeaks-=n)-d)*sizeof(*peaks));
        if (debugLevel > 1) REprintf("x%i %i:%i\n", d_na, d, n);
      }

      if (hasConverged(R+1))
        report.result = MmfirSuccess;

      if (debugLevel) REprintf("%2i %2i %2i % .2f % .2f\n",
        report.iterations, density, extras, maxD(minimality,-9.99), stability);
    }
  }

  if (target && 1/fabs(delta_1) > target)
    returnError(MmfirMissedTarget);

  { // Generate filter coefficients:
    double * a = (double *)prevPeaks; // Reuse, as it's no longer in use.
    double _1 = T>TypeII? -1:1, phi = (T<TypeIII)*.5;

    // Sample the final estimated response; modify for filter type:
    #ifdef _OPENMP
    size_t const sampleIterations = (size_t)(N/2 + 1);
    int const sampleThreads =
      ompThreadCount(sampleIterations, sampleIterations * (size_t)R);
    #pragma omp parallel for if(sampleThreads > 1) num_threads(sampleThreads)
    #endif
    for (i=0; i<=N/2; ++i) {
      double f=2.*i/N;
      a[i] = (double)(A(coefs, R, f) * typeMod(f, T));
    }
    report.FEs += N/2+1;

    // --> time-domain using symmetry-aware IDFT (could also use IFFT):
    a[N/2] /= 1+T/3;
    #ifdef _OPENMP
    size_t const idftIterations = (size_t)((N+1)/2);
    int const idftThreads =
      ompThreadCount(idftIterations, idftIterations * (size_t)(N/2));
    #pragma omp parallel for if(idftThreads > 1) num_threads(idftThreads)
    #endif
    for (i=0; i <(N+1)/2; ++i) {
      double s;
      int j;
      for (s=*a*.5, j=1; j <= N/2; s += sin(PI*((N-1-2.*i)/N*j+phi))*a[j], ++j);
      h[i] = _1*(h[N-1-i] = 2*s/N);
    }
  }

END:
  if (state) { // Needed only for Octave/Matlab:
    free(space0), space0=0;
    int s0=iround(totalWidth*R*pow(4,2+accuracy)), s=0;
    for (j=0, b = bandL; b<=bandR; ++b)
      s += maxI(2, iround(b->portion*s0));
    int S=s+numBands;
    if ((space = (Point *)calloc((size_t)S, sizeof(*space)))) {
      for (j=0, b = bandL; b<=bandR; ++b) {
        int m = maxI(2, iround(b->portion*s0));
        for (i=0; i<m; ++i)
          space[j++].f=b->freqL+(b->freqR-b->freqL)/m*i;
        space[j++].f=b->freqR;
      }
      va_copy(tmpArgs,respFnArgs),respFn(bandL,S,space,tmpArgs),va_end(tmpArgs);
    }
    state->space = space, space=0;
    state->peaks = peaks, peaks=0;
    state->spaceLength=S;
    state->peaksLength=R+1;
  }
  va_end(respFnArgs);
  free(space0);
  free(prevPeaks);
  free(peaks);
  free(coefs);

  // Complete and return the report:
  report.text = report.result >= MmfirInvocationError?
    errorText[report.result-MmfirInvocationError] :
    report.result >= MmfirError?  "numerical-error" :
    report.result >= MmfirWarning?  "not-converged" : "success";
  report.minimax = 1/fabs(delta_1);
  return report;
}
