// remez.cpp
// Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// Version history
// 20200803 GvB       setup for gsignal v0.1.0
// 20261111 GvB       exit on error on Search instead of returning the error
//                    value because R expects a NumericVector not an integer
// 20261111 GvB       corrected several bugs marked by GvB 20261111
//-----------------------------------------------------------------------------

// GvB 20261111
// Evaluate the desired response at the actual grid frequency.
//
// The original Octave implementation interpolated by grid index,
// assuming equally spaced frequencies between the band endpoints.
// However, the grid uses a fixed step delf and explicitly adjusts
// its endpoints, so index-based interpolation can assign incorrect
// desired amplitudes to grid frequencies.
//
// This is particularly important for differentiators and other
// filters with nonconstant desired responses.

#include <Rcpp.h>
using namespace Rcpp;

// GvB 20261111
#include <cmath>
#include <limits>
#include <vector>

/*

 Copyright (C) 1995, 1998 Jake Janovetz <janovetz@uiuc.edu>
 Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
 Copyright (C) 2000 Kai Habel <kahacjde@linux.zrz.tu-berlin.de>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; see the file COPYING.  If not, see
 <https://www.gnu.org/licenses/>.

 */

/**************************************************************************
 *  There appear to be some problems with the routine Search. See comments
 *  therein [search for PAK:].  I haven't looked closely at the rest
 *  of the code---it may also have some problems.
 *************************************************************************/
// GvB 20261111 the present version tried to address those problems


#define CONST const
#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define NEGATIVE       0
#define POSITIVE       1

#define Pi             3.1415926535897932
#define Pi2            6.2831853071795865

#define GRIDDENSITY    16
#define MAXITERATIONS  40


// GvB 20261111 added DensegridLayout structure and GetDensegridLayoutfunction
struct DenseGridLayout
{
  std::vector<int> count;
  int size;
  double delf;
  double grid0;
};


DenseGridLayout GetDenseGridLayout(
    int r,
    int numband,
    const NumericVector& bands,
    int symmetry,
    int griddensity)
{
  if (r <= 0 || griddensity <= 0 || numband <= 0)
    Rcpp::stop("Invalid dense-grid parameters");
  
  if (bands.size() != 2 * static_cast<R_xlen_t>(numband))
    Rcpp::stop("Invalid number of band edges");
  
  DenseGridLayout layout;
  
  layout.count.resize(numband);
  layout.size = 0;
  
  layout.delf =
    0.5 / (static_cast<double>(griddensity) * r);
  
  if (!std::isfinite(layout.delf) || layout.delf <= 0.0)
    Rcpp::stop("Invalid dense-grid spacing");
  
  for (int i = 0; i < 2 * numband; ++i)
  {
    if (!std::isfinite(bands[i]) ||
        bands[i] < 0.0 ||
        bands[i] > 0.5)
    {
      Rcpp::stop("Invalid frequency band edge");
    }
  }
  
  layout.grid0 =
    (symmetry == NEGATIVE && layout.delf > bands[0])
    ? layout.delf
  : bands[0];
  
  for (int band = 0; band < numband; ++band)
  {
    const double lowf =
      (band == 0)
    ? layout.grid0
    : bands[2 * band];
    
    const double highf = bands[2 * band + 1];
    
    if (highf <= lowf)
    {
      Rcpp::stop(
        "Dense-grid band %d has invalid width",
        band + 1
      );
    }
    
    // Preserve the original CreateDenseGrid formula.
    const double raw =
      (highf - lowf) / layout.delf + 0.5;
    
    if (!std::isfinite(raw))
      Rcpp::stop("Nonfinite dense-grid size");
    
    if (raw > static_cast<double>(
        std::numeric_limits<int>::max()))
    {
      Rcpp::stop("Dense-grid size exceeds integer range");
    }
    
    int k = static_cast<int>(raw);
    
    // Preserve the original antisymmetric adjustment.
    if (band == 0 && symmetry == NEGATIVE)
      --k;
    
    // The existing interpolation requires k >= 2.
    if (k < 2)
    {
      Rcpp::stop(
        "Dense-grid band %d contains fewer than "
        "two grid points; increase density or "
        "band width",
        band + 1
      );
    }
    
    if (k > std::numeric_limits<int>::max() -
        layout.size)
    {
      Rcpp::stop("Total dense-grid size exceeds integer range");
    }
    
    layout.count[band] = k;
    layout.size += k;
  }
  
  if (layout.size < r + 1)
  {
    Rcpp::stop(
      "Dense grid has %d points, but at least %d "
      "are required",
      layout.size, r + 1
    );
  }
  
  return layout;
}


/*******************
 * CreateDenseGrid
 *=================
 * Creates the dense grid of frequencies from the specified bands.
 * Also creates the Desired Frequency Response function (D[]) and
 * the Weight function (W[]) on that dense grid
 *
 *
 * INPUT:
 * ------
 * int      r        - 1/2 the number of filter coefficients
 * int      numtaps  - Number of taps in the resulting filter
 * int      numband  - Number of bands in user specification
 * double   bands[]  - User-specified band edges [2*numband]
 * double   des[]    - Desired response per band [2*numband]
 * double   weight[] - Weight per band [numband]
 * int      symmetry - Symmetry of filter - used for grid check
 * int      griddensity
 *
 * OUTPUT:
 * -------
 * int    gridsize   - Number of elements in the dense frequency grid
 * double Grid[]     - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]        - Desired response on the dense grid [gridsize]
 * double W[]        - Weight function on the dense grid [gridsize]
 *******************/

// GvB 20261111 replaced original with corrected CreateDenseGrid function

void CreateDenseGrid(
    int numtaps,
    int numband,
    const NumericVector& bands,
    const NumericVector& des,
    const NumericVector& weight,
    const DenseGridLayout& layout,
    NumericVector& Grid,
    NumericVector& D,
    NumericVector& W,
    int symmetry)
{
  const int gridsize = layout.size;
  
  // Validate output vector lengths.
  if (Grid.size() != gridsize ||
      D.size() != gridsize ||
      W.size() != gridsize)
  {
    Rcpp::stop("Dense-grid vector length mismatch");
  }
  
  // Validate input vector lengths.
  if (bands.size() != 2 * static_cast<R_xlen_t>(numband) ||
      des.size() != 2 * static_cast<R_xlen_t>(numband) ||
      weight.size() != numband)
  {
    Rcpp::stop("Invalid band, desired-response or weight length");
  }
  
  // Validate layout.
  if (static_cast<int>(layout.count.size()) != numband ||
      gridsize <= 0)
  {
    Rcpp::stop("Invalid dense-grid layout");
  }
  
  int j = 0;
  
  /*
   * STEP 1:
   * Construct the complete frequency grid and weights.
   *
   * Do not calculate desired responses yet, because
   * the frequency grid may subsequently be adjusted.
   */
  for (int band = 0; band < numband; ++band)
  {
    double lowf =
      (band == 0)
    ? layout.grid0
    : bands[2 * band];
    
    const double highf = bands[2 * band + 1];
    
    const int k = layout.count[band];
    
    if (k < 2)
    {
      Rcpp::stop(
        "Dense-grid band %d has fewer than two points",
        band + 1
      );
    }
    
    // Check that the current band fits into the allocated grid.
    if (k > gridsize - j)
    {
      Rcpp::stop("Dense-grid construction overflow");
    }
    
    for (int i = 0; i < k; ++i)
    {
      Grid[j] = lowf;
      W[j] = weight[band];
      
      lowf += layout.delf;
      ++j;
    }
    
    /*
     * Preserve the original behavior:
     * force the last point of each band to its exact
     * upper frequency endpoint.
     */
    Grid[j - 1] = highf;
  }
  
  // Verify that all allocated grid points were generated.
  if (j != gridsize)
  {
    Rcpp::stop(
      "Dense-grid construction produced %d "
      "points; expected %d",
      j, gridsize
    );
  }
  
  /*
   * STEP 2:
   * Preserve the existing Type III Nyquist adjustment.
   *
   * Type III filters have odd numtaps and negative
   * symmetry. Their response must vanish at Nyquist.
   */
  if (symmetry == NEGATIVE &&
      (numtaps % 2) &&
      Grid[gridsize - 1] > (0.5 - layout.delf))
  {
    Grid[gridsize - 1] = 0.5 - layout.delf;
  }
  
  /*
   * STEP 3:
   * Calculate the desired response from the ACTUAL
   * grid frequencies, including endpoint adjustments.
   *
   * Linear interpolation within each band:
   *
   * D(f) = a0 + (f - f0) * (a1 - a0) / (f1 - f0)
   *
   * Here f0 and f1 are the original band endpoints,
   * not the shifted first grid frequency.
   */
  j = 0;
  
  for (int band = 0; band < numband; ++band)
  {
    const double f0 = bands[2 * band];
    const double f1 = bands[2 * band + 1];
    
    const double a0 = des[2 * band];
    const double a1 = des[2 * band + 1];
    
    const int k = layout.count[band];
    
    const double bandwidth = f1 - f0;
    
    if (!(bandwidth > 0.0))
    {
      Rcpp::stop(
        "Dense-grid band %d has invalid frequency width",
        band + 1
      );
    }
    
    const double slope = (a1 - a0) / bandwidth;
    
    for (int i = 0; i < k; ++i)
    {
      D[j] = a0 + (Grid[j] - f0) * slope;
      ++j;
    }
  }
}
   
/********************
 * InitialGuess
 *==============
 * Places Extremal Frequencies evenly throughout the dense grid.
 *
 *
 * INPUT:
 * ------
 * int r        - 1/2 the number of filter coefficients
 * int gridsize - Number of elements in the dense frequency grid
 *
 * OUTPUT:
 * -------
 * int Ext[]    - Extremal indexes to dense frequency grid [r+1]
 ********************/

void InitialGuess(int r, NumericVector& Ext, int gridsize)
{
   int i;

   for (i=0; i<=r; i++)
      Ext[i] = i * (gridsize-1) / r;
}


/***********************
 * CalcParms
 *===========
 *
 *
 * INPUT:
 * ------
 * int    r      - 1/2 the number of filter coefficients
 * int    Ext[]  - Extremal indexes to dense frequency grid [r+1]
 * double Grid[] - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double ad[]   - 'b' in Oppenheim & Schafer [r+1]
 * double x[]    - [r+1]
 * double y[]    - 'C' in Oppenheim & Schafer [r+1]
 ***********************/

void CalcParms(int r, NumericVector Ext, NumericVector Grid, NumericVector D, NumericVector W,
               NumericVector& ad, NumericVector& x, NumericVector& y)
{
   int i, j, k, ld;
   double sign, xi, delta, denom, numer;

   /*
    * Find x[]
    */
   for (i=0; i<=r; i++)
      x[i] = cos(Pi2 * Grid[Ext[i]]);

   /*
    * Calculate ad[]  - Oppenheim & Schafer eq 7.132
    */
   ld = (r-1)/15 + 1;         /* Skips around to avoid round errors */
   for (i=0; i<=r; i++)
   {
      denom = 1.0;
      xi = x[i];
      for (j=0; j<ld; j++)
      {
         for (k=j; k<=r; k+=ld)
            if (k != i)
               denom *= 2.0*(xi - x[k]);
      }
      if (fabs(denom)<0.00001)
         denom = 0.00001;
      ad[i] = 1.0/denom;
   }

   /*
    * Calculate delta  - Oppenheim & Schafer eq 7.131
    */
   numer = denom = 0;
   sign = 1;
   for (i=0; i<=r; i++)
   {
      numer += ad[i] * D[Ext[i]];
      denom += sign * ad[i]/W[Ext[i]];
      sign = -sign;
   }
   delta = numer/denom;
   sign = 1;

   /*
    * Calculate y[]  - Oppenheim & Schafer eq 7.133b
    */
   for (i=0; i<=r; i++)
   {
      y[i] = D[Ext[i]] - sign * delta/W[Ext[i]];
      sign = -sign;
   }
}


/*********************
 * ComputeA
 *==========
 * Using values calculated in CalcParms, ComputeA calculates the
 * actual filter response at a given frequency (freq).  Uses
 * eq 7.133a from Oppenheim & Schafer.
 *
 *
 * INPUT:
 * ------
 * double freq - Frequency (0 to 0.5) at which to calculate A
 * int    r    - 1/2 the number of filter coefficients
 * double ad[] - 'b' in Oppenheim & Schafer [r+1]
 * double x[]  - [r+1]
 * double y[]  - 'C' in Oppenheim & Schafer [r+1]
 *
 * OUTPUT:
 * -------
 * Returns double value of A[freq]
 *********************/

double ComputeA(double freq, int r, NumericVector ad, NumericVector x, NumericVector y)
{
   int i;
   double xc, c, denom, numer;

   denom = numer = 0;
   xc = cos(Pi2 * freq);
   for (i=0; i<=r; i++)
   {
      c = xc - x[i];
      if (fabs(c) < 1.0e-7)
      {
         numer = y[i];
         denom = 1;
         break;
      }
      c = ad[i]/c;
      denom += c;
      numer += c*y[i];
   }
   return numer/denom;
}


/************************
 * CalcError
 *===========
 * Calculates the Error function from the desired frequency response
 * on the dense grid (D[]), the weight function on the dense grid (W[]),
 * and the present response calculation (A[])
 *
 *
 * INPUT:
 * ------
 * int    r      - 1/2 the number of filter coefficients
 * double ad[]   - [r+1]
 * double x[]    - [r+1]
 * double y[]    - [r+1]
 * int gridsize  - Number of elements in the dense frequency grid
 * double Grid[] - Frequencies on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the desnse grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double E[]    - Error function on dense grid [gridsize]
 ************************/

void CalcError(int r, NumericVector ad, NumericVector x, NumericVector y,
               int gridsize, NumericVector Grid, NumericVector D, NumericVector W,
               NumericVector& E)
{
   int i;
   double A;

   for (i=0; i<gridsize; i++)
   {
      A = ComputeA(Grid[i], r, ad, x, y);
      E[i] = W[i] * (D[i] - A);
   }
}

/************************
 * Search
 *========
 * Searches for the maxima/minima of the error curve.  If more than
 * r+1 extrema are found, it uses the following heuristic (thanks
 * Chris Hanson):
 * 1) Adjacent non-alternating extrema deleted first.
 * 2) If there are more than one excess extrema, delete the
 *    one with the smallest error.  This will create a non-alternation
 *    condition that is fixed by 1).
 * 3) If there is exactly one excess extremum, delete the smaller
 *    of the first/last extremum
 *
 *
 * INPUT:
 * ------
 * int    r        - 1/2 the number of filter coefficients
 * int    Ext[]    - Indexes to Grid[] of extremal frequencies [r+1]
 * int    gridsize - Number of elements in the dense frequency grid
 * double E[]      - Array of error values.  [gridsize]
 * OUTPUT:
 * -------
 * int    Ext[]    - New indexes to extremal frequencies [r+1]
 ************************/
int Search(int r, NumericVector& Ext, int gridsize, NumericVector E)
{
   int i, j, k, l, extra;     /* Counters */
   int up, alt;
//   int *foundExt;             /* Array of found extremals */

   /*
    * Allocate enough space for found extremals.
    */
   //foundExt = (int *)malloc((2*r) * sizeof(int));
   NumericVector foundExt(2*r);
   k = 0;

   /*
   * Check for extremum at 0.
   */
   if (((E[0]>0.0) && (E[0]>E[1])) ||
       ((E[0]<0.0) && (E[0]<E[1])))
      foundExt[k++] = 0;

   /*
    * Check for extrema inside dense grid
    */
   for (i=1; i<gridsize-1; i++)
   {
      if (((E[i]>=E[i-1]) && (E[i]>E[i+1]) && (E[i]>0.0)) ||
          ((E[i]<=E[i-1]) && (E[i]<E[i+1]) && (E[i]<0.0))) {
         // PAK: we sometimes get too many extremal frequencies
         if (k >= 2*r) return -3;
         foundExt[k++] = i;
      }
   }

   /*
    * Check for extremum at 0.5
    */
   j = gridsize-1;
   if (((E[j]>0.0) && (E[j]>E[j-1])) ||
       ((E[j]<0.0) && (E[j]<E[j-1]))) {
      if (k >= 2*r) return -3;
      foundExt[k++] = j;
   }

   // PAK: we sometimes get not enough extremal frequencies
   if (k < r+1) return -2;


   /*
    * Remove extra extremals
    */
   extra = k - (r+1);
   assert(extra >= 0);

   while (extra > 0)
   {
      if (E[foundExt[0]] > 0.0)
         up = 1;                /* first one is a maxima */
      else
         up = 0;                /* first one is a minima */

      l=0;
      alt = 1;
      for (j=1; j<k; j++)
      {
         if (fabs(E[foundExt[j]]) < fabs(E[foundExt[l]]))
            l = j;               /* new smallest error. */
         if ((up) && (E[foundExt[j]] < 0.0))
            up = 0;             /* switch to a minima */
         else if ((!up) && (E[foundExt[j]] > 0.0))
            up = 1;             /* switch to a maxima */
         else
         {
            alt = 0;
            // PAK: break now and you will delete the smallest overall
            // extremal.  If you want to delete the smallest of the
            // pair of non-alternating extremals, then you must do:
            //
            // if (fabs(E[foundExt[j]]) < fabs(E[foundExt[j-1]])) l=j;
            // else l=j-1;
            break;              /* Ooops, found two non-alternating */
         }                     /* extrema.  Delete smallest of them */
      }  /* if the loop finishes, all extrema are alternating */

      /*
       * If there's only one extremal and all are alternating,
       * delete the smallest of the first/last extremals.
       */
      if ((alt) && (extra == 1))
      {
         if (fabs(E[foundExt[k-1]]) < fabs(E[foundExt[0]]))
            /* Delete last extremal */
            l = k-1;
         // PAK: changed from l = foundExt[k-1];
         else
            /* Delete first extremal */
            l = 0;
         // PAK: changed from l = foundExt[0];
      }

      for (j=l; j<k-1; j++)        /* Loop that does the deletion */
      {
         foundExt[j] = foundExt[j+1];
         assert(foundExt[j]<gridsize);
      }
      k--;
      extra--;
   }

   for (i=0; i<=r; i++)
   {
      assert(foundExt[i]<gridsize);
      Ext[i] = foundExt[i];       /* Copy found extremals to Ext[] */
   }
   return 0;
}


/*********************
 * FreqSample
 *============
 * Simple frequency sampling algorithm to determine the impulse
 * response h[] from A's found in ComputeA
 *
 *
 * INPUT:
 * ------
 * int      N        - Number of filter coefficients
 * double   A[]      - Sample points of desired response [N/2]
 * int      symmetry - Symmetry of desired filter
 *
 * OUTPUT:
 * -------
 * double h[] - Impulse Response of final filter [N]
 *********************/
void FreqSample(int N, NumericVector A, NumericVector& h, int symm)
{
   int n, k;
   double x, val, M;

   M = (N-1.0)/2.0;
   if (symm == POSITIVE)
   {
      if (N%2)
      {
         for (n=0; n<N; n++)
         {
            val = A[0];
            x = Pi2 * (n - M)/N;
            for (k=1; k<=M; k++)
               val += 2.0 * A[k] * cos(x*k);
            h[n] = val/N;
         }
      }
      else
      {
         for (n=0; n<N; n++)
         {
            val = A[0];
            x = Pi2 * (n - M)/N;
            for (k=1; k<=(N/2-1); k++)
               val += 2.0 * A[k] * cos(x*k);
            h[n] = val/N;
         }
      }
   }
   else
   {
      if (N%2)
      {
         for (n=0; n<N; n++)
         {
            val = 0;
            x = Pi2 * (n - M)/N;
            for (k=1; k<=M; k++)
               val += 2.0 * A[k] * sin(x*k);
            h[n] = val/N;
         }
      }
      else
      {
         for (n=0; n<N; n++)
         {
            val = A[N/2] * sin(Pi * (n - M));
            x = Pi2 * (n - M)/N;
            for (k=1; k<=(N/2-1); k++)
               val += 2.0 * A[k] * sin(x*k);
            h[n] = val/N;
         }
      }
   }
}

/*******************
 * isDone
 *========
 * Checks to see if the error function is small enough to consider
 * the result to have converged.
 *
 * INPUT:
 * ------
 * int    r     - 1/2 the number of filter coefficients
 * int    Ext[] - Indexes to extremal frequencies [r+1]
 * double E[]   - Error function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * Returns 1 if the result converged
 * Returns 0 if the result has not converged
 ********************/

int isDone(int r, NumericVector Ext, NumericVector E)
{
   int i;
   double min, max, current;

   min = max = fabs(E[Ext[0]]);
   for (i=1; i<=r; i++)
   {
      current = fabs(E[Ext[i]]);
      if (current < min)
         min = current;
      if (current > max)
         max = current;
   }
   return (((max-min)/max) < 0.0001);
}

/********************
 * remez
 *=======
 * Calculates the optimal (in the Chebyshev/minimax sense)
 * FIR filter impulse response given a set of band edges,
 * the desired response on those bands, and the weight given to
 * the error in those bands.
 *
 * INPUT:
 * ------
 * int     numtaps     - Number of filter coefficients
 * int     numband     - Number of bands in filter specification
 * double  bands[]     - User-specified band edges [2 * numband]
 * double  des[]       - User-specified band responses [numband]
 * double  weight[]    - User-specified error weights [numband]
 * int     type        - Type of filter
 *
 * OUTPUT:
 * -------
 * double h[]      - Impulse response of final filter [numtaps]
 * returns         - true on success, false on failure to converge
 ********************/

// GvB 20261111: replaced allocation and construction sections

// [[Rcpp::export]]
NumericVector remez(NumericVector h, int numtaps,
         int numband, const NumericVector bands,
         const NumericVector des, const NumericVector weight,
         int type, int griddensity)
{
   int i, iter, gridsize, r, symmetry;
   double c;

   if (type == BANDPASS)
      symmetry = POSITIVE;
   else
      symmetry = NEGATIVE;

   r = numtaps/2;                  /* number of extrema */
   if ((numtaps%2) && (symmetry == POSITIVE))
      r++;

   /*
    * Determine exact dense-grid layout.
    * The layout contains the number of grid points per band,
    * the total grid size, and the grid spacing.
    */
   DenseGridLayout layout = GetDenseGridLayout(
     r,
     numband,
     bands,
     symmetry,
     griddensity
   );
   
   gridsize = layout.size;
   
   /*
    * Allocate exactly the required number of elements.
    */
   NumericVector Grid(gridsize), D(gridsize), W(gridsize), E(gridsize);
   NumericVector Ext(r+1), taps(r+1), x(r+1), y(r+1), ad(r+1);
   
   /*
    * Construct the dense frequency grid.
    */
   CreateDenseGrid(
     numtaps, numband, bands, des, weight,
     layout, Grid, D, W, symmetry
   );
   
   /*
    * Validate the generated dense grid.
    */
   for (int i = 0; i < gridsize; ++i)
   {
     if (!std::isfinite(Grid[i]) ||
         Grid[i] < 0.0 ||
         Grid[i] > 0.5)
     {
       Rcpp::stop("Invalid generated grid frequency");
     }
     
     if (!std::isfinite(D[i]) ||
         !std::isfinite(W[i]) ||
         W[i] <= 0.0)
     {
       Rcpp::stop("Invalid generated response or weight");
     }
   }
   
   /*
    * Verify that grid frequencies are strictly increasing.
    */
   for (int i = 1; i < gridsize; ++i)
   {
     if (!(Grid[i] > Grid[i - 1]))
     {
       Rcpp::stop(
         "Dense-grid frequencies are not strictly "
         "increasing at index %d", i
       );
     }
   }
   
   /*
    * Reject the singular Nyquist endpoint for Type II filters.
    */
   if (symmetry == POSITIVE &&
       numtaps % 2 == 0 &&
       Grid[gridsize - 1] >= 0.5)
   {
     Rcpp::stop(
       "Type II filter dense grid includes Nyquist, "
       "where the transformation is singular"
     );
   }
   
   InitialGuess(r, Ext, gridsize);
   /*
    * For Differentiator: (fix grid)
    */
   if (type == DIFFERENTIATOR)
   {
      for (i=0; i<gridsize; i++)
      {
         /* D[i] = D[i]*Grid[i]; */
         if (D[i] > 0.0001)
            W[i] = W[i]/Grid[i];
      }
   }

   /*
    * For odd or Negative symmetry filters, alter the
    * D[] and W[] according to Parks McClellan
    */
   if (symmetry == POSITIVE)
   {
      if (numtaps % 2 == 0)
      {
         for (i=0; i<gridsize; i++)
         {
            c = cos(Pi * Grid[i]);
            D[i] /= c;
            W[i] *= c;
         }
      }
   }
   else
   {
      if (numtaps % 2)
      {
         for (i=0; i<gridsize; i++)
         {
            c = sin(Pi2 * Grid[i]);
            D[i] /= c;
            W[i] *= c;
         }
      }
      else
      {
         for (i=0; i<gridsize; i++)
         {
            c = sin(Pi * Grid[i]);
            D[i] /= c;
            W[i] *= c;
         }
      }
   }

   /*
    * Perform the Remez Exchange algorithm
    */
   for (iter=0; iter<MAXITERATIONS; iter++)
   {
      CalcParms(r, Ext, Grid, D, W, ad, x, y);
      CalcError(r, ad, x, y, gridsize, Grid, D, W, E);
      int err = Search(r, Ext, gridsize, E);
      // GvB 20261111 if (err) return err;
      if (err != 0)
      {
        Rcpp::stop(
          "Remez Search failed: status = %d, r = %d, gridsize = %d",
          err, r, gridsize
        );
      }
      for(i=0; i <= r; i++) assert(Ext[i]<gridsize);
      if (isDone(r, Ext, E))
         break;
   }

   CalcParms(r, Ext, Grid, D, W, ad, x, y);

   /*
    * Find the 'taps' of the filter for use with Frequency
    * Sampling.  If odd or Negative symmetry, fix the taps
    * according to Parks McClellan
    */
   for (i=0; i<=numtaps/2; i++)
   {
      if (symmetry == POSITIVE)
      {
         if (numtaps%2)
            c = 1;
         else
            c = cos(Pi * (double)i/numtaps);
      }
      else
      {
         if (numtaps%2)
            c = sin(Pi2 * (double)i/numtaps);
         else
            c = sin(Pi * (double)i/numtaps);
      }
      taps[i] = ComputeA((double)i/numtaps, r, ad, x, y)*c;
   }

   /*
    * Frequency sampling design with calculated taps
   */
   FreqSample(numtaps, taps, h, symmetry);

   return iter<MAXITERATIONS?h:-1;
}
