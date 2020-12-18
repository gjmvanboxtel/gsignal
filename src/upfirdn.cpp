// upfirdn.cpp
// Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
// Original Octave code:
// Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France) <ecm@apc.univ-paris7.fr>
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
// 20200929  GvB       setup for gsignal v0.1.0
//---------------------------------------------------------------------------------------------------------------------


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix upfirdn (NumericMatrix x, NumericMatrix h, int p, int q)
{

  int rx = x.nrow();
  int cx = x.ncol();
  // assume h has same number of columns as x (check in R code)

  int Lh = h.nrow();
  double r = p/(double(q));
  int Ly = ceil (double((rx-1)*p + Lh) / double(q));

  NumericMatrix y(Ly, cx);

  for (int c = 0; c < cx; c++)
  {
    int m = 0;
    while (m < Ly)
    {
      int n = floor (m/r);
      int lm = (m * q) % p;
      int k = 0;
      double accum = 0.0;
      do
      {
        int ix = n - k;
        if (ix >= rx)
        {
          k ++;
          continue;
        }

        int ih = k * p + lm;
        if ((ih >= Lh) | (ix < 0))
          break;

        accum += h(ih, c) * x (ix, c);
        k++;
      }
      while (1);

      y (m, c) = accum;
      m ++;
    }

  }

  return y;
}

