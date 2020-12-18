// fwht.cpp
// Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
// Original Octave code:
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
// 20201025  GvB       setup for gsignal v0.1.0
//---------------------------------------------------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// Based on pseudocode at
// https://en.wikipedia.org/wiki/Fast_Walsh%E2%80%93Hadamard_transform

// [[Rcpp::export]]
NumericMatrix fwht (NumericMatrix x) {

  int ncols = x.ncol();
  int n = x.nrow();
  NumericMatrix data = clone(x);

  for (int icol = 0; icol < ncols; icol++) {
    int h = 1;
    while (h < n) {
      for (int i = 0; i < n; i += h * 2) {
        for (int j = i; j < i + h; j++) {
          double xx = data(j, icol);
          double yy = data(j + h, icol);
          data(j, icol) = xx + yy;
          data(j + h, icol) = xx - yy;
        }
      }
      h *= 2;
    }
  }
  return (data);
}
