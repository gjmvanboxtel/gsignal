// sosfilt.cpp
// Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
// Original Octave code:
// Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France) <ecm@apc.univ-paris7.fr>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
// See also: http://www.gnu.org/licenses/gpl-2.0.txt
//
// Version history
// 20200413  GvB       setup for gsignal v0.1.0
//---------------------------------------------------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sosfilt (NumericMatrix sos, NumericVector x) {
  
  // Invalid values should be caught by calling function, but just to be safe
  int nSections = sos.nrow();

  int nSosCol = sos.ncol();
  if (nSosCol != 6) {     
    return R_NilValue;
  }
  int nSamp = x.size();
  if (nSamp <= 0) {     
    return R_NilValue;
  }

  NumericVector y(nSamp);

  for (int iSection = 0; iSection < nSections; iSection++) {

    double v0 = 0.0, v1 = 0.0, v2 = 0.0;
    double a0, a1, a2, b0, b1, b2;
    
    a0 =   sos(iSection, 3);
    if (a0 == 0) {
      return rep(R_NaN, nSamp);
    } else {
      a1 = sos(iSection, 4) / a0;
      a2 = sos(iSection, 5) / a0;
      b0 = sos(iSection, 0) / a0;
      b1 = sos(iSection, 1) / a0;
      b2 = sos(iSection, 2) / a0;
    }
    for (int iSamp = 0; iSamp < nSamp; iSamp++) {
      v0 = x(iSamp) - (a1 * v1) - (a2 * v2);
      y(iSamp) = (b0 * v0) + (b1 * v1) + (b2 * v2);
      v2 = v1;
      v1 = v0;
    }
    x = y;
  }

  return y;  
}
