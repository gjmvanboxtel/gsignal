// sosfilt.cpp
// Copyright (C) 2021 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
// 20200413   GvB     setup for gsignal v0.1.0
// 20210328   GvB     v0.3.0; used Python setup to handle initial conditions
//---------------------------------------------------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rsosfilt (NumericMatrix sos, NumericVector x, NumericMatrix zi) {
  
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
  int nZiCol = zi.ncol();
  if (nZiCol != 2) {
    return R_NilValue;
  }
  int nZiRow = zi.nrow();
  if (nZiRow != nSections) {
    return R_NilValue;
  }
  
  NumericVector y = x;
  NumericVector zf = clone(zi);
  double yi = 0.0;
  for (int iSamp = 0; iSamp < nSamp; iSamp++) {
    for (int iSection = 0; iSection < nSections; iSection++) {
      yi = y(iSamp);  // make a temporary copy
      //Use direct II transposed structure:
      y(iSamp) = sos(iSection, 0) * yi + zf(iSection, 0);
      zf(iSection, 0) = sos(iSection, 1) * yi - sos(iSection, 4) * y(iSamp) +  zf(iSection, 1);
      zf(iSection, 1) = sos(iSection, 2) * yi - sos(iSection, 5) * y(iSamp);
    }
  }
  
  List L = List::create(_["y"] = y, _["zf"] = zf);
  return (L);  
}
