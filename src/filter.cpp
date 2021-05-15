// filter.cpp
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
// 20210317  GvB       setup for gsignal v0.3.0 (real-valued filter)
// 20210515  GvB       bugfix: resize a and b vectors; check length zi
//---------------------------------------------------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rfilter (NumericVector b, NumericVector a, NumericVector x, NumericVector zi) {
  
  int lb = b.length();
  int la = a.length();
  int lx = x.length();
  int lzi = zi.length();
  int lab = (la > lb ? la : lb);
  if (lzi != lab - 1) {
    return (R_NilValue);
  }
  NumericVector bb(lab), aa(lab);
  for (int i = 0; i < lb; i++) bb(i) = b(i);
  for (int i = 0; i < la; i++) aa(i) = a(i);

  NumericVector y(lx);
  
  if (la > 1) {
    for (int i = 0; i < lx; i++) {
      y(i) = zi(0) + bb(0) * x(i);
      if (lzi > 0) {
        for (int j = 0; j < lzi - 1; j++) {
          zi(j) = zi(j + 1) - aa(j + 1) * y(i) + bb(j + 1) * x(i);
        }
        zi(lzi - 1) = bb(lzi) * x(i) - aa(lzi) * y(i);
      } else {
        zi(0) = bb(lzi) * x(i) - aa(lzi) * y(i);
      }
    }
  } else if (lzi > 0) {
    for (int i = 0; i < lx; i++) {
      y(i) = zi(0) + bb(0) * x(i);
      if (lzi > 1) {
        for (int j = 0; j < lzi - 1; j++) {
          zi(j) = zi(j + 1) + bb(j + 1) * x(i);
        }
        zi(lzi - 1) = bb(lzi) * x(i);
      } else {
        zi(0) = bb(1) * x(i);
      }
    }
  }
  
  List L = List::create(_("y") = y, _("zf") = zi);
  return (L);
}
