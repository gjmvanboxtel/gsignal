// conv2.R
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
// 20200227  GvB       setup for gsignal v0.1.0
//---------------------------------------------------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix conv2df(NumericMatrix a, NumericMatrix b) {

  int aRows = a.nrow();
  int aCols = a.ncol();
  int bRows = b.nrow();
  int bCols = b.ncol();

  int yRows = aRows + bRows - 1;
  int yCols = aCols + bCols - 1;
  NumericMatrix y(yRows, yCols);

  for(int i=0; i < (yRows + bRows - 1); i++) {
    for(int j=0; j < (yCols + bCols - 1); j++) {
      for(int m=0; m < bRows; m++) {
        int im = i-m;
        for(int n=0; n < bCols; n++) {
          int jn = j-n;
          if (im >= 0 && jn >= 0 && im < aRows && jn < aCols) {
            y(i, j) += a(im, jn) * b(m, n);
          }
        }
      }
    }
  }
  return y;
}

// [[Rcpp::export]]
NumericMatrix conv2ds(NumericMatrix a, NumericMatrix b) {

  int aRows = a.nrow();
  int aCols = a.ncol();
  int bRows = b.nrow();
  int bCols = b.ncol();

  int yRows = aRows;
  int yCols = aCols;
  NumericMatrix y(yRows, yCols);

  int bCX = bCols / 2;
  int bCY = bRows / 2;

  for(int i=0; i < yRows; i++) {
    for(int j=0; j < yCols; j++) {
      for(int m=0; m < bRows; m++) {
        int mm = bRows - 1 - m;
        for(int n=0; n < bCols; n++) {
          int nn = bCols - 1 - n;
          int ii = i + (bCY - mm);
          int jj = j + (bCX - nn);
          if( ii >= 0 && ii < yRows && jj >= 0 && jj < yCols )
            y(i, j) += a(ii, jj) * b(mm, nn);
        }
      }
    }
  }
  return y;
}

// [[Rcpp::export]]
NumericMatrix conv2dv(NumericMatrix a, NumericMatrix b) {

  int aRows = a.nrow();
  int aCols = a.ncol();
  int bRows = b.nrow();
  int bCols = b.ncol();

  int yRows = aRows - bRows + 1;
  int yCols = aCols - bCols + 1;
  NumericMatrix y(yRows, yCols);

  for(int i=0; i < yRows; i++) {
    for(int j=0; j < yCols; j++) {
      for(int m=0; m < bRows; m++) {
        int im = i+m;
        for(int n=0; n < bCols; n++) {
          int jn = j+n;
          int br = bRows - m - 1;
          int bc = bCols - n - 1;
          y(i, j) += a(im, jn) * b(br, bc);
          }
      }
    }
  }
  return y;
}
