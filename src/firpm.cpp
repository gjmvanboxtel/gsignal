// firpm.cpp -- Rcpp interface to Octave Signal's mmfir algorithm.
// Copyright (C) 2026 Geert van Boxtel <gjmvanboxtel@gmail.com>
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
// NOTE: Keep the GPL license notices accompanying mmfir.c/mmfir.h in
// the package.
//
// Version history
// 20261111  GvB      setup for gsignal v0.4.0
//------------------------------------------------------------------------------

#include <Rcpp.h>
#include <vector>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <string>
#include "mmfir.h"

using namespace Rcpp;

// Release buffers transferred by mmfir to MmfirState, including on exceptions.
struct MmfirStateOwner {
  MmfirState state;
  MmfirStateOwner() : state{} {}
  ~MmfirStateOwner() {
    std::free(state.space);
    std::free(state.peaks);
  }
  MmfirStateOwner(const MmfirStateOwner&) = delete;
  MmfirStateOwner& operator=(const MmfirStateOwner&) = delete;
};


// R callback context stays alive on the stack throughout mmfir().
struct RResponseContext {
  SEXP function;
  SEXP extra;
  int order;
  SEXP frequencies;
  SEXP weights;
  bool failed;
  std::string error;
};

static MmfirResult rResponseCallback(
    const MmfirBandSpec* /*bands*/, int len,
    MmfirPoint points[], va_list args)
{
  RResponseContext* ctx = va_arg(args, RResponseContext*);
  if (!ctx || ctx->failed) return MmfirInvalidAmplitude;

  // No Rcpp exceptions or Rf_error() may escape through mmfir's C frames.
  int protected_count = 0;
  SEXP g = PROTECT(Rf_allocVector(REALSXP, len)); ++protected_count;
  for (int i = 0; i < len; ++i) REAL(g)[i] = points[i].f;
  SEXP order = PROTECT(Rf_ScalarInteger(ctx->order)); ++protected_count;
  SEXP call = PROTECT(Rf_lang5(ctx->function, order,
                               ctx->frequencies, g, ctx->weights));
  ++protected_count;
  SEXP tail = call;
  while (CDR(tail) != R_NilValue) tail = CDR(tail);

  const R_xlen_t nextra = XLENGTH(ctx->extra);
  SEXP names = Rf_getAttrib(ctx->extra, R_NamesSymbol);
  for (R_xlen_t i = 0; i < nextra; ++i) {
    SEXP node = PROTECT(Rf_cons(VECTOR_ELT(ctx->extra, i), R_NilValue));
    ++protected_count;
    if (names != R_NilValue && STRING_ELT(names, i) != NA_STRING &&
        CHAR(STRING_ELT(names, i))[0] != '\0')
      SET_TAG(node, Rf_installChar(STRING_ELT(names, i)));
    SETCDR(tail, node);
    tail = node;
  }

  int eval_error = 0;
  SEXP answer = PROTECT(R_tryEvalSilent(call, R_GlobalEnv, &eval_error));
  ++protected_count;
  if (eval_error) {
    ctx->failed = true;
    ctx->error = "R response function raised an error";
    UNPROTECT(protected_count);
    return MmfirInvalidAmplitude;
  }

  SEXP amplitude = answer;
  SEXP weight = R_NilValue;
  if (TYPEOF(answer) == VECSXP) {
    if (XLENGTH(answer) < 1 || XLENGTH(answer) > 2) {
      ctx->failed = true;
      ctx->error = "response list must contain amplitude and optional weight";
      UNPROTECT(protected_count);
      return MmfirInvalidAmplitude;
    }
    amplitude = VECTOR_ELT(answer, 0);
    if (XLENGTH(answer) == 2) weight = VECTOR_ELT(answer, 1);
  }
  auto numeric_vector = [](SEXP x) {
    return TYPEOF(x) == REALSXP || TYPEOF(x) == INTSXP;
  };
  auto value = [](SEXP x, int i) {
    return TYPEOF(x) == REALSXP ? REAL(x)[i] : double(INTEGER(x)[i]);
  };

  if (!numeric_vector(amplitude) || XLENGTH(amplitude) != len) {
    ctx->failed = true;
    ctx->error = "response amplitude must be a numeric vector of grid length";
    UNPROTECT(protected_count);
    return MmfirInvalidAmplitude;
  }
  if (weight != R_NilValue &&
      (!numeric_vector(weight) || XLENGTH(weight) != len)) {
    ctx->failed = true;
    ctx->error = "response weight must be a numeric vector of grid length";
    UNPROTECT(protected_count);
    return MmfirInvalidWeight;
  }
  for (int i = 0; i < len; ++i) {
    double a = value(amplitude, i);
    double w = weight == R_NilValue ? 1.0 : value(weight, i);
    if (!std::isfinite(a) || !std::isfinite(w) || w <= 0) {
      ctx->failed = true;
      ctx->error = "response values must be finite and weights positive";
      UNPROTECT(protected_count);
      return MmfirInvalidWeight;
    }
    points[i].a = a;
    points[i].w = w;
  }
  UNPROTECT(protected_count);
  return MmfirSuccess;
}

// [[Rcpp::export]]
List firpm_mmfir(
    int n,
    NumericVector f,
    NumericVector a,
    NumericVector w,
    int filterClass,
    double accuracy = 0.0,
    double persistence = 0.0,
    double robustness = 0.0,
    bool diagnostics = false,
    SEXP response = R_NilValue,
    SEXP responseArgs = R_NilValue)
{
  if (n < 1 || n == INT_MAX)
    stop("firpm: n must be a positive filter order");
  if (f.size() < 2 || f.size() % 2 != 0)
    stop("firpm: f must contain pairs of band edges");

  const int numBands = f.size() / 2;
  const int numTaps = n + 1;
  const bool custom = response != R_NilValue;
  if (custom && !Rf_isFunction(response))
    stop("firpm: response must be an R function");
  if (custom && TYPEOF(responseArgs) != VECSXP)
    stop("firpm: responseArgs must be a list");
  if (!custom && a.size() != numBands && a.size() != f.size())
    stop("firpm: invalid amplitude vector length");
  if (w.size() != numBands && w.size() != f.size())
    stop("firpm: invalid weight vector length");
  if (filterClass < 0 || filterClass > 2)
    stop("firpm: invalid filter class");

  std::vector<MmfirBandSpec> bands(numBands);
  for (int i = 0; i < numBands; ++i) {
    MmfirBandSpec b = {};
    b.freqL = f[2 * i];
    b.freqR = f[2 * i + 1];
    if (custom) {
      b.ampL = 0.0;
      b.ampR = 0.0;
    } else if (a.size() == numBands) {
      b.ampL = a[i];
      b.ampR = a[i];
    } else {
      b.ampL = a[2 * i];
      b.ampR = a[2 * i + 1];
    }
    if (w.size() == numBands) {
      b.weightL = w[i];
      b.weightR = w[i];
    } else {
      b.weightL = w[2 * i];
      b.weightR = w[2 * i + 1];
    }
    bands[i] = b;
  }

  NumericVector h(numTaps);
  MmfirStateOwner owner;
  MmfirState *state = diagnostics ? &owner.state : nullptr;

  RResponseContext ctx = {response, responseArgs, n, f, w, false, ""};

  MmfirReport report = mmfir(
      REAL(h),
      static_cast<MmfirFilterClass>(filterClass),
      numTaps,
      numBands,
      bands.data(),
      accuracy,
      persistence,
      robustness,
      0.0,       // target
      0,         // flags
      state,
      custom ? rResponseCallback : nullptr,
      &ctx
  );

  if (ctx.failed)
    stop("firpm: %s", ctx.error.c_str());

  if (report.result >= MmfirError) {
    stop("firpm: %s (result = %d)",
         report.text ? report.text : "unknown mmfir error",
         static_cast<int>(report.result));
  }

  // No coefficient reversal: validated against Octave for all FIR types.
  List ans = List::create(
      Named("b") = h,
      Named("minimax") = report.result == MmfirSuccess
          ? report.minimax : -report.minimax,
      Named("result") = static_cast<int>(report.result),
      Named("message") = report.text ? report.text : "",
      Named("iterations") = report.iterations
  );

  if (diagnostics) {
    const MmfirState &s = owner.state;
    if (!s.space || !s.peaks || s.spaceLength <= 0 || s.peaksLength <= 0)
      stop("firpm: mmfir did not return a complete diagnostic state");

    NumericVector fgrid(s.spaceLength);
    NumericVector desInternal(s.spaceLength);
    NumericVector wtInternal(s.spaceLength);
    NumericVector peaks(s.peaksLength);
    for (int i = 0; i < s.spaceLength; ++i) {
      fgrid[i] = s.space[i].f;
      desInternal[i] = s.space[i].a;
      wtInternal[i] = s.space[i].w;
    }
    for (int i = 0; i < s.peaksLength; ++i)
      peaks[i] = s.peaks[i].f;

    ans["fgrid"] = fgrid;
    ans["des_internal"] = desInternal;
    ans["wt_internal"] = wtInternal;
    ans["peaks"] = peaks;
  }
  return ans;
}
