// ultrwin.cc
// Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
// Original code for Octave:
// Copyright (C) 2013 Rob Sykes <robs@users.sourceforge.net>
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
// 20201116  GvB       setup for gsignal v0.1.0
//---------------------------------------------------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

#if !defined M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined (__cplusplus) && __cplusplus > 201402L
#  define ATTR_FALLTHROUGH [[fallthrough]]
#elif defined (__GNUC__) && __GNUC__ >= 7
#  define ATTR_FALLTHROUGH __attribute__ ((__fallthrough__))
#else
#  define ATTR_FALLTHROUGH ((void) 0)
#endif

#if DEBUG_ULTRWIN
#define SHOW1(x) fprintf(stderr, "%c%+.3e", " 1"[(x) > .5], (x) - ((x) > .5))
#endif

#define EPSILON     (1./0x100000) /* For caller rounding error. */
#define BETA_MAX    (12*(1+EPSILON))
#define XMU_MIN     (.99*(1-EPSILON))


static double *
  ultraspherical_win(int n, double mu, double xmu)
  {
    double * w = 0;
    int bad = n < 1 || xmu < XMU_MIN || (!mu && xmu == 1) ||
      (n > BETA_MAX * 2 && xmu * cos(M_PI * BETA_MAX / n) > 1);
    if (!bad && (w = (double *)malloc(sizeof(*w)* n))) {
      int i, j, k, l = 0, m = (n + 1) / 2, met;
      double * divs = w + m - 1, c = 1 - 1 / (xmu * xmu), t, u, v[64], vp, s;
      for (i = 0; i < (int)(sizeof(v) / sizeof(v[0])); v[i++] = 0);
      if (n > 1) for (i = 0; i < m; l = j - (j <= i++)) {
        vp = *v, s = *v = i? (*v + v[1]) * mu * (divs[i] = 1./i) : 1;
        for (met = 0, j = 1, u = 1; ; ++l, v[l] = vp * (i - l) / (mu + l - 1)) {
#define _ t = v[j], v[j] += vp, vp = t, t = s, s += \
          v[j] * (u *= c * (n - i - j) * divs[j]), met = s && s == t, ++j,
        for (k = ((l-j+1) & ~7) + j; j < k && !met; _ _ _ _ _ _ _ _ (void)0) (void)0;
          for (; j <= l && !met; _ (void)0);
#undef _
          if (met || !(j <= i)) break;
        }
        w[i] = s / (n - i - 1);
      }
      else w[0] = 1;
      u = 1 / w[i = m - 1], w[i] = 1;
      for (--i ; i >= 0; u *= (n - 2 - i + mu) / (n - 2 - i), w[i] *= u, --i);
      for (i = 0; i < m; w[n - 1 - i] = w[i], ++i);
    }
    return w;
  }


typedef struct {double f, fp;} uspv_t;


static uspv_t
  ultraspherical_polyval(int n, double mu, double x, double const *divs)
  {
    double fp = n > 0? 2 * x * mu : 1, fpp = 1, f;
    uspv_t result;
    int i, k;
#define _ f = (2*x*(i+mu)*fp - (i+2*mu-1)*fpp) * divs[i+1], fpp=fp, fp=f, ++i,
    for (i = 1, k = i + ((n - i) & ~7); i < k; _ _ _ _ _ _ _ _ (void)0);
    for (; i < n; _ (void)0);
#undef _
    result.f = fp, result.fp = fpp;
    return result;
  }


#define MU_EPSILON  (1./0x4000)
#define EQ(mu,x)    (fabs((mu)-(x)) < MU_EPSILON)


static uspv_t
  ultraspherical_polyval2(     /* With non-+ve integer protection */
int n, double mu, double x, double const * divs)
  {
    int sign = (~(int)floor(mu) & ~(~2u/2))? 1:-1; /* -ve if floor(mu) <0 & odd */
    uspv_t r;
    if (mu < MU_EPSILON && EQ(mu,(int)mu))
      mu = floor(mu + .5) + MU_EPSILON * ((int)mu > mu? -1:1);
    r = ultraspherical_polyval(n, mu, x, divs);
    r.f *= sign, r.fp *= sign;
    return r;
  }


static double
  find_zero(int n, double mu, int l, double extremum_mag, double ripple_ratio,
            double lower_bound, double const *divs)
  {
    double dx, x0, t, x, epsilon = 1e-10, one_over_deriv, target = 0;
    int i, met = 0;
    if (!divs)
      return 0;
    if (!l) {
      double r = ripple_ratio;   /* FIXME: factor in weighted extremum_mag here */
      x = r > 1 ? cosh(acosh(r) / n) : cos(acos(r) / n); /* invert chebpoly-1st */
      x0 = x *= lower_bound / cos(M_PI * .5 / n) + epsilon;
      target = log(extremum_mag * ripple_ratio);
    }
    else {
      double cheb1 = cos(M_PI * (l - .5) / n), cheb2 = cos(M_PI * l / (n + 1));
      if (mu < 1 - l && EQ((int)(mu+.5),mu+.5)) x = met = 1;
      else if (EQ(mu,0)) x = cheb1, met = 1;               /* chebpoly-1st-kind */
      else if (EQ(mu,1)) x = cheb2, met = 1;               /* chebpoly-2nd-kind */
      else x = (cheb1 * cheb2) / (mu * cheb1 + (1 - mu) * cheb2);
      x0 = x;
    }
    for (i = 0; i < 24 && !met; ++i, met = fabs(dx) < epsilon) {/*Newton-Raphson*/
      uspv_t r = ultraspherical_polyval2(n, mu, x, divs);
      if (!(t = ((2*mu + n-1) * r.fp - n*x * r.f)))         /* Fail if div by 0 */
        break;
      one_over_deriv = (1 - x*x) / t;    /* N-R slow for deriv~=1, so take log: */
      if (!l) {                               /* d/dx(f(g(x))) = f'(g(x)).g'(x) */
        one_over_deriv *= r.f;                             /* d/dx(log x) = 1/x */
        if (r.f <= 0)                                 /* Fail if log of non-+ve */
          break;
        if (x + (dx = (target - log(r.f)) * one_over_deriv) <= lower_bound)
          dx = (lower_bound - x) * .875;
        x += dx;
      }
      else x += dx = -r.f * one_over_deriv;
#if DEBUG_ULTRWIN
      fprintf(stderr, "1/deriv=%9.2e dx=%9.2e x=", one_over_deriv, dx);
      SHOW1(x); fprintf(stderr, "\n");
#endif
    }
#if DEBUG_ULTRWIN
    fprintf(stderr, "find_zero(n=%i mu=%g l=%i target=%g r=%g x0=",
            n, mu, l, target, ripple_ratio);
    SHOW1(x0); fprintf(stderr, ") %s ", met? "converged to" : "FAILED at");
    SHOW1(x); fprintf(stderr, " in %i iterations\n", i);
#else
    static_cast<void>(x0);
#endif
    return met? x : 0;
  }


static double *
  make_divs(int n, double **divs)
  {
    int i;
    if (!*divs) {
      *divs = (double *)malloc(n * sizeof(**divs));
      if (*divs)
        for (i = 0; i < n; (*divs)[i] = 1./(i+1), ++i);
    }
    return *divs? *divs - 1 : 0;
  }


#define DIVS make_divs(n, &divs)


typedef enum {uswpt_Xmu, uswpt_Beta, uswpt_AttFirst, uswpt_AttLast} uswpt_t;


double *
  ultraspherical_window(int n, double mu, double par, uswpt_t type, int even_norm,
                        double *xmu_)
  {
    double * w = 0, xmu = 0, * divs = 0, last_extremum_pos = 0;

    if (n > 0 && fabs(mu) <= (8*(1+EPSILON))) switch (type) {
    case uswpt_Beta:
      xmu = mu == 1 && par == 1? 1 : par < .5 || par > BETA_MAX? 0 :
      find_zero(n-1, mu, 1, 0, 0, 0, DIVS) / cos(M_PI * par / n);
      break;

    case uswpt_AttFirst:
      if (par < 0) break;
      ATTR_FALLTHROUGH;

    case uswpt_AttLast:
      if (type == uswpt_AttLast && mu >= 0 && par < 0);
      else if (!EQ(mu,0)) {
        int extremum_num =
          type == uswpt_AttLast? (int)((n-2)/2 +.5) : 1 + EQ(mu,-1.5);
        double extremum_pos =
          find_zero(n-2, mu+1, extremum_num, 0, 0, 0, DIVS);
        double extremum_mag = !extremum_pos? 0 :
          fabs(ultraspherical_polyval2(n-1, mu, extremum_pos, DIVS).f);
        double xmu_lower_bound = !extremum_mag? 0 :
          find_zero(n-1, mu, 1, 0, 0, 0, DIVS); /* 1st null */
        xmu = !xmu_lower_bound? 0 : find_zero(
          n-1, mu, 0, extremum_mag, pow(10, par/20), xmu_lower_bound, DIVS);
        last_extremum_pos =
          type == uswpt_AttLast? extremum_pos : last_extremum_pos;
      }
      else xmu = cosh(acosh(pow(10, par/20))/(n-1)); /* Cheby 1st kind */
      break;

    default: case uswpt_Xmu: xmu = par; break;
    }
#if DEBUG_ULTRWIN
    fprintf(stderr, "n=%i mu=%.3f xmu=%.16g\n", n, mu, xmu);
#endif

    if (xmu > 0)
      w = ultraspherical_win(n, mu, xmu);

    if (w && (~n & !!even_norm) && n > 2 && !(mu == 1 && xmu == 1)) {
      int i = n / 2 - 1, j = 1;
      double * d = DIVS, t = 0, s = -1, x = even_norm == 1? 0 : last_extremum_pos?
      last_extremum_pos : find_zero(n-2, mu+1, i, 0, 0, 0, d);
      x = x? M_PI/2 - acos(x/xmu) : 0;
      for (; i >= 0; t += w[i] * d[j] * (s=-s) * (x?cos(j*x):1), --i, j += 2);
      for (t = M_PI/4 / t, i = 0; t < 1 && i < n; w[i] *= t, ++i);
#if DEBUG_ULTRWIN
      fprintf(stderr, "%snorm DFT(w.sinc πx) @ %g %.16g\n", t<1? "":"NO ", 2*x,t);
#endif
    }
    free(divs);
    if (xmu_)
      *xmu_ = xmu;
    return w;
  }


//----------------------------------------------------------------------------

// [[Rcpp::export]]
Nullable<NumericVector> ultrwin (int m, double mu, double par, int par_type, int even_norm)
{
    double xmu;
    double *w = ultraspherical_window (m, mu, par, static_cast<uswpt_t>(par_type), even_norm, &xmu);

    if (!w)
    {
      return R_NilValue;
    }

    NumericVector ww (m);
    for (int i = 0; i < m; i++)
      ww(i) = w[i];

    return ww;
}
