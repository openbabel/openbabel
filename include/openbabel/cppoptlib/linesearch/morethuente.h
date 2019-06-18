// CppNumericalSolver
#ifndef MORETHUENTE_H_
#define MORETHUENTE_H_

#include "../meta.h"
#include <cmath>

namespace cppoptlib {

template<typename ProblemType, int Ord>
class MoreThuente {

 public:
  using Scalar = typename ProblemType::Scalar;
  using TVector = typename ProblemType::TVector;

  /**
   * @brief use MoreThuente Rule for (strong) Wolfe conditiions
   * @details [long description]
   *
   * @param searchDir search direction for next update step
   * @param objFunc handle to problem
   *
   * @return step-width
   */

  static Scalar linesearch(const TVector &x, const TVector &searchDir, ProblemType &objFunc, const  Scalar alpha_init = 1.0) {
    // assume step width
    Scalar ak = alpha_init;

    Scalar fval = objFunc.value(x);
    TVector  g  = x.eval();
    objFunc.gradient(x, g);

    TVector s = searchDir.eval();
    TVector xx = x.eval();

    cvsrch(objFunc, xx, fval, g, ak, s);

    return ak;
  }

  static int cvsrch(ProblemType &objFunc, TVector &x, Scalar f, TVector &g, Scalar &stp, TVector &s) {
    // we rewrite this from MIN-LAPACK and some MATLAB code
    int info           = 0;
    int infoc          = 1;
    const Scalar xtol   = 1e-15;
    const Scalar ftol   = 1e-4;
    const Scalar gtol   = 1e-2;
    const Scalar stpmin = 1e-15;
    const Scalar stpmax = 1e15;
    const Scalar xtrapf = 4;
    const int maxfev   = 20;
    int nfev           = 0;

    Scalar dginit = g.dot(s);
    if (dginit >= 0.0) {
      // no descent direction
      // TODO: handle this case
      return -1;
    }

    bool brackt      = false;
    bool stage1      = true;

    Scalar finit      = f;
    Scalar dgtest     = ftol * dginit;
    Scalar width      = stpmax - stpmin;
    Scalar width1     = 2 * width;
    TVector wa = x.eval();

    Scalar stx        = 0.0;
    Scalar fx         = finit;
    Scalar dgx        = dginit;
    Scalar sty        = 0.0;
    Scalar fy         = finit;
    Scalar dgy        = dginit;

    Scalar stmin;
    Scalar stmax;

    while (true) {

      // make sure we stay in the interval when setting min/max-step-width
      if (brackt) {
        stmin = std::min<Scalar>(stx, sty);
        stmax = std::max<Scalar>(stx, sty);
      } else {
        stmin = stx;
        stmax = stp + xtrapf * (stp - stx);
      }

      // Force the step to be within the bounds stpmax and stpmin.
      stp = std::max<Scalar>(stp, stpmin);
      stp = std::min<Scalar>(stp, stpmax);

      // Oops, let us return the last reliable values
      if (
      (brackt && ((stp <= stmin) || (stp >= stmax)))
      || (nfev >= maxfev - 1 ) || (infoc == 0)
      || (brackt && ((stmax - stmin) <= (xtol * stmax)))) {
        stp = stx;
      }

      // test new point
      x = wa + stp * s;
      f = objFunc.value(x);
      objFunc.gradient(x, g);
      nfev++;
      Scalar dg = g.dot(s);
      Scalar ftest1 = finit + stp * dgtest;

      // all possible convergence tests
      if ((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0))
        info = 6;

      if ((stp == stpmax) & (f <= ftest1) & (dg <= dgtest))
        info = 5;

      if ((stp == stpmin) & ((f > ftest1) | (dg >= dgtest)))
        info = 4;

      if (nfev >= maxfev)
        info = 3;

      if (brackt & (stmax - stmin <= xtol * stmax))
        info = 2;

      if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit)))
        info = 1;

      // terminate when convergence reached
      if (info != 0)
        return -1;

      if (stage1 & (f <= ftest1) & (dg >= std::min<Scalar>(ftol, gtol)*dginit))
        stage1 = false;

      if (stage1 & (f <= fx) & (f > ftest1)) {
        Scalar fm = f - stp * dgtest;
        Scalar fxm = fx - stx * dgtest;
        Scalar fym = fy - sty * dgtest;
        Scalar dgm = dg - dgtest;
        Scalar dgxm = dgx - dgtest;
        Scalar dgym = dgy - dgtest;

        cstep( stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);

        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;
      } else {
        // this is ugly and some variables should be moved to the class scope
        cstep( stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
      }

      if (brackt) {
        if (fabs(sty - stx) >= 0.66 * width1)
          stp = stx + 0.5 * (sty - stx);
        width1 = width;
        width = fabs(sty - stx);
      }
    }

    return 0;
  }

  static int cstep(Scalar& stx, Scalar& fx, Scalar& dx, Scalar& sty, Scalar& fy, Scalar& dy, Scalar& stp,
  Scalar& fp, Scalar& dp, bool& brackt, Scalar& stpmin, Scalar& stpmax, int& info) {
    info = 0;
    bool bound = false;

    // Check the input parameters for errors.
    if ((brackt & ((stp <= std::min<Scalar>(stx, sty) ) | (stp >= std::max<Scalar>(stx, sty)))) | (dx * (stp - stx) >= 0.0)
    | (stpmax < stpmin)) {
      return -1;
    }

    Scalar sgnd = dp * (dx / fabs(dx));

    Scalar stpf = 0;
    Scalar stpc = 0;
    Scalar stpq = 0;

    if (fp > fx) {
      info = 1;
      bound = true;
      Scalar theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
      Scalar s = std::max<Scalar>(theta, std::max<Scalar>(dx, dp));
      Scalar gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
      if (stp < stx)
        gamma = -gamma;
      Scalar p = (gamma - dx) + theta;
      Scalar q = ((gamma - dx) + gamma) + dp;
      Scalar r = p / q;
      stpc = stx + r * (stp - stx);
      stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
      if (fabs(stpc - stx) < fabs(stpq - stx))
        stpf = stpc;
      else
        stpf = stpc + (stpq - stpc) / 2;
      brackt = true;
    } else if (sgnd < 0.0) {
      info = 2;
      bound = false;
      Scalar theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      Scalar s = std::max<Scalar>(theta, std::max<Scalar>(dx, dp));
      Scalar gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));
      if (stp > stx)
        gamma = -gamma;

      Scalar p = (gamma - dp) + theta;
      Scalar q = ((gamma - dp) + gamma) + dx;
      Scalar r = p / q;
      stpc = stp + r * (stx - stp);
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (fabs(stpc - stp) > fabs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      brackt = true;
    } else if (fabs(dp) < fabs(dx)) {
      info = 3;
      bound = 1;
      Scalar theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      Scalar s = std::max<Scalar>(theta, std::max<Scalar>( dx, dp));
      Scalar gamma = s * sqrt(std::max<Scalar>(static_cast<Scalar>(0.), (theta / s) * (theta / s) - (dx / s) * (dp / s)));
      if (stp > stx)
        gamma = -gamma;
      Scalar p = (gamma - dp) + theta;
      Scalar q = (gamma + (dx - dp)) + gamma;
      Scalar r = p / q;
      if ((r < 0.0) & (gamma != 0.0)) {
        stpc = stp + r * (stx - stp);
      } else if (stp > stx) {
        stpc = stpmax;
      } else {
        stpc = stpmin;
      }
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (brackt) {
        if (fabs(stp - stpc) < fabs(stp - stpq)) {
          stpf = stpc;
        } else {
          stpf = stpq;
        }
      } else {
        if (fabs(stp - stpc) > fabs(stp - stpq)) {
          stpf = stpc;
        } else {
          stpf = stpq;
        }

      }
    } else {
      info = 4;
      bound = false;
      if (brackt) {
        Scalar theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
        Scalar s = std::max<Scalar>(theta, std::max<Scalar>(dy, dp));
        Scalar gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
        if (stp > sty)
          gamma = -gamma;

        Scalar p = (gamma - dp) + theta;
        Scalar q = ((gamma - dp) + gamma) + dy;
        Scalar r = p / q;
        stpc = stp + r * (sty - stp);
        stpf = stpc;
      } else if (stp > stx)
        stpf = stpmax;
      else {
        stpf = stpmin;
      }
    }

    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    } else {
      if (sgnd < 0.0) {
        sty = stx;
        fy = fx;
        dy = dx;
      }

      stx = stp;
      fx = fp;
      dx = dp;
    }

    stpf = std::min<Scalar>(stpmax, stpf);
    stpf = std::max<Scalar>(stpmin, stpf);
    stp = stpf;

    if (brackt & bound) {
      if (sty > stx) {
        stp = std::min<Scalar>(stx + static_cast<Scalar>(0.66) * (sty - stx), stp);
      } else {
        stp = std::max<Scalar>(stx + static_cast<Scalar>(0.66) * (sty - stx), stp);
      }
    }

    return 0;

  }

};

}

#endif /* MORETHUENTE_H_ */
