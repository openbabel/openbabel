// CppNumericalSolver
#ifndef NELDERMEADSOLVER_H_
#define NELDERMEADSOLVER_H_
#include <cmath>
#include <Eigen/Core>
#include "isolver.h"
#include "../meta.h"

namespace cppoptlib {

template<typename ProblemType>
class NelderMeadSolver : public ISolver<ProblemType, 0> {
 public:
  using Superclass = ISolver<ProblemType, 0>;
  using typename Superclass::Scalar;
  using typename Superclass::TVector;
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  MatrixType x0;
  SimplexOp lastOp = SimplexOp::Place;
  Status stop_condition;
  bool initialSimplexCreated = false;

  MatrixType makeInitialSimplex(TVector &x) {
    size_t DIM = x.rows();

    MatrixType s = MatrixType::Zero(DIM, DIM + 1);
    for (int c = 0; c < int(DIM) + 1; ++c) {
      for (int r = 0; r < int(DIM); ++r) {
        s(r, c) = x(r);
        if (r == c - 1) {
          if (x(r) == 0) {
            s(r, c) = 0.00025;
          } else {
            s(r, c) = (1 + 0.05) * x(r);
          }
        }
      }
    }

    return s;
  }

  /**
   * @brief minimize
   * @details [long description]
   *
   * @param objFunc [description]
   */
  void minimize(ProblemType &objFunc, TVector &x) {

    const Scalar rho = 1.;    // rho > 0
    const Scalar xi  = 2.;    // xi  > max(rho, 1)
    const Scalar gam = 0.5;   // 0 < gam < 1

    const size_t DIM = x.rows();

    // create initial simplex
    if (not initialSimplexCreated) {
      x0 = makeInitialSimplex(x);
    }

    // compute function values
    std::vector<Scalar> f; f.resize(DIM + 1);
    std::vector<int> index; index.resize(DIM + 1);
    for (int i = 0; i < int(DIM) + 1; ++i) {
      f[i] = objFunc(static_cast<TVector >(x0.col(i)));
      index[i] = i;
    }

    sort(index.begin(), index.end(), [&](int a, int b)-> bool { return f[a] < f[b]; });

    int iter = 0;
    const int maxIter = this->m_stop.iterations * DIM;
    while (
      objFunc.callback(this->m_current, x0.col(index[0])) and
      (iter < maxIter)
    ) {
      // conv-check
      Scalar max1 = fabs(f[index[1]] - f[index[0]]);
      Scalar max2 = (x0.col(index[1]) - x0.col(index[0]) ).array().abs().maxCoeff();
      for (int i = 2; i < int(DIM) + 1; ++i) {
        Scalar tmp1 = fabs(f[index[i]] - f[index[0]]);
        if (tmp1 > max1)
          max1 = tmp1;

        Scalar tmp2 = (x0.col(index[i]) - x0.col(index[0]) ).array().abs().maxCoeff();
        if (tmp2 > max2)
          max2 = tmp2;
      }
      const Scalar tt1 = std::max<Scalar>(Scalar(1.e-04), 10 * std::nextafter(f[index[0]], std::numeric_limits<Scalar>::epsilon()) - f[index[0]]);
      const Scalar tt2 = std::max<Scalar>(Scalar(1.e-04), 10 * (std::nextafter(x0.col(index[0]).maxCoeff(), std::numeric_limits<Scalar>::epsilon())
                    - x0.col(index[0]).maxCoeff()));

      // User-defined stopping criteria
      this->m_current.iterations = iter;
      this->m_current.fDelta = max1;
      this->m_current.xDelta = max2;
      stop_condition = checkConvergence(this->m_stop, this->m_current);
      if (this->m_stop.iterations != 0 and stop_condition != Status::Continue) {
        break;
      }

      // Allow stopping in the callback. This callback gets the correct current
      // state unlike the simple one in while(), which get previous state.
      if (objFunc.detailed_callback(this->m_current, lastOp, index[0], x0, f) == false) {
        stop_condition = Status::UserDefined;
        break;
      }

      // max(||x - shift(x) ||_inf ) <= tol,
      if (max1 <=  tt1) {
        // values to similar
        if (max2 <= tt2) {
          stop_condition = Status::FDeltaTolerance;
          break;
        }
      }

      //////////////////////////

      // midpoint of the simplex opposite the worst point
      TVector x_bar = TVector::Zero(DIM);
      for (int i = 0; i < int(DIM); ++i) {
        x_bar += x0.col(index[i]);
      }
      x_bar /= Scalar(DIM);

      // Compute the reflection point
      const TVector x_r   = ( 1. + rho ) * x_bar - rho   * x0.col(index[DIM]);
      const Scalar f_r = objFunc(x_r);
      lastOp = SimplexOp::Reflect;

      if (f_r < f[index[0]]) {
        // the expansion point
        const TVector x_e = ( 1. + rho * xi ) * x_bar - rho * xi   * x0.col(index[DIM]);
        const Scalar f_e = objFunc(x_e);
        if ( f_e < f_r ) {
          // expand
          lastOp = SimplexOp::Expand;
          x0.col(index[DIM]) = x_e;
          f[index[DIM]] = f_e;
        } else {
          // reflect
          lastOp = SimplexOp::Reflect;
          x0.col(index[DIM]) = x_r;
          f[index[DIM]] = f_r;
        }
      } else {
        if ( f_r < f[index[DIM - 1]] ) {
          x0.col(index[DIM]) = x_r;
          f[index[DIM]] = f_r;
        } else {
          // contraction
          if (f_r < f[index[DIM]]) {
            const TVector x_c = (1 + rho * gam) * x_bar - rho * gam * x0.col(index[DIM]);
            const Scalar f_c = objFunc(x_c);
            if ( f_c <= f_r ) {
              // outside
              x0.col(index[DIM]) = x_c;
              f[index[DIM]] = f_c;
              lastOp = SimplexOp::ContractOut;
            } else {
              shrink(x0, index, f, objFunc);
              lastOp = SimplexOp::Shrink;
            }
          } else {
            // inside
            const TVector x_c = ( 1 - gam ) * x_bar + gam   * x0.col(index[DIM]);
            const Scalar f_c = objFunc(x_c);
            if (f_c < f[index[DIM]]) {
              x0.col(index[DIM]) = x_c;
              f[index[DIM]] = f_c;
              lastOp = SimplexOp::ContractIn;
            } else {
              shrink(x0, index, f, objFunc);
              lastOp = SimplexOp::Shrink;
            }
          }
        }
      }
      sort(index.begin(), index.end(), [&](int a, int b)-> bool { return f[a] < f[b]; });
      iter++;
      if (iter >= maxIter) {
        stop_condition = Status::IterationLimit;
      }
      else {
        stop_condition = Status::UserDefined; // if stopped in the callback in while()
      }
    } // while loop

    // report the last result
    objFunc.detailed_callback(this->m_current, lastOp, index[0], x0, f);
    x = x0.col(index[0]);
  }

  void shrink(MatrixType &x, std::vector<int> &index, std::vector<Scalar> &f, ProblemType &objFunc) {
    const Scalar sig = 0.5;   // 0 < sig < 1
    const int DIM = x.rows();
    f[index[0]] = objFunc(x.col(index[0]));
    for (int i = 1; i < DIM + 1; ++i) {
      x.col(index[i]) = sig * x.col(index[i]) + (1. - sig) * x.col(index[0]);
      f[index[i]] = objFunc(x.col(index[i]));
    }
  }

  // Need our own checker here to get rid of the gradient test used in other solvers
  template<typename T>
  Status checkConvergence(const Criteria<T> &stop, const Criteria<T> &current) {
    if ((stop.iterations > 0) && (current.iterations > stop.iterations)) {
      return Status::IterationLimit;
    }
    if ((stop.xDelta > 0) && (current.xDelta < stop.xDelta)) {
      return Status::XDeltaTolerance;
    }
    if ((stop.fDelta > 0) && (current.fDelta < stop.fDelta)) {
      return Status::FDeltaTolerance;
    }
    return Status::Continue;
  }

}; /* class NelderMeadSolver */

} /* namespace cppoptlib */

#endif /* NELDERMEADSOLVER_H_ */
