#ifndef PROBLEM_H
#define PROBLEM_H

#include <array>
#include <vector>
#include <Eigen/Core>

#include "meta.h"

namespace cppoptlib {

template<typename Scalar_, int Dim_ = Eigen::Dynamic>
class Problem {
 public:
  static const int Dim = Dim_;
  typedef Scalar_ Scalar;
  using TVector   = Eigen::Matrix<Scalar, Dim, 1>;
  using THessian  = Eigen::Matrix<Scalar, Dim, Dim>;
  using TCriteria = Criteria<Scalar>;
  using TIndex = typename TVector::Index;
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

 public:
  Problem() {}
  virtual ~Problem()= default;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
  virtual bool callback(const Criteria<Scalar> &state, const TVector &x) {
    return true;
  }

  virtual bool detailed_callback(const Criteria<Scalar> &state, SimplexOp op, int index, const MatrixType &x, std::vector<Scalar> f) {
    return true;
  }
#pragma GCC diagnostic pop

  /**
   * @brief returns objective value in x
   * @details [long description]
   *
   * @param x [description]
   * @return [description]
   */
  virtual Scalar value(const  TVector &x) = 0;
  /**
   * @brief overload value for nice syntax
   * @details [long description]
   *
   * @param x [description]
   * @return [description]
   */
  Scalar operator()(const  TVector &x) {
    return value(x);
  }
  /**
   * @brief returns gradient in x as reference parameter
   * @details should be overwritten by symbolic gradient
   *
   * @param grad [description]
   */
  virtual void gradient(const  TVector &x,  TVector &grad) {
    finiteGradient(x, grad);
  }

  /**
   * @brief This computes the hessian
   * @details should be overwritten by symbolic hessian, if solver relies on hessian
   */
  virtual void hessian(const TVector &x, THessian &hessian) {
    finiteHessian(x, hessian);
  }

  virtual bool checkGradient(const TVector &x, int accuracy = 3) {
    // TODO: check if derived class exists:
    // int(typeid(&Rosenbrock<double>::gradient) == typeid(&Problem<double>::gradient)) == 1 --> overwritten
    const TIndex D = x.rows();
    TVector actual_grad(D);
    TVector expected_grad(D);
    gradient(x, actual_grad);
    finiteGradient(x, expected_grad, accuracy);
    for (TIndex d = 0; d < D; ++d) {
      Scalar scale = std::max(static_cast<Scalar>(std::max(fabs(actual_grad[d]), fabs(expected_grad[d]))), Scalar(1.));
      if(fabs(actual_grad[d]-expected_grad[d])>1e-2 * scale)
        return false;
    }
    return true;

  }

  virtual bool checkHessian(const TVector &x, int accuracy = 3) {
    // TODO: check if derived class exists:
    // int(typeid(&Rosenbrock<double>::gradient) == typeid(&Problem<double>::gradient)) == 1 --> overwritten
    const TIndex D = x.rows();

    THessian actual_hessian = THessian::Zero(D, D);
    THessian expected_hessian = THessian::Zero(D, D);
    hessian(x, actual_hessian);
    finiteHessian(x, expected_hessian, accuracy);
    for (TIndex d = 0; d < D; ++d) {
      for (TIndex e = 0; e < D; ++e) {
        Scalar scale = std::max(static_cast<Scalar>(std::max(fabs(actual_hessian(d, e)), fabs(expected_hessian(d, e)))), Scalar(1.));
        if(fabs(actual_hessian(d, e)- expected_hessian(d, e))>1e-1 * scale)
          return false;
      }
    }
    return true;
  }

  void finiteGradient(const  TVector &x, TVector &grad, int accuracy = 0) {
    // accuracy can be 0, 1, 2, 3
    const Scalar eps = 2.2204e-6;
    static const std::array<std::vector<Scalar>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    static const std::array<std::vector<Scalar>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    static const std::array<Scalar, 4> dd = {2, 12, 60, 840};

    grad.resize(x.rows());
    TVector& xx = const_cast<TVector&>(x);

    const int innerSteps = 2*(accuracy+1);
    const Scalar ddVal = dd[accuracy]*eps;

    for (TIndex d = 0; d < x.rows(); d++) {
      grad[d] = 0;
      for (int s = 0; s < innerSteps; ++s)
      {
        Scalar tmp = xx[d];
        xx[d] += coeff2[accuracy][s]*eps;
        grad[d] += coeff[accuracy][s]*value(xx);
        xx[d] = tmp;
      }
      grad[d] /= ddVal;
    }
  }

  void finiteHessian(const TVector &x, THessian &hessian, int accuracy = 0) {
    const Scalar eps = std::numeric_limits<Scalar>::epsilon()*10e7;

    hessian.resize(x.rows(), x.rows());
    TVector& xx = const_cast<TVector&>(x);

    if(accuracy == 0) {
      for (TIndex i = 0; i < x.rows(); i++) {
        for (TIndex j = 0; j < x.rows(); j++) {
          Scalar tmpi = xx[i];
          Scalar tmpj = xx[j];

          Scalar f4 = value(xx);
          xx[i] += eps;
          xx[j] += eps;
          Scalar f1 = value(xx);
          xx[j] -= eps;
          Scalar f2 = value(xx);
          xx[j] += eps;
          xx[i] -= eps;
          Scalar f3 = value(xx);
          hessian(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);

          xx[i] = tmpi;
          xx[j] = tmpj;
        }
      }
    } else {
      /*
        \displaystyle{{\frac{\partial^2{f}}{\partial{x}\partial{y}}}\approx
        \frac{1}{600\,h^2} \left[\begin{matrix}
          -63(f_{1,-2}+f_{2,-1}+f_{-2,1}+f_{-1,2})+\\
          63(f_{-1,-2}+f_{-2,-1}+f_{1,2}+f_{2,1})+\\
          44(f_{2,-2}+f_{-2,2}-f_{-2,-2}-f_{2,2})+\\
          74(f_{-1,-1}+f_{1,1}-f_{1,-1}-f_{-1,1})
        \end{matrix}\right] }
      */
      for (TIndex i = 0; i < x.rows(); i++) {
        for (TIndex j = 0; j < x.rows(); j++) {
          Scalar tmpi = xx[i];
          Scalar tmpj = xx[j];

          Scalar term_1 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += -2*eps;  term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += -1*eps;  term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += 1*eps;   term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += 2*eps;   term_1 += value(xx);

          Scalar term_2 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += -2*eps;  term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += -1*eps;  term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += 2*eps;   term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += 1*eps;   term_2 += value(xx);

          Scalar term_3 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += -2*eps;  term_3 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += 2*eps;   term_3 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += -2*eps;  term_3 -= value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += 2*eps;   term_3 -= value(xx);

          Scalar term_4 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += -1*eps;  term_4 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += 1*eps;   term_4 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += -1*eps;  term_4 -= value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += 1*eps;   term_4 -= value(xx);

          xx[i] = tmpi;
          xx[j] = tmpj;

          hessian(i, j) = (-63 * term_1+63 * term_2+44 * term_3+74 * term_4)/(600.0 * eps * eps);
        }
      }
    }

  }

};
}

#endif /* PROBLEM_H */
