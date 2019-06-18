// CppNumericalSolver
#ifndef CONJUGATEDGRADIENTDESCENTSOLVER_H_
#define CONJUGATEDGRADIENTDESCENTSOLVER_H_

#include <Eigen/Core>
#include "isolver.h"
#include "../linesearch/armijo.h"

namespace cppoptlib {

template<typename ProblemType>
class ConjugatedGradientDescentSolver : public ISolver<ProblemType, 1> {

 public:
  using Superclass = ISolver<ProblemType, 1>;
  using typename Superclass::Scalar;
  using typename Superclass::TVector;

  /**
   * @brief minimize
   * @details [long description]
   *
   * @param objFunc [description]
   */
  void minimize(ProblemType &objFunc, TVector &x0) {
    TVector grad(x0.rows());
    TVector grad_old(x0.rows());
    TVector Si(x0.rows());
    TVector Si_old(x0.rows());

    this->m_current.reset();
    do {
      objFunc.gradient(x0, grad);

      if (this->m_current.iterations == 0) {
        Si = -grad;
      } else {
        const double beta = grad.dot(grad) / (grad_old.dot(grad_old));
        Si = -grad + beta * Si_old;
      }

      const double rate = Armijo<ProblemType, 1>::linesearch(x0, Si, objFunc) ;

      x0 = x0 + rate * Si;

      grad_old = grad;
      Si_old = Si;

      this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
      // std::cout << "iter: "<<iter<< " f = " <<  objFunc.value(x0) << " ||g||_inf "<<gradNorm   << std::endl;
      ++this->m_current.iterations;
      this->m_status = checkConvergence(this->m_stop, this->m_current);
    } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue) );

  }

};

} /* namespace cppoptlib */

#endif /* CONJUGATEDGRADIENTDESCENTSOLVER_H_ */
