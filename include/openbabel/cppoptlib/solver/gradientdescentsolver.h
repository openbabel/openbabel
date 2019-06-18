// CppNumericalSolver
#ifndef GRADIENTDESCENTSOLVER_H_
#define GRADIENTDESCENTSOLVER_H_

#include <Eigen/Core>
#include "isolver.h"
#include "../linesearch/morethuente.h"

namespace cppoptlib {

template<typename ProblemType>
class GradientDescentSolver : public ISolver<ProblemType, 1> {

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

    TVector direction(x0.rows());
    this->m_current.reset();
    do {
      ;
      objFunc.gradient(x0, direction);
      const Scalar rate = MoreThuente<ProblemType, 1>::linesearch(x0, -direction, objFunc) ;
      x0 = x0 - rate * direction;
      this->m_current.gradNorm = direction.template lpNorm<Eigen::Infinity>();
      // std::cout << "iter: "<<iter<< " f = " <<  objFunc.value(x0) << " ||g||_inf "<<gradNorm  << std::endl;
      ++this->m_current.iterations;
      this->m_status = checkConvergence(this->m_stop, this->m_current);
    } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));
    if (this->m_debug > DebugLevel::None) {
        std::cout << "Stop status was: " << this->m_status << std::endl;
        std::cout << "Stop criteria were: " << std::endl << this->m_stop << std::endl;
        std::cout << "Current values are: " << std::endl << this->m_current << std::endl;
    }
  }

};

} /* namespace cppoptlib */

#endif /* GRADIENTDESCENTSOLVER_H_ */
