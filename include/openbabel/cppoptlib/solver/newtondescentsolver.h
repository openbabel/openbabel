// CppNumericalSolver
#include <iostream>
#include <Eigen/LU>
#include "isolver.h"
#include "../linesearch/armijo.h"

#ifndef NEWTONDESCENTSOLVER_H_
#define NEWTONDESCENTSOLVER_H_

namespace cppoptlib {

template<typename ProblemType>
class NewtonDescentSolver : public ISolver<ProblemType, 2> {
  public:
    using Superclass = ISolver<ProblemType, 2>;
    using typename Superclass::Scalar;
    using typename Superclass::TVector;
    using typename Superclass::THessian;

    void minimize(ProblemType &objFunc, TVector &x0) {
        const int DIM = x0.rows();
        TVector grad = TVector::Zero(DIM);
        THessian hessian = THessian::Zero(DIM, DIM);
        this->m_current.reset();
        do {
            objFunc.gradient(x0, grad);
            objFunc.hessian(x0, hessian);
            hessian += (1e-5) * THessian::Identity(DIM, DIM);
            TVector delta_x = hessian.lu().solve(-grad);
            const double rate = Armijo<ProblemType, 1>::linesearch(x0, delta_x, objFunc) ;
            x0 = x0 + rate * delta_x;
            // std::cout << "iter: "<<iter<< ", f = " <<  objFunc.value(x0) << ", ||g||_inf "<<gradNorm  << std::endl;
            ++this->m_current.iterations;
            this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
            this->m_status = checkConvergence(this->m_stop, this->m_current);
        } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));
    }
};

}
/* namespace cppoptlib */

#endif /* NEWTONDESCENTSOLVER_H_ */
