// CppNumericalSolver
#ifndef WOLFERULE_H_
#define WOLFERULE_H_

#include "../meta.h"

namespace cppoptlib {

/**
 * @brief this tries to guess the correct stepwith in a bisection-style
 * @details WARNING: THIS IS QUITE HACKY. TEST ARMIJO before.
 *
 * @tparam T scalar type
 * @tparam P problem type
 * @tparam Ord order of solver
 */
template<typename T, typename P, int Ord>
class WolfeHeuristic {

 public:

  static T linesearch(const Vector<T> & x0, const Vector<T> & searchDir, P &objFunc, const T alpha_init = 1) {

    Vector<T> x = x0;

    // evaluate phi(0)
    T phi0 = objFunc.value(x0);
    // evaluate phi'(0)
    Vector<T> grad(x.rows());
    objFunc.gradient(x, grad);

    T phi0_dash = searchDir.dot(grad);

    T alpha = alpha_init;
    bool decrease_direction = true;

    // 200 guesses
    for(size_t iter = 0; iter < 200; ++iter) {

      // new guess for phi(alpha)
      Vector<T> x_candidate = x + alpha * searchDir;
      const T phi = objFunc.value(x_candidate);

      // decrease condition invalid --> shrink interval
      if (phi > phi0 + 0.0001 * alpha * phi0_dash) {
        alpha *= 0.5;
        decrease_direction = false;
      } else {

        // valid decrease --> test strong wolfe condition
        Vector<T> grad2(x.rows());
        objFunc.gradient(x_candidate, grad2);
        const T phi_dash = searchDir.dot(grad2);

        // curvature condition invalid ?
        if ((phi_dash < 0.9 * phi0_dash) || !decrease_direction) {
          // increase interval
          alpha *= 4.0;
        } else {
          // both condition are valid --> we are happy
          x = x_candidate;
          grad = grad2;
          phi0 = phi;
          return alpha;
        }
      }
    }

    return alpha;
  }

};
}

#endif /* WOLFERULE_H_ */
