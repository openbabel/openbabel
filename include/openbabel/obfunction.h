/**********************************************************************
obfunction.h - Base class for force fields and scroring functions which
               depend on atom positions.
 
Copyright (C) 2008 by Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OPENBABEL_OBFUNCTION_H
#define OPENBABEL_OBFUNCTION_H

#include <vector>
#include <Eigen/Core>

namespace OpenBabel {

  /** @class OBFunction
   *  @brief Base class for force fields and scoring functions which 
   *  depend on atom coordinates. 
   */
  class OBFunction {
    protected:
      std::vector<Eigen::Vector3d> m_positions;
      std::vector<Eigen::Vector3d> m_gradients;
    public:
      /**
       * Constructor
       */
      OBFunction() {}
      /**
       * Destructor
       */
      virtual ~OBFunction() {}
      /**
       * Evaluate the function. This method is implemented by subclasses.
       */
      virtual double Eval(bool gradients = true) = 0;
      /**
       * @return True if this function has analytical gradients. 
       */
      virtual bool HasAnalyticalGradients() const { return false; }
      /** 
       * @brief Get the atom positions.
       */
      std::vector<Eigen::Vector3d>&  GetPositions() { return m_positions; }
      /** 
       * @brief Get the atom gradients.
       */
      std::vector<Eigen::Vector3d>&  GetGradients() { return m_gradients; } 
      /** 
       * @brief Calculate the potential energy function derivative numerically with 
       * repect to the coordinates of atom with index @p idx (this vector is the gradient)
       *
       * @code
       *         f(1) - f(0)
       * f'(0) = -----------      f(1) = f(0+h)
       *              h
       * @endcode
       *
       * @param idx Index for the atom starting at 0.
       * 
       * @return The gradient for atom with index @p idx.
       */
      Eigen::Vector3d NumericalDerivative(unsigned int idx);
      /** 
       * @brief Calculate the potential energy function 2nd derivative numerically with 
       * repect to the coordinates of atom with index @p idx.
       *
       * @code
       *         f(2) - 2f(1) + f(0)
       * f'(0) = -------------------      f(1) = f(0+h)
       *                 h^2              f(1) = f(0+2h)
       * @endcode
       *
       * @param idx Index for the atom starting at 0.
       * 
       * @return The 2nd derivative for atom with index @p idx.
       */ 
      Eigen::Vector3d NumericalSecondDerivative(unsigned int idx);
 
  };

} // namespace OpenBabel
#endif
