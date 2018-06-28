/**********************************************************************************
gaussian_integrals.h - coulombintegrals for Gaussian s-type orbitals

Copyright (C) 2018 by Mohammad M Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
                      David van der Spoel Group
                      Uppsala University

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************************/

#ifndef GAUSSIAN_INTEGRALS_H
#define GAUSSIAN_INTEGRALS_H

#include <cmath>

/*! \brief 
 * Compute the Coulomb overlap integral between two gaussian distributed charges.
 * The widths xi and xj may both be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */  
double Coulomb_GG(double r,double xi,double xj);

/*! \brief 
 * Compute the derivative of the Coulomb overlap integral between two gaussian 
 * distributed charges with respect to the distance.
 * The widths xi and xj may both be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double DCoulomb_GG(double r,double xi,double xj);


/*! \brief 
 * Compute the Coulomb overlap integral between a gaussian distributed charge
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */  
double Nuclear_GG(double r,double xi);

/*! \brief 
 * Compute the derivative of the Coulomb overlap integral between a gaussian 
 * distributed charge and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */   
double DNuclear_GG(double r,double xi);


#endif /*GAUSSIAN_INTEGRALS_H*/
