/**********************************************************************************
slater_integrals.h - coulombintegrals for Slater s-type orbitals

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
 
#ifndef SLATER_INTEGRALS_H
#define SLATER_INTEGRALS_H

#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief 
 * Compute the Slater overlap integral between two 1S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */ 
double Slater_1S_1S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between 1S and 2S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_1S_2S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between 1S and 2S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_2S_1S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between 1S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_1S_3S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between 1S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_3S_1S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral two 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_2S_2S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between 2S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_2S_3S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between 2S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_3S_2S(double r, double xi, double xj);

/*! \brief 
 * Compute the Slater overlap integral between two 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double Slater_3S_3S(double r, double xi, double xj);


/*! \brief 
 * Compute the derivative of Slater overlap integral between two 1S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_1S_1S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between 1S and 2S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_1S_2S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between 1S and 2S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_2S_1S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between 1S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_1S_3S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between 1S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_3S_1S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between two 2S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_2S_2S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between 2S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_2S_3S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between 2S and 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_3S_2S(double r, double xi, double xj);

/*! \brief 
 * Compute the derivative of Slater overlap integral between two 3S Slater Orbitals.
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */
double DSlater_3S_3S(double r, double xi, double xj);

/*! \brief 
 * Compute the Coulomb overlap integral between 1S Slater orbital
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */
double Nuclear_1S(double r, double xi);

/*! \brief 
 * Compute the Coulomb overlap integral between 2S Slater orbital
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */
double Nuclear_2S(double r, double xi);

/*! \brief 
 * Compute the Coulomb overlap integral between 3S Slater orbital
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */
double Nuclear_3S(double r, double xi);

/*! \brief 
 * Compute the derivative of Coulomb overlap integral between 1S Slater orbital
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */
double DNuclear_1S(double r, double xi);

/*! \brief 
 * Compute the derivative of Coulomb overlap integral between 2S Slater orbital
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */
double DNuclear_2S(double r, double xi);

/*! \brief 
 * Compute the derivative of Coulomb overlap integral between 2S Slater orbital
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */
double DNuclear_3S(double r, double xi);

/*! \brief 
 * Compute the Slater overlap integral between two Slater distributed charges.
 * The widths xi and xj may both be zero. The Slater function types i and j
 * should be in the range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] j  Slater function type
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */    
double Coulomb_SS(double r,int i,int j,double xi,double xj);

/*! \brief 
 * Compute the derivative of the Slater overlap integral between two Slater 
 * distributed charges with respect to r.
 * The widths xi and xj may both be zero. The Slater function types i and j
 * should be in the range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] j  Slater function type
 * \param[in] xi orbital exponent in 1/nm 
 * \param[in] xj orbital exponent in 1/nm 
 * \return    Integral value
 */    
double DCoulomb_SS(double r,int i,int j,double xi,double xj);


/*! \brief 
 * Compute the Slater overlap integral between a Slater distributed charge and a
 * point charge. The width xi may be zero. The Slater function type i should be in the
 * range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */   
double Nuclear_SS(double r,int i, double xi);

/*! \brief 
 * Compute the derivative of the Slater overlap integral between a Slater 
 * distributed charge and a point charge with respect to r.
 * The width xi may be zero. The Slater function type i should be in the
 * range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] xi orbital exponent in 1/nm 
 * \return    Integral value
 */  
double DNuclear_SS(double r,int i, double xi);

#ifdef __cplusplus
}
#endif /*SLATER_INTEGRALS_H*/

#endif
