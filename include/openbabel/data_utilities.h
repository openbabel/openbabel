/**********************************************************************
data_utilities.h - Global data and resource file parsers.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Copyright (C) 2015 by David van der Spoel

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_DATA_UTILITIES_H
#define OB_DATA_UTILITIES_H

#include <openbabel/babelconfig.h>

#include <vector>

namespace OpenBabel {
 class OBMol;

 const double HARTEE_TO_KCALPERMOL = 627.509469;
 const double HARTREE_TO_KJPERMOL = 2625.49962;
 const double KJPERMOL_TO_KCALPERMOL = 1.0/4.184;
 const double RYDBERG_TO_KCALPERMOL = 313.755026;
 const double ELECTRONVOLT_TO_KCALPERMOL = 23.060538;
  
/*! \brief
 * Convenience function to extract thermochemistry from a molecule structure
 *
 * \param[in] mol          The molecule structure
 * \param[in] bVerbose     If true will print information 
 * \param[inout] Nsymm     If not zero and differing from the rotational symmetry
 *                         in the input molecule, corrections to the entropy and
 *                         free energy will be applied. If zero will hold the symmetry
 *                         number from the input molecule on return.
 * \param[out] temperature The temperature
 * \param[out] DeltaHf0    Enthalpy of formation at T = 0
 * \param[out] DeltaHfT    Enthalpy of formation at T
 * \param[out] DeltaGfT    Gibbs energy of formation at T
 * \param[out] DeltaSfT    Entropy of formation at T
 * \param[out] S0T         Standard entropy at T
 * \param[out] CVT         Heat capacity at T and constant Volume
 * \param[out] Scomponents Translational, Rotational and Vibrational components of S0
 * \return true if all values were found, false otherwise.
 */
 OBAPI bool extract_thermochemistry(OpenBabel::OBMol  &mol,
				    bool               bVerbose,
				    int               *Nsymm,
				    int                Nrotbonds,
				    double             dbdt,
				    double            *temperature,
				    double            *DeltaHf0,
				    double            *DeltaHfT,
				    double            *DeltaGfT,
				    double            *DeltaSfT,
				    double            *S0T,
				    double            *CVT,
				    double            *CPT,
				    std::vector<double> &Scomponents,
				    double            *ZPVE);

}

#endif //DATA_UTILITIES_H

//! \file data_utilities.h
//! \brief Data related tools and utilities
