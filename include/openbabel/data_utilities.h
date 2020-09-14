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
#include <mutex>

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

class OBAPI OBTranslator
{
	int _from, _to;

public:
	//! Constructor
	OBTranslator();
	OBTranslator(const char*, const char*);
	//! Destructor
	~OBTranslator() {};

	//! Set the initial atom type to be translated
	bool SetFromType(const char*);
	//! Set the destination atom type for translation
	bool SetToType(const char*);
	//! Translate atom types
	bool Translate(char *to, const char *from) const; // to, from
												//! Translate atom types
												//! \return whether the translation was successful
	bool Translate(std::string &to, const std::string &from) const; // to, from
															  //! Translate atom types
															  //! \return the translated atom type, or an empty string if not possible
	std::string Translate(const std::string &from) const;

	//! \return the initial atom type to be translated
	std::string GetFromType() const;
	//! \return the destination atom type for translation
	std::string GetToType() const;
};

class OBAPI OBResidueObserver
{
	int _resnum;

public:
	//! Sets the table to access the residue information for a specified
	//!  residue name
	//! \return whether this residue name is in the table
	bool SetResName(const std::string &);
	//! \return the bond order for the bond specified in the current residue
	//! \deprecated Easier to use the two-argument form
	int  LookupBO(const std::string &);
	//! \return the bond order for the bond specified between the two specified
	//! atom labels
	int  LookupBO(const std::string &, const std::string&);
	//! Look up the atom type and hybridization for the atom label specified
	//! in the first argument for the current residue
	//! \return whether the atom label specified is found in the current residue
	bool LookupType(const std::string &,std::string&,int&);
	//! Assign bond orders, atom types and residues for the supplied OBMol
	//! based on the residue information assigned to atoms
};

class OBGlobalMutex: public std::mutex
{
public:

#ifndef HAS_NOEXCEPT
 #if defined(__clang__)
  #if __has_feature(cxx_noexcept)
   #define HAS_NOEXCEPT
  #endif
 #else
  #if defined(__GXX_EXPERIMENTAL_CXX0X__) && __GNUC__ * 10 + __GNUC_MINOR__ >= 46 || \
  defined(_MSC_FULL_VER) && _MSC_FULL_VER >= 190023026
   #define HAS_NOEXCEPT
  #endif
 #endif
#endif

#ifdef HAS_NOEXCEPT
 #define NOEXCEPT noexcept
#else
 #define NOEXCEPT
#endif

	constexpr OBGlobalMutex() NOEXCEPT = default;

	OBGlobalMutex(const OBGlobalMutex &) { OBGlobalMutex(); } // Copy constructor

	~OBGlobalMutex() = default;

	OBGlobalMutex& operator=(const OBGlobalMutex &) { return *this; } // Copy assignment operator
};
}

#endif //DATA_UTILITIES_H

//! \file data_utilities.h
//! \brief Data related tools and utilities
