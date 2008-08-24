/**********************************************************************
molchrg.h - Assign Gasteiger partial charges.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#ifndef OB_MOLCHRG_H
#define OB_MOLCHRG_H

namespace OpenBabel
{

  class GasteigerState;

  // class introduction in molchrg.cpp
  class OBAPI OBGastChrg
  {
    protected:
      std::vector <GasteigerState*> m_gsv; //!< vector of internal GasteigerState (for each atom)
      
      /** 
       * @brief Set initial partial charges in @p mol.
       * 
       * @param mol The molecule.
       *
       * Carbonyl O => -0.5
       * Phosphate O => -0.666
       * Sulfate O => -0.5
       * All other atoms are set to have their initial charge from their formal charge.
       */
      void InitialPartialCharges(OBMol &mol);
      /** 
       * @brief 
       * 
       * @param atom
       * @param a
       * @param b
       * @param c
       * 
       * @return 
       */
      bool GasteigerSigmaChi(OBAtom *atom, double &a, double &b, double &c);

    public:
      /** 
       * @brief Constructor.
       */
      OBGastChrg() {}
      /** 
       * @brief Destructor.
       */
      ~OBGastChrg();
      /** 
       * @brief Assign partial charges to this OBMol, setting each OBAtom partial charge.
       * 
       * @param mol The molecule
       */
      bool AssignPartialCharges(OBMol &mol);
      /** 
       * @brief Resize the internal GasteigerState vector to match the size of 
       * the molecule.
       * 
       * @param size The new size
       */
      void GSVResize(int size);
  };

  /** @class GasteigerState molchrg.h <openbabel/molchrg.h>
      @brief Helper class for OBGastChrg which stores the Gasteiger states 
      of a given atom
   */
  class OBAPI GasteigerState
  {
    public:
      /** 
       * @brief Constructor.
       */
      GasteigerState();
      /** 
       * @brief Destructor.
       */
      ~GasteigerState() {}
      /** 
       * @brief Set the values. 
       * 
       * @param _a
       * @param _b
       * @param _c
       * @param _q
       */
      void SetValues(double _a,double _b,double _c,double _q)
      {
        a = _a;
        b = _b;
        c = _c;
        denom=a+b+c;
        q = _q;
      }

      double a, b, c;
      double denom;
      double chi;
      double q;
  };

}

#define OB_GASTEIGER_DENOM  20.02
#define OB_GASTEIGER_DAMP   0.5
#define OB_GASTEIGER_ITERS  6

#endif // OB_MOLCHRG_H

//! @file molchrg.h
//! @brief Assign Gasteiger partial charges.
