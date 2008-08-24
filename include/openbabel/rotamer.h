/**********************************************************************
rotamer.h - Handle rotamer list data.
 
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2008 Tim Vandermeersch
 
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

#ifndef OB_ROTAMER_H
#define OB_ROTAMER_H

#include <vector>
#include <map>

#include <openbabel/mol.h>
#include <openbabel/rotor.h>
#include <openbabel/generic.h>

namespace OpenBabel
{

  /** @class OBRotamerList rotamer.cpp <openbabel/rotamer.h> 
      @brief Supports a set of rotamer coordinate sets for some number of potentially rotatable bonds
      Further class introduction in rotamer.cpp
   */
  class OBAPI OBRotamerList : public OBGenericData
  {
    private:
      //! Number of atoms in the base coordinate set (i.e., OBMol::NumAtoms())
      unsigned int m_numAtoms;
      //! Base coordinates
      double*      m_coords;
      //! Individual bond rotors (from an OBRotor object or other)
      //! unsigned int* should point to 4 atom indexes for the rotor
      //! vector<int> is used to store the atoms that should be rotated for the rotor
      std::vector<std::pair<unsigned int*,std::vector<int> > > m_rotors;
      //! \brief Index of each rotor's different sampling states
      //! Usually from OBRotor::GetTorsionValues()
      std::vector<std::vector<double> >        m_torvals;
      //! Individual rotamer states (i.e., the array of rotor settings)
      std::vector<std::vector<unsigned char> > m_keys;
      
      /** 
       * Set up a rotamer list based on the supplied reference atoms and the number of rotors
       * @param mol The molecule to evaluate
       * @param ref An array of the 4 dihedral atoms for each rotor
       * @param nrotors The number of rotors (i.e., the size of ref / 4)
       */
      void Setup(OBMol &mol, unsigned int *ref, int nrotors);
      /**
       * Set all rotamer keys, does no error checking.
       */
      void SetRotamers(const std::vector<std::vector<unsigned char> > &keys) { m_keys = keys; }
      /**
       * @return A reference array (as used by AddRotamer() as a configuration of the individual rotor bonds
       */
      void GetReferenceArray(unsigned int*) const;
      /** 
       * Copies the coordinates, NOT the pointer, into this object
       * @param coords The conformer set for the molecule
       * @param N  The number of atoms in the molecule
       */
      void SetBaseCoordinates(double* coords, unsigned int N);

    public:
      /** 
       * Constructor
       */
      OBRotamerList();
      /**
       * Destructor.
       */
      ~OBRotamerList();
      /**
       * overload OBGenericData::Clone()
       */
      virtual OBGenericData* Clone(OBBase* parent) const;
      //! \name OBRotamerList modification methods
      //@{
      /**
       * Set up a rotamer list based on an already created OBRotorList. Setup() will copy
       * the molecule's current coordinates to be used as base coordinates. 
       */
      void Setup(OBMol&, OBRotorList&);
      /**
       * Add a rotamer to the list based on the supplied coordinate set as a double*. This
       * functions will calculate the torsion angles for the rotors and add these to internal
       * rotamer keys.
       * @param coords The coordinates from which a rotor key will be generated.
       * @return True if coords != 0.
       * 
       * This function can be used when you want to copy some rotamers (conformers) from your 
       * mol:
       * @code
       * ...
       * OBRotorList rl;
       * rl.Setup(mol)
       *
       * OBRotamerList rml;
       * rml.Setup(mol, rl);
       *
       * mol.SetConformer(3);
       * if (!rml.AddRotamer(mol.GetCoordinates()))
       *   cout << "Invalid key!" << endl;
       *  
       * mol.SetConformer(4);
       * if (!rml.AddRotamer(mol.GetCoordinates()))
       *   cout << "Invalid key!" << endl;
       *
       * // add some more rotamers
       * vector<int> key(rml.NumRotors(), 0);
       * if (!rml.AddRotamer(key))
       *   cout << "Invalid key!" << endl;
       *
       * key[2] = 1;
       * if (!rml.AddRotamer(key))
       *   cout << "Invalid key!" << endl;
       * ...
       * mol.ExpandConformerList(mol, mol.GetConformers());
       * @endcode
       */
      bool AddRotamer(double *coords);
      /**
       * Add a rotamer to the list based on the torsion angles in @p angles.
       */
      bool AddRotamer(std::vector<double> &angles);
      /** 
       * Add a rotamer to the list based on @p key as a configuration of the individual rotor bonds.
       * The key should have a size of NumRotors(). 
       *
       * The elements are indexes into the torsion values from OBRotor::GetTorsionValues().
       *
       * @param key The key for the rotamer.
       * @return True if the key is valid.
       */ 
      bool AddRotamer(std::vector<int> &key);
      /** 
       * Create a conformer list using the internal base set of coordinates.
       * @return The set of coordinates by rotating the bonds in each rotamer.
       */
      std::vector<double*> CreateConformerList(OBMol& mol);
      /** 
       * Create a conformer list using the internal base set of coordinates.
       * @return The set of coordinates as a reference in @p confs.
       */
      void ExpandConformerList(OBMol&mol, std::vector<double*> &confs);
      /**
       * Set the current coordinates for @p mol to the specified @p key.
       * @param mol The molecule.
       * @param key The rotamer key.
       */ 
      void SetCurrentCoordinates(OBMol &mol, std::vector<int> &key);
      //@}
      
      //! \name OBRotamerList data request methods
      //@{
      /**
       * @return The number of atoms in the base OBMol
       */
      unsigned int NumAtoms() const { return m_numAtoms; }
      /**
       * @return The number of rotatable bonds considered.
       */
      unsigned int NumRotors() const { return m_rotors.size(); }
      /** 
       * @return the number of rotamer (conformation) coordinate sets.
       */
      unsigned int NumRotamers() const { return m_keys.size(); }
      /** 
       * @return Pointer to a the base coordinates.
       */
      double* GetBaseCoordinates() const { return m_coords; }
      //@}
      
      //! \name Iterator methods
      //@{
      std::vector<std::vector<unsigned char> >::iterator BeginRotamer() { return m_keys.begin(); }
      std::vector<std::vector<unsigned char> >::iterator EndRotamer()   { return m_keys.end();   }
      //@}

  };

  /// @cond DEV
  class rotor_digit {
    public:
      rotor_digit(unsigned int rs) : resolution_size(rs), state(0) {}
      int get_state() { return state; }
      unsigned int size() { return resolution_size; }
      bool next()
      {
        if (state < (int)(resolution_size - 1)) {
          state++;
          return false;
        } else
          state = 0;
        
	return true;
      }
    private:
      unsigned int resolution_size;
      int state;
  } typedef rotor_digit;
  /// @endcond 
  
  /** @class OBRotamerKeys
   *  @brief A class to generate all possible rotamer keys
   */
  class OBAPI OBRotamerKeys
  {
      /** 
      @class OBRotamerKeys rotamer.h <openbabel/rotamer.h>
      @brief A class to generate all possible rotamer keys.

      This class can generate all possible rotor keys for a set of OBRotors 
      which can all have their own resolution. Thanks to Yongjin Xu for this
      patch. 

      the code blow is taken from  OBForceField::SystematicRotorSearch():
      \code
      #include <openbabel/rotamer.h>
      #include <openbabel/mol.h>

      // See OBConversion class to fill the mol object.
      OBMol mol;
      OBRotorList rl;
      OBRotamerList rotamers;

      rl.Setup(mol);
      rotamers.Setup(mol, rl);
    
      cout << "number of rotatable bonds: " <<  rl.Size() << endl;

      if (!rl.Size()) { // only one conformer
        cout << "generated only one conformer" << endl;
        // exit here 
      }

      OBRotamerKeys keys;
      OBRotorIterator ri;
      OBRotor *rotor = rl.BeginRotor(ri);
      for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
        rotorKeys.AddRotor(rotor->GetTorsionValues().size());
      }
    
      while (rotorKeys.Next()) {
        std::vector<int> rotorKey = rotorKeys.GetKey();
        rotamers.AddRotamer(rotorKey);
      }

      rotamers.ExpandConformerList(mol, mol.GetConformers());
      \endcode 
      **/
    public:
      //! Constructor
      OBRotamerKeys()
      {
        _vr.clear();
      }
      //! Clear all rotors
      void Clear(){
        _vr.clear();
      }
      //! Number of rotor keys (= number of possible conformers)
      unsigned int NumKeys()
      {
	unsigned int numKeys = 0;
	
	while (Next())
	  numKeys++;
	
	return numKeys; 
      }
      /**
       * Add a rotor
       * @param size the rotor resolution
       */
      void AddRotor(unsigned int size)
      {
        rotor_digit *rd;
        rd = new rotor_digit(size);
        _vr.push_back(*rd);
      }
      /**
       * Select the next rotor key
       * @return true if there are more rotor keys
       */
      bool Next()
      {
        if(_vr.size() == 0)
          return false;
        
	bool carry = _vr[0].next();
        unsigned int i = 1;
        while (carry) {
          if(i == _vr.size())
            return false;
          
	  carry = _vr[i].next();
          i++;
        }
        
        return true;
      }
      /** 
       * Get the currently selected rotor key
       * @return current rotor key
       */
      std::vector<int>& GetKey()
      {
        rt.clear();
        for(unsigned int i = 0; i < _vr.size(); i++){
          rt.push_back(_vr[i].get_state());
        }
        
	return rt;
      }
    private:
      std::vector<rotor_digit> _vr;
      std::vector<int> rt;
  };


}

#endif // OB_ROTAMER_H

//! \file rotamer.h
//! \brief Handle rotamer list data.
