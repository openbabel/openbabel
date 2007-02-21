/**********************************************************************
rotor.h - Rotate torsional according to rotor rules.
 
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
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

#ifndef OB_ROTOR_H
#define OB_ROTOR_H

#include <openbabel/parsmart.h>
#include <openbabel/typer.h>

namespace OpenBabel
{

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

  //! \class OBRotorRule rotor.h <openbabel/rotor.h>
  //! \brief A rule for torsional conformer searching, defined by a SMARTS pattern
  //!
  //! Rules define a SMARTS pattern to match and a set of 4 reference atoms
  //! defining the dihedral angle. The rule can either define a set of possible
  //! dihedral angles in degrees and/or a "delta" (i.e., all multiples of delta will
  //! be considered)
  class OBAPI OBRotorRule
  {
    int                 _ref[4]; //!< Reference atoms specifying the dihedral angle (as integers), numbered from 1 inside the SMARTS pattern
    double              _delta;  //!< (optional) the resolution of a dihedral step in degrees
    std::string         _s;      //!< Text of the SMARTS pattern
    OBSmartsPattern*    _sp;     //!< The SMARTS pattern for the rotation rule
    std::vector<double> _vals;   //!< At least one torsion angle (in radians) to evaluate
  public:
  
  OBRotorRule(char *buffer,int ref[4],std::vector<double> &vals,double d):
    _delta(d), _s(buffer), _vals(vals)
    {
      _sp = new OBSmartsPattern;
      _sp->Init(buffer);
      memcpy(_ref,ref,sizeof(int)*4);
    }

    ~OBRotorRule()
      {
        if (_sp)
          {
            delete _sp;
            _sp = NULL;
          }
      }
  
    //! \return whether this rotor rule is valid (i.e., is the SMARTS pattern valid)
    bool    IsValid()    {        return(_sp->IsValid());       }
    //! \return a copy of the reference atom indexes inside the SMARTS pattern
    //! 
    //!  These should be freed after use.
    void    GetReferenceAtoms(int ref[4]) { memcpy(ref,_ref,sizeof(int)*4); }
    //! Set the resolution (delta) of a torsional step in degrees
    void    SetDelta(double d)    {       _delta = d;           }
    //! \return the resolution (delta) of a torsional step in degrees
    double  GetDelta()            {       return(_delta);       }
    //! \return a reference to the dihedral angles to evaluate (in radians)
    std::vector<double>   &GetTorsionVals()    { return(_vals); }
    //! \return the text of the SMARTS pattern for this rule
    std::string  &GetSmartsString(){      return(_s);           }
    //! \return the exact OBSmartsPattern object for this rule
    OBSmartsPattern *GetSmartsPattern() {  return(_sp);         }
  };

  //! \class OBRotorRules rotor.h <openbabel/rotor.h>
  //! \brief Database of default hybridization torsional rules and SMARTS-defined OBRotorRule objects
  //! 
  //! Use to automatically evaluate potentially rotatable bonds to generate
  //! lists of dihedral angles to consider.
  //! e.g., rotamer/conformer energy calculations
 class OBAPI OBRotorRules : public OBGlobalDataBase
  {
    bool                       _quiet;  //!< Control debugging output from GetRotorIncrements()
    std::vector<OBRotorRule*>  _vr;     //!< Database of specific OBRotorRules defined by SMARTS patterns
    std::vector<double>        _sp3sp3; //!< Default dihedral angles to check for generic sp3 - sp3 hybridized rotatable bonds (in radians)
    std::vector<double>        _sp3sp2; //!< Default dihedral angles to check for generic sp3 - sp2 hybridized rotatable bonds (in radians)
    std::vector<double>        _sp2sp2; //!< Default dihedral angles to check for generic sp2 - sp2 hybridized rotatable bonds (in radians)
  public:
    OBRotorRules();
    ~OBRotorRules();
  
    void ParseLine(const char*);
    //! \return the number of rotor rules
    unsigned int GetSize()                 { return _vr.size();}

    //! Set the filename to be used for the database. Default = torlib.txt
    void SetFilename(std::string &s)       { _filename = s;    }

    //! Determine the torsional angles to evaluate based on the database
    //! \param mol molecule to evaluate
    //! \param bond rotatable bond to evaluate
    //! \param refs set to be the atom indexes (in mol) of the dihedral angle
    //! \param vals set to be the list of angles to evaluate (in radians)
    //! \param delta potential dihedral angle steps (in degrees)
    void GetRotorIncrements(OBMol& mol,OBBond* bond,int refs[4],
                            std::vector<double> &vals,double &delta);
    //! Turn off debugging output from GetRotorIncrements()
    void Quiet()                           { _quiet=true;      }
  };

  //! \class OBRotor rotor.h <openbabel/rotor.h>
  //! \brief A single rotatable OBBond as part of rotamer searching
  class OBAPI OBRotor
  {
    int _idx,_ref[4];
    int *_rotatoms,_size,_numcoords;
    double _delta;
    double _imag,_refang;
    OBBond *_bond;
    std::vector<int> _torsion;
    OBBitVec _fixedatoms,_evalatoms;
    std::vector<double> _res;  //!< torsion resolution
    std::vector<double> _invmag;
    std::vector<std::vector<double> > _sn,_cs,_t;
  public:
    OBRotor();
    ~OBRotor()
      {
        if (_rotatoms)
          delete [] _rotatoms;
      }
    int     Size()
    {
      return((_res.empty())?0:_res.size());
    }
    int     GetIdx() const
    {
      return(_idx);
    }
    void    SetNumCoords(int nc)
    {
      _numcoords = nc;
    }
    void    SetBond(OBBond *bond)
    {
      _bond = bond;
    }
    void    SetEvalAtoms(OBBitVec &bv)
    {
      _evalatoms = bv;
    }
    void    SetDihedralAtoms(std::vector<int> &vi)
    {
      _torsion = vi;
    }
    void    SetDelta(double d)
    {
      _delta = d;
    }
    void    SetDihedralAtoms(int ref[4]);
    void    SetRotAtoms(std::vector<int>&);

    inline void SetToAngle(double *c,double setang)
    {
      double dx,dy,dz,sn,cs,t,ang,mag;
      ang = setang - CalcTorsion(c);
      if (fabs(ang) < 1e-5)
        return;

      sn = sin(ang);
      cs = cos(ang);
      t = 1 - cs;
      dx = c[_torsion[1]]   - c[_torsion[2]];
      dy = c[_torsion[1]+1] - c[_torsion[2]+1];
      dz = c[_torsion[1]+2] - c[_torsion[2]+2];
      mag = sqrt(SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
      Set(c,sn,cs,t,1.0/mag);
    }

    void    SetRotor(double *,int,int prev=-1);
    void    Set(double*,int);
    void    Precompute(double*);
    void    Set(double *c,int ridx,int cidx)
    {
      Set(c,_sn[cidx][ridx],_cs[cidx][ridx],_t[cidx][ridx],_invmag[cidx]);
    }
    void    Set(double*,double,double,double,double);
    void    Precalc(std::vector<double*>&);
    void    SetIdx(int idx)
    {
      _idx = idx;
    }
    void    SetFixedAtoms(OBBitVec &bv)
    {
      _fixedatoms = bv;
    }
    void    SetTorsionValues(std::vector<double> &tmp)
    {
      _res = tmp;
    }
    void    RemoveSymTorsionValues(int);
    void    GetDihedralAtoms(int ref[4])
    {
      for (int i=0;i<4;++i)
        ref[i]=_ref[i];
    }
    void    *GetRotAtoms()
    {
      return(_rotatoms);
    }
    double   CalcTorsion(double *);
    double   CalcBondLength(double*);
    double   GetDelta()
    {
      return(_delta);
    }
    OBBond *GetBond()
    {
      return(_bond);
    }
    std::vector<int> &GetDihedralAtoms()
      {
        return(_torsion);
      }
    std::vector<double> &GetResolution()
      {
        return(_res);
      }
    std::vector<double>::iterator BeginTorIncrement()
      {
        return(_res.begin());
      }
    std::vector<double>::iterator EndTorIncrement()
      {
        return(_res.end());
      }
    OBBitVec &GetEvalAtoms()
      {
        return(_evalatoms);
      }
    OBBitVec &GetFixedAtoms()
      {
        return(_fixedatoms);
      }
  };


  //! A standard iterator over a vector of rotors
  typedef std::vector<OBRotor*>::iterator OBRotorIterator;

  //! \class OBRotorList
  //! \brief Given an OBMol, set up a list of possibly rotatable torsions, 
  class OBAPI OBRotorList
  {
    bool _quiet;                    //!< Control debugging output
    bool _removesym;                //!< Control removal of symmetric rotations
    OBBitVec _fix;                  //!< Bit vector of fixed (i.e., invariant) atoms
    OBRotorRules _rr;               //!< Database of rotatable bonds and dihedral angles to test
    std::vector<int> _dffv;         //!< Distance from fixed
    std::vector<OBRotor*> _rotor;   //!< List of individual OBRotor torsions
    //!
    std::vector<std::pair<OBSmartsPattern*,std::pair<int,int> > > _vsym2;
    //! 
    std::vector<std::pair<OBSmartsPattern*,std::pair<int,int> > > _vsym3;
  public:

    OBRotorList();
    ~OBRotorList();

    //! Clear the internal list of rotors and reset
    void   Clear();

    //! \return the number of rotors in this list
    int    Size()
    {
      return((_rotor.empty()) ? 0: _rotor.size());
    }
    //! Intialize the private OBRotorRules database from a specific file
    void   Init(std::string &fname)
    {
      _rr.SetFilename(fname);
      _rr.Init();
    }
    //! Turn off debugging output
    void   SetQuiet()      { _quiet=true; _rr.Quiet();     }

    //! Set the list of fixed (invariant) atoms to the supplied OBBitVec
    void   SetFixAtoms(OBBitVec &fix) { _fix = fix;        }

    //! \brief Return whether this bond is fixed and thus not rotatable
    //! \return true if the bond and at least one neighboring bond has fixed atoms
    bool   IsFixedBond(OBBond*);
    //! \return whether this rotor list has any fixed (invariant) atoms
    bool   HasFixedAtoms()
    {
      return(!_fix.Empty());
    }

    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Insures the fixed atoms are respected, but otherwise functions like
    //! SetRotAtoms()
    void   SetRotAtomsByFix(OBMol&);

    //! \brief Set the atoms to rotate from the dihedral atoms for each rotor
    //! Uses OBRotor->GetDihedralAtoms() to call OBRotor->SetRotAtoms()
    //! and standarizes the dihedral angles via OBRotor->SetDihedralAtoms()
    //! \return True
    bool   SetRotAtoms(OBMol&);

    //! Setup this rotor list for the supplied molecule
    //! Calls FindRotors(), SetEvalAtoms(), and AssignTorVals()
    //! \return True if rotatable bonds were found
    bool   Setup(OBMol &);
    //! Find all potentially rotatable bonds in the molecule
    //! Uses OBBond::IsRotor() for initial evaluation
    //! \return True
    bool   FindRotors(OBMol &);
    //! Determines which atoms should be used to calculate the internal energy
    //! if the dihedral angle of the rotor is modified
    //! \return True
    bool   SetEvalAtoms(OBMol&);
    //! Using the OBRotorRules database, set the torsion values (and delta)
    //! to be evaluated and tested
    bool   AssignTorVals(OBMol &);

    //! Has no effect
    //! \deprecated Currently has no effect
    void   IgnoreSymmetryRemoval()    { _removesym = false;}
    //! Rotates each bond to zero and 180 degrees and tests
    //! if the 2 conformers are duplicates.  if so - the symmetric torsion
    //! values are removed from consideration during a search
    void   RemoveSymVals(OBMol&);

    //! \name Iterator methods
    //@{
    OBRotor *BeginRotor(OBRotorIterator &i)
      { i = _rotor.begin(); return((i ==_rotor.end()) ? NULL:*i); }
    OBRotor *NextRotor(OBRotorIterator &i)
      { ++i; return((i ==_rotor.end()) ? NULL:*i); }
    OBRotorIterator BeginRotors()   { return(_rotor.begin()); }
    OBRotorIterator EndRotors()     { return(_rotor.end());   }
    //@}

    // Not declared
    //! \deprecated Not declared. Use Setup() for top-level functionality
    bool   IdentifyEvalAtoms(OBMol &mol) { return SetEvalAtoms(mol); } 
  };

} // end namespace OpenBabel

#endif // OB_ROTOR_H

//! \file rotor.h
//! \brief Rotate torsional according to rotor rules.
