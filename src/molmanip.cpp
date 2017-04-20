/********************************************************************** 
molmanip.cpp - Manipulate the geometries of molecules.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
Some portions Copyright (C) 2017 by Torsten Sachse

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

#include <set>

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <openbabel/molmanip.h>

#define PI 3.1415926535897932384626433832795

using namespace std;

namespace OpenBabel
{

  /** \class OBMolManip molmanip.h <openbabel/molmanip.h>
      \brief Molecule Manipulation Class

      The OBMolManip class is meant to streamline modifications to the
      geometry of molecules, such as appending another molecule, bonding
      two molecules together, and setting and computing angles.
  */

  //
  // OBMolManip member functions
  //

  //! \brief Set a torsion angle (i.e., dihedral angle).
  //!
  //! Set the dihedral angle around the bond between atoms 2 and 3 to angle in degrees.
  //! The direction is defined by atoms 1 and 4. A dihedral angle of 0 means a
  //! cis-configuration whereas PI means a trans configuration.
  bool OBMolManip::SetDihedralAngle(const int idxa1, const int idxa2, const int idxa3, const int idxa4, const double angle)
  {
    MOLMANIPSANITY(false)
    double ang = angle*DEG_TO_RAD;
    OBBond *bond;
    OBAtom *a1, *a2, *a3, *a4;

    a1=_manip->GetAtom(idxa1);
    a2=_manip->GetAtom(idxa2);
    a3=_manip->GetAtom(idxa3);
    a4=_manip->GetAtom(idxa4);

    bond=a2->GetBond(a3);

    if (bond == NULL){cerr << "Warning: atoms " << idxa2 << " and " << idxa3 << " form no bond.\n";};
    if (!(bond->IsRotor())){cerr << "Warning: bond between atoms " << idxa2 << " and " << idxa3 << " is no rotor.\n";};

    vector<int> tor;
    vector<int> atoms;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetDihedralAngle", obAuditMsg);

    tor.push_back(idxa1);
    tor.push_back(idxa2);
    tor.push_back(idxa3);
    tor.push_back(idxa4);

    _manip->FindChildren(atoms, a2->GetIdx(), a3->GetIdx());
    int j;
    for (j = 0 ; (unsigned)j < atoms.size() ; j++ )
      atoms[j] = (atoms[j] - 1) * 3;

    double v2x,v2y,v2z;
    double radang,m[9];
    double x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

    //calculate the torsion angle
    radang = CalcTorsionAngle(a1->GetVector(),
                              a2->GetVector(),
                              a3->GetVector(),
                              a4->GetVector()) / RAD_TO_DEG;
    //
    // now we have the torsion angle (radang) - set up the rot matrix
    //

    //find the difference between current and requested
    rotang = ang - radang;

    double* _c = _manip->_c;

    sn = sin(rotang);
    cs = cos(rotang);
    t = 1 - cs;

    v2x = _c[tor[1]]   - _c[tor[2]];
    v2y = _c[tor[1]+1] - _c[tor[2]+1];
    v2z = _c[tor[1]+2] - _c[tor[2]+2];

   //normalize the rotation vector
    mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
    x = v2x/mag;
    y = v2y/mag;
    z = v2z/mag;

    //set up the rotation matrix
    m[0]= t*x*x + cs;
    m[1] = t*x*y + sn*z;
    m[2] = t*x*z - sn*y;
    m[3] = t*x*y - sn*z;
    m[4] = t*y*y + cs;
    m[5] = t*y*z + sn*x;
    m[6] = t*x*z + sn*y;
    m[7] = t*y*z - sn*x;
    m[8] = t*z*z + cs;

    //
    //now the matrix is set - time to rotate the atoms
    //
    tx = _c[tor[1]];
    ty = _c[tor[1]+1];
    tz = _c[tor[1]+2];
    vector<int>::iterator i;
    for (i = atoms.begin(); i != atoms.end(); ++i)
      {
        j = *i;

        _c[j] -= tx;
        _c[j+1] -= ty;
        _c[j+2]-= tz;
        x = _c[j]*m[0] + _c[j+1]*m[1] + _c[j+2]*m[2];
        y = _c[j]*m[3] + _c[j+1]*m[4] + _c[j+2]*m[5];
        z = _c[j]*m[6] + _c[j+1]*m[7] + _c[j+2]*m[8];
        _c[j] = x;
        _c[j+1] = y;
        _c[j+2] = z;
        _c[j] += tx;
        _c[j+1] += ty;
        _c[j+2] += tz;
      }
  }


  double OBMolManip::GetDihedralAngle(const int idxa1, const int idxa2, const int idxa3, const int idxa4)
  {
    MOLMANIPSANITY(360.0)
    OBAtom *a1, *a2, *a3, *a4;

    a1=_manip->GetAtom(idxa1);
    a2=_manip->GetAtom(idxa2);
    a3=_manip->GetAtom(idxa3);
    a4=_manip->GetAtom(idxa4);

    return(CalcTorsionAngle(a1->GetVector(),
                            a2->GetVector(),
                            a3->GetVector(),
                            a4->GetVector()));
  }

  //! \brief Bond to molecules together
  //!
  //! One bond in each molecule is broken (i1-i2 and m1-m2, respectively)
  //! and i1 and m1 are then bonded together. The bonds are properly aligned
  //! beforehand.
  bool OBMolManip::BondMolecule(OBMol &mol, int i1, int i2, int m1, int m2){
    MOLMANIPSANITY(false)
    OBMol bondmol;;
    if (_manip == &mol){
      bondmol = OBMol(mol);
    }
    else{
      bondmol = (OBMol &)(mol);
    }
    OBMolManip bondmanip(&bondmol);
    //Rotate the to-be-bonded molecule so that the two bonds (i1-i2 and m2-m2) align.
    //  Obtain specific points in space and specific directions for first molecule
    OBAtom *a1, *a2;
    a1 = _manip->GetAtom(i1);
    a2 = _manip->GetAtom(i2);
    vector3 bond_dir1    = (a1->GetVector() - a2->GetVector());
    vector3 bond_center1 = (a1->GetVector() + a2->GetVector()) / 2.0;
    //  Extract bond information
    OBBond* bond = _manip->GetBond(i1,i2);
    if (!bond){
      std::cerr << "ERROR: atoms " << i1 << " and " << i2 << " of primary molecule form no bond." << std::endl;
      return false;
    }
    unsigned int order = bond->GetBondOrder();
    unsigned int flags = bond->GetFlags();
    //  Obtain specific points in space and specific directions for second molecule
    a1 = bondmol.GetAtom(m1);   a2 = bondmol.GetAtom(m2);
    vector3 bond_dir2    = (a1->GetVector() - a2->GetVector());
    vector3 bond_center2 = (a1->GetVector() + a2->GetVector()) / 2.0;
    //  Check whether both bonds have the same flags and bond order
    bond = bondmol.GetBond(m1,m2);
    if (!bond){
      std::cerr << "ERROR: atoms " << m1 << " and " << m2 << " of secondary molecule form no bond." << std::endl;
      return false;
    }
    if (bond->GetBondOrder() != order){
      std::cerr << "WARNING: bond orders of the to-be-split bonds do not agreee, will reuse value from primary one." << std::endl;
      std::cerr << "         However, I recommend writing this to file and reading it back in again." << std::endl;
    }
    if (bond->GetFlags() != flags){
      std::cerr << "WARNING: flas of the to-be-split bonds do not agreee, will reuse value from primary one." << std::endl;
      std::cerr << "         However, I recommend writing this to file and reading it back in again." << std::endl;
    }
    //  Only proceed if both vectors involved are OK
    if (!bond_dir1.CanBeNormalized()) {
      std::cerr << "ERROR: cannot get bond direction of the primary molecule." << std::endl;
      return false;
    }
    if (!bond_dir2.CanBeNormalized()) {
      std::cerr << "ERROR: cannot get bond direction of the to-be-bonded molecule." << std::endl;
      return false;
    }
    // find angles and axes by which to rotate the to-be-bonded molecule
    bond_dir1 = bond_dir1.normalize();
    bond_dir2 = bond_dir2.normalize();
    vector3 tempvec = cross(bond_dir1, bond_dir2);
    double angle;
    if (tempvec.CanBeNormalized())
      angle = acos(dot(bond_dir1,bond_dir2))*RAD_TO_DEG - 180.0;
    else{
      tempvec = cross(vector3(1.0,0.0,0.0),bond_dir1);
      if (!(tempvec.CanBeNormalized()))
        tempvec = cross(vector3(0.0,1.0,0.0),bond_dir1);
      if (dot(bond_dir1,bond_dir2) > 0.0)
        angle = 180.0;
      else
        angle = 0.0;
    }
    // rotate (therefore a proper displacement is required since rotation is always aroung an axes through the origin)
    bondmanip.Translate(-bond_center2);
    bondmanip.Rotate(tempvec, angle);
    bondmanip.Translate(bond_center1);
    //save the atoms that shall be bonded afterwards
    a1 = _manip->GetAtom(i1);
    a2 = bondmol.GetAtom(m1);
    //Cleave off the bits of each molecule that should be removed when bonding them.
    if (!CleaveOff(i1,i2)){
      std::cerr << "ERROR cleaving off part of the primary molecule." << std::endl;
      return false;
    }
    if (!bondmanip.CleaveOff(m1,m2)){
      std::cerr << "ERROR cleaving off part of the to-be-bonded molecule." << std::endl;
      return false;
    }
    //Copy over atom positions and bonds from bondmol
    int nr_atoms = _manip->NumAtoms();
    AppendMolecule(bondmol);
    //Add final bond
    unsigned int begin_atom_idx = a1->GetIdx();
    unsigned int end_atom_idx = a2->GetIdx() + nr_atoms;
    _manip->AddBond(begin_atom_idx, end_atom_idx, order, flags, -1);
    return true;
  }

  //! \brief Copy the data in a molecule over to this one.
  void OBMolManip::AppendMolecule(OBMol &src)
  {
    MOLMANIPSANITY(void())
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    OBAtom *atom, *a;
    OBBond *bond;
    unsigned int at_count=1;
    const unsigned int nr_atoms=_manip->NumAtoms();
    const unsigned int nr_bonds=_manip->NumBonds();
    unsigned int begin_atom_idx, end_atom_idx, order, flags;

    _manip->BeginModify();
    _manip->_vatom.reserve(nr_atoms+src.NumAtoms());
    _manip->_atomIds.reserve(nr_atoms+src.NumAtoms());
    _manip->_vbond.reserve(nr_bonds+src.NumBonds());
    _manip->_bondIds.reserve(nr_bonds+src.NumBonds());

    for (atom = src.BeginAtom(i);atom;atom = src.NextAtom(i), ++at_count)
    {
      a = _manip->CreateAtom();
      a->Duplicate(atom);
      _manip->AddAtom(*a,true);
    }
  
    for (bond = src.BeginBond(j);bond;bond = src.NextBond(j))
    {
      begin_atom_idx = bond->GetBeginAtomIdx() + nr_atoms;
      end_atom_idx = bond->GetEndAtomIdx() + nr_atoms;
      order = bond->GetBondOrder();
      flags = bond->GetFlags();
      _manip->AddBond(begin_atom_idx, end_atom_idx, order, flags, -1);
    }
    _manip->EndModify();
  }

  //! \brief Remove atoms and bonds on one side of a bond.
  //!
  //! Find all atoms that are to be cleaved off.
  //! Those are connected to the atom indexed by i2.
  //! That atom will also be removed byt i1 will stay.
  bool OBMolManip::CleaveOff(int i1, int i2){
    MOLMANIPSANITY(false)
    std::vector<OBAtom*> cleave_atoms;
    OBAtom *a1,*a2;
    a1 = _manip->GetAtom(i1);
    a2 = _manip->GetAtom(i2);
    _manip->FindChildren(cleave_atoms, a1 , a2);
    //FindChildren does not include a1 or a2 in the vector
    cleave_atoms.push_back(a2);
    if (cleave_atoms.size()+1 >= _manip->NumAtoms()){
      std::cerr << "ERROR: the cleave-off requested would remove the entire molecule: i1: " << i1 << ", i2: " << i2 << std::endl;
      return false;
    }
    //Find all bonds that are to be removed.
    std::set<OBBond*> cleave_bonds;
    for (std::vector<OBAtom*>::iterator atomit = cleave_atoms.begin(); atomit!=cleave_atoms.end(); ++atomit){
      OBBondIterator bondit;
      OBBond* b = (*atomit)->BeginBond(bondit);
      while (b){
        if (cleave_bonds.find(b) == cleave_bonds.end()){
          cleave_bonds.insert(b);
        }
        b = (*atomit)->NextBond(bondit);
      }
    }
    //Remove bonds
    _manip->BeginModify();
    for (std::set<OBBond*>::iterator bondit = cleave_bonds.begin(); bondit!=cleave_bonds.end(); ++bondit){
      //Delete and destroy bond
      _manip->DeleteBond(*bondit,true);
    }
    _manip->EndModify();
    //Remove atoms 
    _manip->BeginModify();
    for (std::vector<OBAtom*>::iterator atomit = cleave_atoms.begin(); atomit!=cleave_atoms.end(); ++atomit){
      //Delete and destroy atom
      _manip->DeleteAtom(*atomit,true);
    }
    _manip->EndModify();
    return true;
  }

  //! \brief Partition a molecule by a plane.
  //!
  //! Copy all atoms on one side (distance in Hessian normal form >0)
  //! of a plane to another molecule.
  void OBMolManip::PartMolecule(OBMol &dest,  const double direction[3], const double point[3])
  {
    MOLMANIPSANITY(void())
    vector<OBAtom*>::iterator i;
    OBAtom *atom, *a;
    unsigned int at_count=1;

    for (atom = _manip->BeginAtom(i);atom;atom = _manip->NextAtom(i))
    {
      double x=atom->GetX();
      double y=atom->GetY();
      double z=atom->GetZ();
      if ( (x-point[0])*direction[0] + (y-point[1])*direction[1] + (z-point[2])*direction[2] >= 0 ){
        dest.AddAtom(*atom,true);
        a=dest.GetAtom(at_count);
        a->Duplicate(atom);
        ++at_count;
      }
    }
  }

  //! \brief Rotate a molecule around a given axis (3 doubles) by an angle in degrees.
  void OBMolManip::Rotate(const double axis[3], const double angle)
  {
    MOLMANIPSANITY(void())
    matrix3x3 matrix;
    matrix.RotAboutAxisByAngle(vector3(axis[0],axis[1],axis[2]),angle);
    Rotate(_manip->_vatom, matrix);
  }

  //! \brief Rotate a molecule around a given axis (vector3) by an angle in degrees.
  void OBMolManip::Rotate(const vector3 axis, const double angle)
  {
    MOLMANIPSANITY(void())
    matrix3x3 matrix;
    matrix.RotAboutAxisByAngle(axis,angle);
    Rotate(_manip->_vatom, matrix);
  }

  //! \brief Rotate a subset of atoms of a molecule using a rotation matrix.
  void OBMolManip::Rotate(std::vector<OBAtom*> &atom_vec, const matrix3x3 m)
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    vector3 tempvec;

    for (i = atom_vec.begin(); i != atom_vec.end(); ++i)
    {
      if (!(*i)){
        break;
      }
      atom = *i;
      tempvec = atom->GetVector();
      tempvec *= m;
      atom->SetVector(tempvec);
    }
  }

  //! \brief Rotate a molecule around one of its three main axes by an angle in degrees.
  void OBMolManip::Rotate(const int main_axis_nr, const double angle)
  {
    MOLMANIPSANITY(void())
    Rotate(_manip->_vatom, main_axis_nr, angle);
  }

  //! \brief Rotate a subset of atoms of a molecule one of its three main axes by an angle in degrees.
  void OBMolManip::Rotate(std::vector<OBAtom*> &atom_vec, const int main_axis_nr, const double angle)
  {
    if(main_axis_nr > 3 || main_axis_nr < 1){cerr << "there are only three main axes."; return;};
    vector3 axis, tempvec, eigenvals;
    matrix3x3 matrix, tempmatrix;
    vector3 center = vector3(0,0,0);
    vector<OBAtom*>::iterator i;
    OBAtom *atom;
    double d[3][3];
    double x,y,z;

    int elements = 0;
    for (i = atom_vec.begin(); i != atom_vec.end(); ++i)
    {
      if (!(*i)){
        break;
      }
      ++elements;
      atom = *i;
      center += atom->GetVector();
    }
    center /= elements;
    d[0][0]=0;d[0][1]=0;d[1][1]=0;
    d[0][2]=0;d[1][2]=0;d[2][2]=0;
    d[1][0]=0;d[2][0]=0;d[2][1]=0;

    for (i = atom_vec.begin(); i != atom_vec.end(); ++i)
    {
      if (!(*i)){
        break;
      }
      atom = *i;
      tempvec=atom->GetVector()-center;
      x=tempvec.GetX();
      y=tempvec.GetY();
      z=tempvec.GetZ();
      d[0][0]+=y*y+z*z;
      d[0][1]-=x*y;
      d[1][1]+=x*x+z*z;
      d[0][2]-=x*z;
      d[1][2]-=y*z;
      d[2][2]+=x*x+y*y;
    }

    d[1][0]=d[0][1];
    d[2][0]=d[0][2];
    d[2][1]=d[1][2];

    matrix=matrix3x3(d);

    tempmatrix = matrix.findEigenvectorsIfSymmetric(eigenvals);

    axis = tempmatrix.GetColumn(3-main_axis_nr);
    matrix.RotAboutAxisByAngle(axis,angle);
    Rotate(atom_vec, matrix);
  }

  //! \brief Translate a molecule by a vector (thee doubles)
  void OBMolManip::Translate(const double vector[3])
  {
    const vector3 v = vector3(vector[0],vector[1],vector[2]);
    Translate(vector3(v));
  }

  //! \brief Translate a subset of atoms of a molecule by a vector (thee doubles)
  void OBMolManip::Translate(std::vector<OBAtom*> &atom_vec, const double vec[3])
  {
    const vector3 v = vector3(vec[0],vec[1],vec[2]);
    OBAtom* atom;
    std::vector<OBAtom*>::iterator i;
    vector3 tempvec;

    int elements = 0;
    for (i = atom_vec.begin(); i != atom_vec.end(); ++i)
    {
      if (!(*i)){
        break;
      }
      ++elements;
      atom = *i;
      tempvec = atom->GetVector();
      tempvec += v;
      atom->SetVector(tempvec);
    }
  }

  //! \brief Align a center of a molecule as well as its third and second main axes.
  //!
  //! The molecule's center is first aligned to a point and it's then rotated so that
  //! its third and second main axes are aligned in the specified directions.
  //! All coordinates have to be given as triples of doubles.
  void OBMolManip::Align(const double p[3], const double v1[3], const double v2[3])
  {
    MOLMANIPSANITY(void())
    Align(vector3(p[0],p[1],p[2]),vector3(v1[0],v1[1],v1[2]),vector3(v2[0],v2[1],v2[2]));
  }

  //! \brief Align a center of a molecule as well as its third and second main axes.
  //!
  //! The molecule's center is first aligned to a point and it's then rotated so that
  //! its third and second main axes are aligned in the specified directions.
  //! All coordinates have to be given as vector3.
  void OBMolManip::Align(const vector3 p, const vector3 v1, const vector3 v2)
  {
    MOLMANIPSANITY(void())
    Align(_manip->_vatom, p, v1, v2);
  }

  //! \brief Align the center of a subset of atoms of a molecule as well as the associated third and second main axes.
  //!
  //! The molecule's center is first aligned to a point and it's then rotated so that
  //! its third and second main axes are aligned in the specified directions.
  //! All coordinates have to be given as vector3.
  void OBMolManip::Align(vector<OBAtom*> &atom_vec, const vector3 p, const vector3 v1, const vector3 v2)
  {
    matrix3x3 matrix, tempmatrix; 
    vector3 eigenvals, tempvec, main1, main2;
    vector3 center = vector3(0,0,0);
    double d[3][3];
    double x,y,z,angle;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    int elements = 0;
    for (i = atom_vec.begin(); *i != NULL && i != atom_vec.end(); ++i)
    {
      center += (*i)->GetVector();
      ++elements;
    }
    center /= elements;
    d[0][0]=0;d[0][1]=0;d[1][1]=0;
    d[0][2]=0;d[1][2]=0;d[2][2]=0;
    d[1][0]=0;d[2][0]=0;d[2][1]=0;

    for (i = atom_vec.begin(); *i != NULL && i != atom_vec.end(); ++i)
    {
      atom = *i;
      tempvec=atom->GetVector()-center;
      x=tempvec.GetX();
      y=tempvec.GetY();
      z=tempvec.GetZ();
      d[0][0]+=y*y+z*z;
      d[0][1]-=x*y;
      d[1][1]+=x*x+z*z;
      d[0][2]-=x*z;
      d[1][2]-=y*z;
      d[2][2]+=x*x+y*y;
    }
    d[1][0]=d[0][1];
    d[2][0]=d[0][2];
    d[2][1]=d[1][2];
    
    matrix=matrix3x3(d);

    tempmatrix = matrix.findEigenvectorsIfSymmetric(eigenvals);

    main1 = tempmatrix.GetColumn(0);
    tempvec = cross(v1,main1);
    if ( tempvec.CanBeNormalized() ) {
      angle = vectorAngle(v1,main1);
      matrix.RotAboutAxisByAngle(tempvec,angle);
    }
    else{
      matrix.RotAboutAxisByAngle(main1,0.0);
    }

    main2 = tempmatrix.GetColumn(1);
    main2 *= matrix;
    tempvec = cross(v2,main2);
    if ( tempvec.CanBeNormalized() ) {
      angle = vectorAngle(v2,main2);
      tempmatrix.RotAboutAxisByAngle(tempvec,angle);
    }
    else{
      tempmatrix.RotAboutAxisByAngle(main2,0.0);
    }

    matrix=tempmatrix*matrix;
    
    for (i = atom_vec.begin(); *i != NULL && i != atom_vec.end(); ++i)
    {
      atom = *i;
      tempvec = atom->GetVector();
      tempvec -= center;
      tempvec *= matrix;
      tempvec += p;
      atom->SetVector(tempvec);
    }
    
  }

  //! \brief Set an angle.
  //!
  //! Change the bond defined by the three atom indices to the value given by angle.
  //! The first atom is kept fixed and everything attached to the second one is moved.
  bool OBMolManip::SetAngle(const int idxa1, const int idxa2, const int idxa3, const double angle)
  {
    MOLMANIPSANITY(false)
    vector<int> children;
    vector3 tempvec, translation, v1, v2, n;
    matrix3x3 m;
    OBBond *bond1, *bond2;
    OBAtom *a1, *a2, *a3, *atom;
    vector<int>::iterator i;
    double currentangle;

    if (idxa1 == idxa2 || idxa1 == idxa3 || idxa2 == idxa3 ){cerr << "some atoms are the same\a";return false;};

    a1=_manip->GetAtom(idxa1);
    a2=_manip->GetAtom(idxa2);
    a3=_manip->GetAtom(idxa3);
    v1=a1->GetVector()-a2->GetVector();
    v2=a3->GetVector()-a2->GetVector();

    bond1=a2->GetBond(a1);
    bond2=a2->GetBond(a3);

    if (bond1 == NULL){cerr << "Warning: atoms " << idxa1 << " and " << idxa2 << " form no bond.\n";};
    if (bond2 == NULL){cerr << "Warning: atoms " << idxa2 << " and " << idxa3 << " form no bond.\n";};

    n=cross(v1,v2);
    currentangle=vectorAngle(v1,v2);
#define streched 180
    if (currentangle == streched){cerr << "angle is 180deg, changing it would be arbitrary\a";return false;}
#undef streched

    //find which atoms to rotate
    _manip->FindChildren(children,a2->GetIdx(),a3->GetIdx());
    children.push_back(a3->GetIdx());

    if (children.size() == _manip->NumAtoms()-1){cerr << "all atoms connected\n";return false;};

    translation=a2->GetVector();

    m.RotAboutAxisByAngle(n,currentangle-angle);
    
    _manip->BeginModify();
    //rotate atoms
    for (i = children.begin();i != children.end();++i)
    {
      atom = _manip->GetAtom(*i);
      tempvec = atom->GetVector();
      tempvec-=translation;
      tempvec *= m;   //rotate the point
      tempvec+=translation;
      atom->SetVector(tempvec);
    }
    _manip->EndModify();
    return true;
  }

  //! \brief Get an angle.
  //!
  //! The bond is defined by the three atom indices.
  double OBMolManip::GetAngle(const int idxa1, const int idxa2, const int idxa3)
  {
    MOLMANIPSANITY(false)
    vector3 v1, v2, n;
    OBAtom *a1, *a2, *a3;
    vector<int>::iterator i;
    double currentangle;

    if (idxa1 == idxa2 || idxa1 == idxa3 || idxa2 == idxa3 ){cerr << "some atoms are the same\n";return false;};

    a1=_manip->GetAtom(idxa1);
    a2=_manip->GetAtom(idxa2);
    a3=_manip->GetAtom(idxa3);
    v1=a1->GetVector()-a2->GetVector();
    v2=a3->GetVector()-a2->GetVector();

    n=cross(v1,v2);
    currentangle=vectorAngle(v1,v2);

    return currentangle;
  }

  //! \brief Mirror a molecule by a plane or point using triples of doubles.
  //!
  //! The mirrir plane is defined by a point in the plane and a normal vector. If
  //! the normal vector's legths vanishes, point inversion is performed.
  //! If center_it==true, the molecule's center is temporarily moved to the origin
  //! before applying the transformation and moved to its original position afterwards.
  void OBMolManip::Mirror(const double normal[3], const double point[3], bool center_it)
  {
    Mirror(vector3(normal[0],normal[1],normal[2]),vector3(point[0],point[1],point[2]), center_it);
  }

  //! \brief Mirror a molecule by a plane or point using vector3.
  //!
  //! The mirrir plane is defined by a point in the plane and a normal vector. If
  //! the normal vector's legths vanishes, point inversion is performed.
  //! If center_it==true, the molecule's center is temporarily moved to the origin
  //! before applying the transformation and moved to its original position afterwards.
  void OBMolManip::Mirror(const vector3 normal, const vector3 point, bool center_it){
    MOLMANIPSANITY(void())
    Mirror(_manip->_vatom, normal, point, center_it);
  }

  //! \brief Mirror a subset of atoms of a molecule by a plane or point using vector3.
  //!
  //! The mirrir plane is defined by a point in the plane and a normal vector. If
  //! the normal vector's legths vanishes, point inversion is performed.
  //! If center_it==true, the molecule's center is temporarily moved to the origin
  //! before applying the transformation and moved to its original position afterwards.
  void OBMolManip::Mirror(std::vector<OBAtom*> &atom_vec, const vector3 normal, const vector3 point, bool center_it)
  {
    
    double dist, pointdist;
    vector3 center = vector3(0,0,0);
    vector3 tempvec, norm;
    matrix3x3 matrix; 
    double d[3][3];

    vector<OBAtom*>::iterator i;
    OBAtom *atom;

    pointdist = (point[0] * normal[0] +
          point[1] * normal[1] +
          point[2] * normal[2]);

    if (center_it){
      int elements = 0;
      for (i = atom_vec.begin(); *i != NULL && i != atom_vec.end(); ++i)
      {
        atom = *i;
        center += atom->GetVector();
        ++elements;
      }
      center /= elements;
    }

    if ( !(normal.CanBeNormalized()) ){
      center += point;
      for (i = atom_vec.begin(); *i != NULL && i != atom_vec.end(); ++i)
      {
        atom = *i;
        tempvec = atom->GetVector();
        tempvec -= center;
        tempvec *= -1;
        tempvec += center;
        atom->SetVector(tempvec);
      }
    }
    else{
      norm = normal / sqrt( (normal[0] * normal[0]) +
                 (normal[1] * normal[1]) +
                 (normal[2] * normal[2])
                );
      double a,b,c;
      a = norm[0];
      b = norm[1];
      c = norm[2];

      d[0][0] = 1-(2*a*a);
      d[1][1] = 1-(2*b*b);
      d[2][2] = 1-(2*c*c);
      d[1][0] = -2*a*b; d[0][1] = d[1][0];
      d[2][0] = -2*a*c; d[0][2] = d[2][0];
      d[2][1] = -2*b*c; d[1][2] = d[2][1];
      matrix=matrix3x3(d);

      for (i = atom_vec.begin(); *i != NULL && i != atom_vec.end(); ++i)
      {
        atom = *i;
        tempvec = atom->GetVector();
        tempvec -= center;
        tempvec -= pointdist*norm;
        tempvec *= matrix;
        tempvec += pointdist*norm;
        tempvec += center;
        atom->SetVector(tempvec);
      }
    }
  }

  /*! this method adds the vector v to all atom positions in all conformers */
  void OBMolManip::Translate(const vector3 &v)
  {
    MOLMANIPSANITY(void())
    for (int i = 0;i < _manip->NumConformers();++i)
      Translate(v,i);
  }

  /*! this method adds the vector v to all atom positions in the
    conformer nconf. If nconf == OB_CURRENT_CONFORMER, then the atom
    positions in the current conformer are translated. */
  void OBMolManip::Translate(const vector3 &v, int nconf)
  {
    MOLMANIPSANITY(void())
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Translate", obAuditMsg);

    int i,size;
    double x,y,z;
    double *c = (nconf == OB_CURRENT_CONFORMER)? _manip->_c : _manip->GetConformer(nconf);

    x = v.x();
    y = v.y();
    z = v.z();
    size = _manip->NumAtoms();
    for (i = 0;i < size;++i)
      {
        c[i*3  ] += x;
        c[i*3+1] += y;
        c[i*3+2] += z;
      }
  }

  void OBMolManip::Rotate(const double u[3][3])
  {
    int i,j,k;
    double m[9];
    for (k=0,i = 0;i < 3;++i)
      for (j = 0;j < 3;++j)
        m[k++] = u[i][j];

    for (i = 0;i < _manip->NumConformers();++i)
      Rotate(m,i);
  }

  void OBMolManip::Rotate(const double m[9])
  {
    MOLMANIPSANITY(void())
    for (int i = 0;i < _manip->NumConformers();++i)
      Rotate(m,i);
  }

  void OBMolManip::Rotate(const double m[9],int nconf)
  {
    MOLMANIPSANITY(void())
    int i,size;
    double x,y,z;
    double *c = (nconf == OB_CURRENT_CONFORMER)? _manip->_c : _manip->GetConformer(nconf);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Rotate", obAuditMsg);

    size = _manip->NumAtoms();
    for (i = 0;i < size;++i)
      {
        x = c[i*3  ];
        y = c[i*3+1];
        z = c[i*3+2];
        c[i*3  ] = m[0]*x + m[1]*y + m[2]*z;
        c[i*3+1] = m[3]*x + m[4]*y + m[5]*z;
        c[i*3+2] = m[6]*x + m[7]*y + m[8]*z;
      }
  }

  void OBMolManip::AddConformer(double* f, bool deepcopy)
  {
      MOLMANIPSANITY(void())
      if (deepcopy){
          double* new_f = new double[3*_manip->NumAtoms()];
          memcpy((char*)new_f,(char*)f,sizeof(double)*3*_manip->NumAtoms());
          f = new_f;
      }
      _manip->_vconf.push_back(f);
  }

  void OBMolManip::DeleteConformers(int start_idx, int end_idx)
  {
    MOLMANIPSANITY(void())
    if (start_idx < 0 || start_idx >= (signed)_manip->_vconf.size() || end_idx < 0 ||
            end_idx >= (signed)_manip->_vconf.size() || start_idx >= end_idx)
      return;
    for (std::vector<double*>::iterator it = _manip->_vconf.begin()+start_idx;
            it!=_manip->_vconf.begin()+end_idx+1; ++it){
        delete [] *it;
    }
    //delete elements between first and last, not including the one pointed to by last
    _manip->_vconf.erase((_manip->_vconf.begin()+start_idx), (_manip->_vconf.begin()+end_idx+1));
  }

  //nothing to be done so far
  OBMolManip::~OBMolManip(){
  }

} // end namespace OpenBabel

//! \file molmanip.cpp
//! \brief Manipulate the geometries of molecules. Implementation of OBMolManip
