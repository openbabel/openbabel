/**********************************************************************
rotamer.cpp - Handle rotamer list data.
 
Copyright (C) 1998, 1999, 2000-2002 OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2008 by Tim Vandermeersch
 
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

#include <openbabel/babelconfig.h>
#include <openbabel/rotamer.h>

using namespace std;

namespace OpenBabel
{

  /** \class OBRotamerList rotamer.h <openbabel/rotamer.h>

      A high-level class for rotamer / conformer generation, intended mainly
      for use with the related class OBRotorList and the OBRotorRules database

      Rotamers represent conformational isomers formed simply by rotation of
      dihedral angles. On the other hand, conformers may include geometric
      relaxation (i.e., slight modification of bond lengths, bond angles, etc.)

      The following shows an example of generating 2 conformers using different
      rotor states. Similar code could be used for systematic or Monte Carlo
      conformer sampling when combined with energy evaluation (molecular
      mechanics or otherwise).

      \code
      OBRotorList rl; // used to sample all rotatable bonds via the OBRotorRules data
      // If you want to "fix" any particular atoms (i.e., freeze them in space)
      // then set up an OBBitVec of the fixed atoms and call
      // rl.SetFixAtoms(bitvec);
      rl.Setup(mol);

      // How many rotatable bonds are there?
      cerr << " Number of rotors: " << rl.Size() << endl;

      std::vector<int> rotorKey(rl.Size(), 0);

      // each entry represents the configuration of a rotor
      // e.g. indexes into OBRotor::GetTorsionValues() -- the different angles
      //   to sample for a rotamer search
      for (unsigned int i = 0; i < rl.Size(); ++i)
      rotorKey[i] = 0; // could be anything from 0 .. OBRotor->GetTorsionValues().size()
      // -1 is for no rotation

      // The OBRotamerList can generate conformations (i.e., coordinate sets)
      OBRotamerList rotamers;
      rotamers.Setup(mol, rl);

      // keys can be invalid because they don't have the right size
      // or because one of their values is not a valid index for the 
      // matching OBRotor::GetTorsionValues()
      if (!rotamers.AddRotamer(rotorKey))
        cout << "Invalid key!" << endl;
      
      rotorKey[1] = 2; // switch one rotor
      if (!rotamers.AddRotamer(rotorKey))
        cout << "Invalid key!" << endl;

      rotamers.ExpandConformerList(mol, mol.GetConformers());

      // change the molecule conformation
      mol.SetConformer(0); // rotorKey 0, 0, ...
      conv.Write(&mol);

      mol.SetConformer(1); // rotorKey 0, 2, ...
      conv.Write(&mol);

      \endcode

  **/

#if !HAVE_RINT
  inline double rint(double x)
  {
    return ( (x < 0.0) ? ceil(x-0.5) : floor(x+0.5));
  }
#endif

  //////////////////////////////////////////////////////
  //
  // Private member functions
  //
  //////////////////////////////////////////////////////
  
  void OBRotamerList::Setup(OBMol &mol, unsigned int *ref, int nrotors)
  {
    //clear the old stuff out if necessary
    m_torvals.clear();
    m_keys.clear();

    vector<pair<unsigned int*,vector<int> > >::iterator k;
    for (k = m_rotors.begin();k != m_rotors.end();++k)
      delete [] k->first;
    m_rotors.clear();

    //create the new list
    vector<int> children;
    unsigned int *atomlist;
    for (unsigned int i = 0; i < nrotors; ++i) {
      atomlist = new unsigned int [4];
      atomlist[0] = ref[i*4  ];
      atomlist[1] = ref[i*4+1];
      atomlist[2] = ref[i*4+2];
      atomlist[3] = ref[i*4+3];
      mol.FindChildren(children, atomlist[1], atomlist[2]);
      m_rotors.push_back(pair<unsigned int*,vector<int> > (atomlist, children));
    }
  }
 
  void OBRotamerList::GetReferenceArray(unsigned int *ref) const
  {
    vector<pair<unsigned int*,vector<int> > >::const_iterator i;
    unsigned int j;
    for (j = 0, i = m_rotors.begin(); i != m_rotors.end(); ++i) {
      ref[j++] = i->first[0];
      ref[j++] = i->first[1];
      ref[j++] = i->first[2];
      ref[j++] = i->first[3];
    }
  }

  //////////////////////////////////////////////////////
  //
  // Puplic member functions
  //
  //////////////////////////////////////////////////////
 
  OBGenericData* OBRotamerList::Clone(OBBase* newparent) const
  {
    //Since the class contains OBAtom pointers, the new copy use
    //these from the new molecule, newparent
    OBMol* newmol = static_cast<OBMol*>(newparent);
	
    OBRotamerList *new_rml = new OBRotamerList;
    new_rml->_attr = _attr;
    new_rml->_type = _type;

    //Set the base coordinates
    new_rml->SetBaseCoordinates(GetBaseCoordinates(), NumAtoms());

    //Set reference array
    unsigned int *ref = new unsigned int [NumRotors()*4];
    if (ref) {
      GetReferenceArray(ref);
      new_rml->Setup(*newmol, ref, NumRotors()); 
      delete [] ref;
    }

    //Set Rotamers
    new_rml->SetRotamers(m_keys);
    
    return new_rml;
  }

  OBRotamerList::OBRotamerList() : m_numAtoms(0), m_coords(0)
  {
    _type = OBGenericDataType::RotamerList;
    _attr = "RotamerList";
  }
 
  OBRotamerList::~OBRotamerList()
  {
    vector<pair<unsigned int*,vector<int> > >::iterator j;
    for (j = m_rotors.begin();j != m_rotors.end();++j)
      delete [] j->first;

    //Delete the interal base coordinate list
    if (m_coords)
      delete [] m_coords;
  }

  void OBRotamerList::Setup(OBMol &mol, OBRotorList &rl)
  {
    //clear the old stuff out if necessary
    m_torvals.clear();
    m_keys.clear();

    // set the base coordinates
    SetBaseCoordinates(mol.GetCoordinates(), mol.NumAtoms());
    
    // delete the atoms references
    vector<pair<unsigned int*, vector<int> > >::iterator k;
    for (k = m_rotors.begin(); k != m_rotors.end(); ++k)
      delete [] k->first;
    m_rotors.clear();

    //create the new list
    OBRotor *rotor;
    vector<OBRotor*>::iterator i;
    vector<int> children;

    unsigned int ref[4];
    unsigned int *atomlist;
    for (rotor = rl.BeginRotor(i);rotor;rotor = rl.NextRotor(i)) {
      atomlist = new unsigned int [4];
      rotor->GetDihedralAtoms(ref);
      atomlist[0] = ref[0];
      atomlist[1] = ref[1];
      atomlist[2] = ref[2];
      atomlist[3] = ref[3];
      mol.FindChildren(children, ref[1], ref[2]);
      m_rotors.push_back(pair<unsigned int*,vector<int> > (atomlist, children));
      m_torvals.push_back(rotor->GetTorsionValues());
    }

    vector<double>::iterator n;
    vector<vector<double> >::iterator m;
    for (m = m_torvals.begin();m != m_torvals.end();++m)
      for (n = m->begin();n != m->end();++n)
        *n *= RAD_TO_DEG;
  }


  bool OBRotamerList::AddRotamer(double *c)
  {
    if (!c)
      return false;

    unsigned int idx, size;
    double angle;
    const double res = 255.0 / 360.0; // make sure we can fit 360 in a byte
    Eigen::Vector3d v1, v2, v3, v4;

    vector<unsigned char> key(m_rotors.size());

    vector<pair<unsigned int*,vector<int> > >::iterator i;
    for (size = 0, i = m_rotors.begin(); i != m_rotors.end(); ++i, ++size)
      {
        idx = (i->first[0] - 1) * 3;
        v1 = Eigen::Vector3d(c[idx],c[idx+1],c[idx+2]);
        idx = (i->first[1] - 1) * 3;
        v2 = Eigen::Vector3d(c[idx],c[idx+1],c[idx+2]);
        idx = (i->first[2] - 1) * 3;
        v3 = Eigen::Vector3d(c[idx],c[idx+1],c[idx+2]);
        idx = (i->first[3] - 1) * 3;
        v4 = Eigen::Vector3d(c[idx],c[idx+1],c[idx+2]);

        angle = VectorTorsion(v1, v2, v3, v4);
        while (angle < 0.0)
          angle += 360.0;
        while (angle > 360.0)
          angle -= 360.0;
        key[size] = (unsigned char)rint(angle*res);
      }

    m_keys.push_back(key);

    return true;
  }

  bool OBRotamerList::AddRotamer(std::vector<double> &angles)
  {
    double angle;
    const double res = 255.0 / 360.0; // make sure we can fit 360 in a byte
    
    if (angles.size() != m_rotors.size())
      return false; // wrong size 

    vector<unsigned char> rot(m_rotors.size());
    for (unsigned int i = 0; i < m_rotors.size(); ++i) { // for rotors
      angle = angles[i];

      while (angle < 0.0)
        angle += 360.0;
      while (angle > 360.0)
        angle -= 360.0;

      rot[i] = (unsigned char)rint(angle*res);
    }
  
    m_keys.push_back(rot);
    
    return true;
  }
  
  bool OBRotamerList::AddRotamer(std::vector<int> &key)
  {
    double angle;
    const double res = 255.0 / 360.0; // make sure we can fit 360 in a byte
    
    if (key.size() != m_rotors.size())
      return false; // wrong size key

    vector<unsigned char> rot(m_rotors.size());
    for (unsigned int i = 0; i < m_rotors.size(); ++i) { // for rotors
      if (key[i] < 0) {
        unsigned int idx;
        Eigen::Vector3d v1,v2,v3,v4;

        // this rotor should not be rotated, so we calculate the 
        // angle and store it.
        idx = (m_rotors[i].first[0] - 1) * 3;
        v1 = Eigen::Vector3d( m_coords[idx], m_coords[idx+1], m_coords[idx+2] );
        idx = (m_rotors[i].first[1] - 1) * 3;
        v2 = Eigen::Vector3d( m_coords[idx], m_coords[idx+1], m_coords[idx+2] );
        idx = (m_rotors[i].first[2] - 1) * 3;
        v3 = Eigen::Vector3d( m_coords[idx], m_coords[idx+1], m_coords[idx+2] );
        idx = (m_rotors[i].first[3] - 1) * 3;
        v4 = Eigen::Vector3d( m_coords[idx], m_coords[idx+1], m_coords[idx+2] );
        
        angle = VectorTorsion(v1, v2, v3, v4);
      } else if (key[i] >= (int)m_torvals[i].size()) {
        return false; // invalid index
      } else {
        angle = m_torvals[i][key[i]]; // valid key
      }

      while (angle < 0.0)
        angle += 360.0;
      while (angle > 360.0)
        angle -= 360.0;

      rot[i] = (unsigned char)rint(angle*res);
    }

    m_keys.push_back(rot);
    
    return true;
  }

  void OBRotamerList::ExpandConformerList(OBMol &mol, vector<double*> &clist)
  {
    vector<double*> tmpclist = CreateConformerList(mol);

    //transfer the conf list
    vector<double*>::iterator k;
    for (k = clist.begin();k != clist.end();++k)
      delete [] *k;
    clist = tmpclist;
  }

  //! Create a conformer list using the internal base set of coordinates
  vector<double*> OBRotamerList::CreateConformerList(OBMol& mol)
  {
    double angle;
    const double invres = 360.0 / 255.0;
    vector<unsigned char> key;
    vector<double*> tmpclist;
    vector<vector<unsigned char> >::iterator i;

    for (i = m_keys.begin(); i != m_keys.end(); ++i) {
      key = *i;
      double *c = new double [mol.NumAtoms() * 3];
      // copy the base coordinates to c
      memcpy(c, m_coords, sizeof(double) * mol.NumAtoms() * 3);

      // set the torsions based on the keys
      for (unsigned int j = 0; j < m_rotors.size(); ++j) {
        // convert angle from char (0 ... 255) to double (-180.0 ... 180.0)
        angle = invres * ((double) key[j]);
        if (angle > 180.0)
          angle -= 360.0;
           
        SetTorsion(c, m_rotors[j].first, DEG_TO_RAD * angle, m_rotors[j].second);
      }
      
      tmpclist.push_back(c);
    }

    return tmpclist;
  }
  
  //! Change the current coordinate set
  //! \since version 2.2
  void OBRotamerList::SetCurrentCoordinates(OBMol &mol, std::vector<int> &key)
  {
    unsigned int i;
    double angle;
    
    if (key.size() != (m_rotors.size()))
      return; // wrong size key
    
    double *rot = new double [m_rotors.size()];
    
    double *c = mol.GetCoordinates();
    for (i = 0;i < m_rotors.size();++i) {
      if (key[i] == -1) // skip this rotor
        continue;
      else {
    	angle = m_torvals[i][key[i]];
      	
        while (angle < 0.0)
       	  angle += 360.0;
    	while (angle > 360.0)
     	  angle -= 360.0;
        
        SetTorsion(c, m_rotors[i].first, DEG_TO_RAD * angle, m_rotors[i].second);
      } // set an angle
    } // for rotors
  }

  //Copies the coordinates in bc, NOT the pointers, into the object
  void OBRotamerList::SetBaseCoordinates(double* coords, unsigned int N)
  {
    //Clear out old data
    if (m_coords)
      delete [] m_coords;
    
    if (coords) {
      m_numAtoms = N;
      //Copy new data
      m_coords = new double [3*N];
      for (unsigned int j = 0; j < 3*N; ++j)
        m_coords[j] = coords[j];
    } else {
      m_numAtoms = 0;
      m_coords = 0;
    }
  }

} //namespace OpenBabel

//! \file rotamer.cpp
//! \brief Handle rotamer list data.
