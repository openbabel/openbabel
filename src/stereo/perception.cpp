/**********************************************************************
  perception.cpp - Stereochemistry perception

  Copyright (C) 2009 by Tim Vandermeersch
 
  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/canon.h>
#include <openbabel/oberror.h>

namespace OpenBabel {
 
  struct StereoPerception
  {
    StereoPerception(OBMol *mol) : m_mol(mol) 
    {
    }

    /**
     * Perform symmetry analysis.
     */
    void FindSymmetry()
    {
      // FIXME: replace with OBSym
      std::vector<unsigned int> canlbls;
      OBBitVec allbits(m_mol->NumAtoms());
      FOR_ATOMS_OF_MOL(a, m_mol)
        allbits.SetBitOn(a->GetIdx());
      CanonicalLabels(m_mol, allbits, symClasses, canlbls);
    }

    /**
     * Find all tetrahedral atoms using the symmetry classes.
     */
    std::vector<unsigned long> FindTetrahedral()
    {
      std::vector<unsigned long> centers;

      //do quick test to see if there are any possible chiral centers
      bool mayHaveChiralCenter = false;
      OBAtom *atom, *nbr;
      std::vector<OBAtom*>::iterator i;
      for (atom = m_mol->BeginAtom(i); atom; atom = m_mol->NextAtom(i))
        if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3) {
          mayHaveChiralCenter = true;
          break;
        }

      if (!mayHaveChiralCenter)
        return centers;

      if (symClasses.size() != m_mol->NumAtoms())
        FindSymmetry();
      std::vector<unsigned int> tlist;
      std::vector<unsigned int>::iterator k;

      bool ischiral;
      for (atom = m_mol->BeginAtom(i); atom; atom = m_mol->NextAtom(i)) {
        if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3) {
          tlist.clear();
          ischiral = true;

          std::vector<OBBond*>::iterator j;
          for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j)) {
            for (k = tlist.begin(); k != tlist.end(); ++k)
              if (symClasses[nbr->GetIdx()-1] == *k)
                ischiral = false;

            if (ischiral)
              tlist.push_back(symClasses[nbr->GetIdx()-1]);
            else
              break;
          }

          if (ischiral)
            centers.push_back(atom->GetId());
        }
      }

      return centers;
    }

    //! Calculate a signed volume given a set of 4 coordinates
    double signed_volume(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d)
    {
      vector3 A,B,C;
      A = b-a;
      B = c-a;
      C = d-a;
      matrix3x3 m(A,B,C);
      return(m.determinant());
    }
 
    void TetrahedralFrom3D()
    {
      // find all tetrahedral centers
      std::vector<unsigned long> centers = FindTetrahedral();
      
      std::vector<unsigned long>::iterator i;
      for (i = centers.begin(); i != centers.end(); ++i) {
        OBAtom *center = m_mol->GetAtomById(*i);
 
        // make sure we have at least 3 heavy atom neighbors
        if (center->GetHvyValence() < 3) {
          std::stringstream errorMsg;
          errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of " 
                   << center->GetHvyValence() << std::endl;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
          continue;
        }
        
        OBTetrahedralStereo::Config config;
        config.center = *i;
        FOR_NBORS_OF_ATOM(nbr, center) {
          if (config.from == OBStereo::NoId)
            config.from = nbr->GetId();
          else
            config.refs.push_back(nbr->GetId());
        }

        bool use_central_atom = false;
           
        // Create a vector with neighbor atom ids
        /*
        OBStereo::Refs nbrIds;
        FOR_NBORS_OF_ATOM(nbr, center) {
          nbrIds.push_back(nbr->GetId());
        }

        // sort the neighbor atoms to insure a consistent ordering
        sort(nbrIds.begin(), nbrIds.end());
        */
        // Create a vector with the coordinates of the neighbor atoms
        std::vector<vector3> nbrCoords;
        OBAtom *from = m_mol->GetAtomById(config.from);
        nbrCoords.push_back(from->GetVector());
        for (OBStereo::RefIter id = config.refs.begin(); id != config.refs.end(); ++id) {
          OBAtom *nbr = m_mol->GetAtomById(*id);
          nbrCoords.push_back(nbr->GetVector());
        }
    
        // Checks for a neighbour having 0 co-ords (added hydrogen etc)
        for (std::vector<vector3>::iterator coord = nbrCoords.begin(); coord != nbrCoords.end(); ++coord) { 
          // are the coordinates zero to 6 or more significant figures
          if (coord->IsApprox(VZero, 1.0e-6)) {
            if (!use_central_atom) {
              use_central_atom = true;
            } else {
              obErrorLog.ThrowError(__FUNCTION__, 
                  "More than 2 neighbours have 0 co-ords when attempting 3D chiral calculation", obInfo);
            }
          }
        }

        // If we have three heavy atoms we can use the chiral center atom itself for the fourth
        // will always give same sign (for tetrahedron), magnitude will be smaller.
        if ((config.refs.size() == 2) || use_central_atom) {
          nbrCoords.push_back(center->GetVector());
          config.refs.push_back(OBStereo::ImplicitId); // need to add largest number on end to work
        }

        double V = signed_volume(nbrCoords[0], nbrCoords[1], nbrCoords[2], nbrCoords[3]);
        if (V < 0.0)
          config.winding = OBStereo::AntiClockwise;

        OBTetrahedralStereo *th = new OBTetrahedralStereo(m_mol);
        th->SetConfig(config);
        m_mol->SetData(th);
      }
    
    }

    void CisTransFrom3D()
    {
    
    }

    OBMol *m_mol;
    std::vector<unsigned int> symClasses;
  };



  void StereoFrom3D(OBMol *mol)
  {
    if (mol->HasChiralityPerceived())
      return;
    mol->SetChiralityPerceived();
      
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom3D", obAuditMsg);

    StereoPerception perceive(mol);
    perceive.TetrahedralFrom3D();
    perceive.CisTransFrom3D();
   
  
  
  }


}

