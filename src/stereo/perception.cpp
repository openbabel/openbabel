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

    ////////////////////////////////////////////////////////////////////////////
    //
    //  General (0D)
    //
    ////////////////////////////////////////////////////////////////////////////

    /**
     * Perform symmetry analysis.
     */
    void FindSymmetry()
    {
      // FIXME: replace with OBSymmetry
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

      // do quick test to see if there are any possible chiral centers
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

    /**
     * Find all cis/trans double bonds and add a OBCisTransStereo object for each one.
     */
    std::vector<unsigned long> FindCisTrans()
    {
      std::vector<unsigned long> bonds;
          
      //do quick test to see if there are any possible chiral centers
      bool mayHaveCisTransBond = false;
      std::vector<OBBond*>::iterator i;
      for (OBBond *bond = m_mol->BeginBond(i); bond; bond = m_mol->NextBond(i))
        if (bond->GetBO() == 2) {
          mayHaveCisTransBond = true;
          break;
        }

      if (!mayHaveCisTransBond)
        return bonds;

      if (symClasses.size() != m_mol->NumAtoms())
        FindSymmetry();

      bool isCisTrans;
      for (OBBond *bond = m_mol->BeginBond(i); bond; bond = m_mol->NextBond(i)) {
        if (bond->GetBO() == 2) {
          OBAtom *begin = bond->GetBeginAtom();
          OBAtom *end = bond->GetEndAtom();
          if (!begin || !end) 
            continue;

          // two implicit atoms on either begin or end?
          if ((begin->GetValence() == 1) || (end->GetValence() == 1))
            continue;
          
          isCisTrans = true;
          std::vector<OBBond*>::iterator j;
          
          if (begin->GetValence() == 2) {
            // begin atom has two neighbors, the first is the end atom. The second should 
            // be a heavy atom in which case the thirth will be assumed to be implicit
            // (hydrogen, lone pair on N, ...)
            for (OBAtom *nbr = begin->BeginNbrAtom(j); nbr; nbr = begin->NextNbrAtom(j)) {
              if (nbr->GetId() == end->GetId())
                continue;
              // two implicit atoms?
              if (nbr->IsHydrogen())
                isCisTrans = false;
            }
          } else { 
            // valence == 3
            std::vector<unsigned int> tlist;
            
            for (OBAtom *nbr = begin->BeginNbrAtom(j); nbr; nbr = begin->NextNbrAtom(j)) {
              // skip end atom
              if (nbr->GetId() == end->GetId())
                continue;
              // do we already have an atom with this symmetry class?
              if (tlist.size()) {
                // compare second with first
                if (symClasses[nbr->GetIdx()-1] == tlist.at(0))
                  isCisTrans = false;
                break;
              }
                
              // save first summetry class
              tlist.push_back(symClasses[nbr->GetIdx()-1]);
            }
          }

          if (end->GetValence() == 2) {
            // begin atom has two neighbors, the first is the end atom. The second should 
            // be a heavy atom in which case the thirth will be assumed to be implicit
            // (hydrogen, lone pair on N, ...)
            for (OBAtom *nbr = end->BeginNbrAtom(j); nbr; nbr = end->NextNbrAtom(j)) {
              if (nbr->GetId() == begin->GetId())
                continue;
              // two implicit atoms?
              if (nbr->IsHydrogen())
                isCisTrans = false;
            }
          } else { 
            // valence == 3
            std::vector<unsigned int> tlist;
            
            for (OBAtom *nbr = end->BeginNbrAtom(j); nbr; nbr = end->NextNbrAtom(j)) {
              // skip end atom
              if (nbr->GetId() == begin->GetId())
                continue;
              // do we already have an atom with this symmetry class?
              if (tlist.size()) {
                // compare second with first
                if (symClasses[nbr->GetIdx()-1] == tlist.at(0))
                  isCisTrans = false;
                break;
              }
                
              // save first summetry class
              tlist.push_back(symClasses[nbr->GetIdx()-1]);
            }
          }

          if (isCisTrans)
            bonds.push_back(bond->GetId());
        }
      }

      return bonds;
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //  From3D
    //
    ////////////////////////////////////////////////////////////////////////////

    //! Calculate the "sign of a volume" given by a set of 4 coordinates
    double VolumeSign(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d)
    {
      vector3 A, B, C;
      A = b - a;
      B = c - a;
      C = d - a;
      matrix3x3 m(A, B, C);
      return m.determinant();
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
           
        // Create a vector with the coordinates of the neighbor atoms
        std::vector<vector3> nbrCoords;
        OBAtom *from = m_mol->GetAtomById(config.from);
        nbrCoords.push_back(from->GetVector());
        for (OBStereo::RefIter id = config.refs.begin(); id != config.refs.end(); ++id) {
          OBAtom *nbr = m_mol->GetAtomById(*id);
          nbrCoords.push_back(nbr->GetVector());
        }
    
        // Checks for a neighbour having 0 co-ords (added hydrogen etc)
        /* FIXME: needed? if the molecule has 3D coords, additional 
         * hydrogens will get coords using OBAtom::GetNewBondVector
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
        */

        // If we have three heavy atoms we can use the chiral center atom itself for the fourth
        // will always give same sign (for tetrahedron), magnitude will be smaller.
        if ((config.refs.size() == 2) || use_central_atom) {
          nbrCoords.push_back(center->GetVector());
          config.refs.push_back(OBStereo::ImplicitId); // need to add largest number on end to work
        }

        double sign = VolumeSign(nbrCoords[0], nbrCoords[1], nbrCoords[2], nbrCoords[3]);
        if (sign < 0.0)
          config.winding = OBStereo::AntiClockwise;

        OBTetrahedralStereo *th = new OBTetrahedralStereo(m_mol);
        th->SetConfig(config);
        m_mol->SetData(th);
      }
    
    }


    void CisTransFrom3D()
    {
      // find all cis/trans bonds
      std::vector<unsigned long> bonds = FindCisTrans();
    
      std::vector<unsigned long>::iterator i;
      for (i = bonds.begin(); i != bonds.end(); ++i) {
        OBBond *bond = m_mol->GetBondById(*i);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();

        // Create a vector with the coordinates of the neighbor atoms
        std::vector<vector3> bondVecs;
        OBCisTransStereo::Config config;
        // begin
        config.begin = begin->GetId();
        FOR_NBORS_OF_ATOM (nbr, begin) {
          if (nbr->GetId() == end->GetId())
            continue;
          config.refs.push_back(nbr->GetId());
          bondVecs.push_back(nbr->GetVector() - begin->GetVector());
        }
        if (config.refs.size() == 1) {
          config.refs.push_back(OBStereo::ImplicitId);
          vector3 pos;
          m_mol->GetAtomById(config.refs.at(0))->GetNewBondVector(pos, 1.0);
          bondVecs.push_back(pos - begin->GetVector());
        }
        // end
        config.end = end->GetId();
        FOR_NBORS_OF_ATOM (nbr, end) {
          if (nbr->GetId() == begin->GetId())
            continue;
          config.refs.push_back(nbr->GetId());
          bondVecs.push_back(nbr->GetVector() - end->GetVector());
        }
        if (config.refs.size() == 3) {
          config.refs.push_back(OBStereo::ImplicitId);
          vector3 pos;
          m_mol->GetAtomById(config.refs.at(2))->GetNewBondVector(pos, 1.0);
          bondVecs.push_back(pos - end->GetVector());
        }

        // 0      3       0   3
        //  \    /         \ /      angle 0-*-3 & 1-*-2: 60 degrees (cis)
        //   C==C    -->    *       angle 0-*-1 & 2-*-3: 120 degrees (same bond atom)
        //  /    \         / \      angle 0-*-2 & 1-*-3: 180 degrees (trans)
        // 1      2       1   2
        double angle = vectorAngle(bondVecs[0], bondVecs[2]);
        if (IsNear(angle, 60.0, 10.0))
          config.shape = OBStereo::ShapeZ;
        if (IsNear(angle, 180.0, 10.0))
          config.shape = OBStereo::ShapeU;

        OBCisTransStereo *ct = new OBCisTransStereo(m_mol);
        ct->SetConfig(config);
        m_mol->SetData(ct);
      }

    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //  From2D
    //
    //  Reference:
    //  [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the 
    //  Unambiguous Identification of the Stereochemical Characteristics of 
    //  Compounds During Their Registration in Databases. Molecules 2000, 6,
    //  915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
    ////////////////////////////////////////////////////////////////////////////

    //! Calculate the "sign of a triangle" given by a set of 4 coordinates
    double TriangleSign(const vector3 &a, const vector3 &b, const vector3 &c)
    {
      // equation 6 from [1]
      return (a.x() - c.x()) * (b.y() - c.y()) - (a.y() - c.y()) * (b.x() - c.x());
    }

    void TetrahedralFrom2D()
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
        
        // find the hash, wedge and 2 plane atoms
        OBAtom *plane1 = 0;
        OBAtom *plane2 = 0;
        OBAtom *hash = 0;
        OBAtom *wedge = 0;
        FOR_BONDS_OF_ATOM(bond, center) {
          OBAtom *nbr = bond->GetNbrAtom(center);
          // hash bonds
          if (bond->IsHash()) {
            if (bond->GetBeginAtom()->GetId() == center->GetId()) {
              // this is a 'real' hash bond going from center to nbr
              if (hash) {
                // we already have a 'real' hash
                plane1 = plane2 = hash = wedge = 0;
                break;              
              }
              hash = nbr;
            } else {
              // this is an 'inverted' hash bond going from nbr to center
              if (wedge) {
                // we already have a 'real' wedge
                plane1 = plane2 = hash = wedge = 0;
                break;              
              }
              wedge = nbr;
            }
            continue;
          }
          // wedge bonds
          if (bond->IsWedge()) {
            if (bond->GetBeginAtom()->GetId() == center->GetId()) {
              // this is a 'real' wedge bond going from center to nbr
              if (wedge) {
                // we already have a 'real' hash
                plane1 = plane2 = hash = wedge = 0;
                break;              
              }
              wedge = nbr;
            } else {
              // this is an 'inverted' hash bond going from nbr to center
              if (hash) {
                // we already have a 'real' wedge
                plane1 = plane2 = hash = wedge = 0;
                break;              
              }
              hash = nbr;
            }
            continue;
          }
          // plane bonds
          if (!plane1)
            plane1 = nbr;
          else
            plane2 = nbr;
        }
        
        OBTetrahedralStereo::Config config;
        config.center = *i;
        using namespace std;

        // plane1 + plane2, hash, wedge
        if (plane1 && plane2 && hash && wedge) {
          config.from = wedge->GetId();
          config.refs.resize(3);
          config.refs[0] = plane1->GetId();
          config.refs[1] = plane2->GetId();
          config.refs[2] = hash->GetId();
          double sign = TriangleSign(plane1->GetVector(), plane2->GetVector(), hash->GetVector());
          cout << "sign = " << sign << endl;
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else
        // plane1 + plane2, hash
        if (plane1 && plane2 && hash) {
          config.from = OBStereo::ImplicitId;
          config.refs.resize(3);
          config.refs[0] = plane1->GetId();
          config.refs[1] = plane2->GetId();
          config.refs[2] = hash->GetId();
          double sign = TriangleSign(plane1->GetVector(), plane2->GetVector(), hash->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else
        // plane1 + plane2, wedge
        if (plane1 && plane2 && wedge) {
          config.towards = OBStereo::ImplicitId;
          config.view = OBStereo::ViewTowards;
          config.refs.resize(3);
          config.refs[0] = plane1->GetId();
          config.refs[1] = plane2->GetId();
          config.refs[2] = wedge->GetId();
          double sign = TriangleSign(plane1->GetVector(), plane2->GetVector(), wedge->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else {
          std::stringstream errorMsg;
          errorMsg << "Symmetry analysis found atom with id " << center->GetId() << 
                   " to be a tetrahedral atom but the wedge/hash bonds can't be interpreted."
                   << std::endl;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
          continue;
        }

        OBTetrahedralStereo *th = new OBTetrahedralStereo(m_mol);
        th->SetConfig(config);
        m_mol->SetData(th);
      }
    
    }


    void CisTransFrom2D()
    {
      // find all cis/trans bonds
      std::vector<unsigned long> bonds = FindCisTrans();
    
      std::vector<unsigned long>::iterator i;
      for (i = bonds.begin(); i != bonds.end(); ++i) {
        OBBond *bond = m_mol->GetBondById(*i);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();

        // Create a vector with the coordinates of the neighbor atoms
        std::vector<vector3> bondVecs;
        OBCisTransStereo::Config config;
        // begin
        config.begin = begin->GetId();
        FOR_NBORS_OF_ATOM (nbr, begin) {
          if (nbr->GetId() == end->GetId())
            continue;
          config.refs.push_back(nbr->GetId());
          bondVecs.push_back(nbr->GetVector());
        }
        if (config.refs.size() == 1) {
          config.refs.push_back(OBStereo::ImplicitId);
          vector3 pos;
          m_mol->GetAtomById(config.refs.at(0))->GetNewBondVector(pos, 1.0);
          bondVecs.push_back(pos);
        }
        // end
        config.end = end->GetId();
        FOR_NBORS_OF_ATOM (nbr, end) {
          if (nbr->GetId() == begin->GetId())
            continue;
          config.refs.push_back(nbr->GetId());
          bondVecs.push_back(nbr->GetVector());
        }
        if (config.refs.size() == 3) {
          config.refs.push_back(OBStereo::ImplicitId);
          vector3 pos;
          m_mol->GetAtomById(config.refs.at(2))->GetNewBondVector(pos, 1.0);
          bondVecs.push_back(pos);
        }

        // 0      3       
        //  \    /        2 triangles: 0-1-b & 2-3-a
        //   a==b    -->  same sign: U
        //  /    \        opposite sign: Z
        // 1      2       
        double sign1 = TriangleSign(bondVecs[0], bondVecs[1], end->GetVector());
        double sign2 = TriangleSign(bondVecs[2], bondVecs[3], begin->GetVector());
        double sign = sign1 * sign2;

        if (sign < 0.0) // opposite sign
          config.shape = OBStereo::ShapeZ;

        OBCisTransStereo *ct = new OBCisTransStereo(m_mol);
        ct->SetConfig(config);
        m_mol->SetData(ct);
      }

    }


    OBMol *m_mol;
    std::vector<unsigned int> symClasses;
  };

  void StereoFrom3D(OBMol *mol, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;
    mol->SetChiralityPerceived();
      
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom3D", obAuditMsg);

    StereoPerception perceive(mol);
    perceive.TetrahedralFrom3D();
    perceive.CisTransFrom3D();
  }

  void StereoFrom2D(OBMol *mol, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;
    mol->SetChiralityPerceived();
      
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom2D", obAuditMsg);

    StereoPerception perceive(mol);
    perceive.TetrahedralFrom2D();
    perceive.CisTransFrom2D();
  }



}

