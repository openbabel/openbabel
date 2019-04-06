/**********************************************************************
Copyright (C) 2014 by Mathias Laurin

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

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <cstdlib>

namespace OpenBabel {

class AoforceFormat : public OBMoleculeFormat {
  public:
    // Register this format type ID
    AoforceFormat() { OBConversion::RegisterFormat("aoforce", this); }

    virtual const char* Description() {  // required
      return
          "Turbomole AOFORCE output format\n"
          "Read vibrational frequencies and intensities\n";
    }

    virtual const char* SpecificationURL() {
      return "http://www.turbomole-gmbh.com/manuals/";
    }

    virtual unsigned int Flags() { return READONEONLY | NOTWRITABLE; }

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
};

// Instantiate
AoforceFormat theAoforceFormat;

bool AoforceFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv) {
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if (pmol == NULL) return false;

  std::istream &ifs = *pConv->GetInStream();
  std::string line;
  OBMol &mol = *pmol;

  std::vector<double> Frequencies;
  std::vector<double> Intensities;
  std::vector< std::vector<vector3> > Lx;
  mol.BeginModify();
  while (std::getline(ifs, line)) {
    std::vector<std::string> vs;
    if (line.find("atomic coordinates") != std::string::npos) {
      // The next lines contain 8 fields:
      // x y z atom shells charge pseudo isotop
      while (std::getline(ifs, line) && line.length()) {
        tokenize(vs, line);
        OBAtom *atom = mol.NewAtom();
        vector3 coords(atof(vs[0].c_str()),
                       atof(vs[1].c_str()),
                       atof(vs[2].c_str()));
        coords *= 0.529177249;  // Bohr to Angstrom
        atom->SetVector(coords);
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[3].c_str()));
        atom->SetPartialCharge(atof(vs[5].c_str()));
        atom->SetIsotope(atoi(vs[7].c_str()));
      }
    } else if (line.find("   mode   ") != std::string::npos) {
      // Normal modes and vibrational frequencies
      std::getline(ifs, line);  // empty line
      std::getline(ifs, line);  // frequency
      tokenize(vs, line);
      // for each frequency
      for (std::vector<std::string>::const_iterator
          iter = vs.begin() + 1; iter != vs.end(); ++iter) {
        Frequencies.push_back(atof(iter->c_str()));
      }
      std::getline(ifs, line);  // empty line
      std::getline(ifs, line);  // symmetry
      std::getline(ifs, line);  // empty line
      std::getline(ifs, line);  // IR
      std::getline(ifs, line);  // |dDIP/dQ|   (a.u.)
      std::getline(ifs, line);  // intensity (km/mol)
      tokenize(vs, line);
      // for each intensity
      for (std::vector<std::string>::const_iterator
          iter = vs.begin() + 2; iter != vs.end(); ++iter) {
        Intensities.push_back(atof(iter->c_str()));
      }
      std::getline(ifs, line);  // intensity (%)
      std::getline(ifs, line);  // empty line
      std::getline(ifs, line);  // RAMAN
      std::getline(ifs, line);  // empty line
      Lx.resize(Frequencies.size());
      // normal modes for each atom
      for (unsigned int atom = 0; atom != mol.NumAtoms(); ++atom) {
        std::vector<double> xs;
        std::getline(ifs, line);  // idx, element, "x", [list]
        tokenize(vs, line);
        for (std::vector<std::string>::const_iterator
            iter = vs.begin() + 3; iter != vs.end(); ++iter) {
          xs.push_back(atof(iter->c_str()));
        }
        std::vector<double> ys;
        std::getline(ifs, line);  // "y", [list]
        tokenize(vs, line);
        for (std::vector<std::string>::const_iterator
            iter = vs.begin() + 1; iter != vs.end(); ++iter) {
          ys.push_back(atof(iter->c_str()));
        }
        std::vector<double> zs;
        std::getline(ifs, line);  // "z", [list]
        tokenize(vs, line);
        for (std::vector<std::string>::const_iterator
            iter = vs.begin() + 1; iter != vs.end(); ++iter) {
          zs.push_back(atof(iter->c_str()));
        }
        // for each new frequency
        std::vector< std::vector<vector3> >::iterator
        lxIter = Lx.end() - xs.size();
        std::vector<double>::const_iterator
        xIter = xs.begin(),
        yIter = ys.begin(),
        zIter = zs.begin();
        for (; xIter != xs.end() && yIter != ys.end() && zIter != zs.end()
            && lxIter != Lx.end(); ++xIter, ++yIter, ++zIter, ++lxIter) {
          // push normal modes
          lxIter->push_back(vector3(*xIter, *yIter, *zIter));
        }
      }
    }
  }
  OBVibrationData *vd = new OBVibrationData;
  vd->SetData(Lx, Frequencies, Intensities);
  mol.SetData(vd);
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  return true;
}
}  // namespace OpenBabel
