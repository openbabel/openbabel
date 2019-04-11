/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2009 by David C. Lonie for CASTEP

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <cstdlib>

#define EV_TO_KCAL_PER_MOL 23.060538
#define GPA_A3_TO_KCAL_PER_MOL 0.14383639

using namespace std;
namespace OpenBabel {
  class CASTEPFormat : public OBMoleculeFormat
  {
  public:

    CASTEPFormat()
    {
      OBConversion::RegisterFormat("castep",this);
    }

    virtual const char* Description()
    {
      return
        "CASTEP format\n"
        "The format used by CASTEP.\n\n";
    };

    virtual const char* SpecificationURL(){return "http://www.castep.org/";};

    /* Flags() can return be any of the following combined by |
       or be omitted if none apply
       NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
       READBINARY  WRITEBINARY  READXML  ZEROATOMSOK */
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    };

    ////////////////////////////////////////////////////
    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    /* Add declarations for any local function or member variables used.
       Generally only a single instance of a format class is used. Keep this in
       mind if you employ member variables. */
  };
  ////////////////////////////////////////////////////

  //Make an instance of the format class
  CASTEPFormat theCASTEPFormat;

  /////////////////////////////////////////////////////////////////

  bool CASTEPFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();

    bool coordsAreFractional = false;
    // H = U + PV breakdown is not given, must calculate it
    bool hasPressureData = false;
    bool hasVolumeData = false;
    bool hasEnthalpyData = false;
    double pressure, volume, enthalpy;
    char buffer[BUFF_SIZE], tag[BUFF_SIZE];
    double x,y,z,a,b,c,alpha,beta,gamma;
    vector<string> vs;
    matrix3x3 ortho;
    int atomicNum;
    OBUnitCell *cell = new OBUnitCell();
    pmol->BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE)) {
      // Unit cell info
      if (strstr(buffer, "Lattice parameters(A)       Cell Angles")) {
        ifs.getline(buffer,BUFF_SIZE); // a, alpha
        tokenize(vs, buffer);
        a = atof(vs.at(2).c_str());
        alpha = atof(vs.at(5).c_str());

        ifs.getline(buffer,BUFF_SIZE); // b, beta
        tokenize(vs, buffer);
        b = atof(vs.at(2).c_str());
        beta = atof(vs.at(5).c_str());

        ifs.getline(buffer,BUFF_SIZE); // c, gamma
        tokenize(vs, buffer);
        c = atof(vs.at(2).c_str());
        gamma = atof(vs.at(5).c_str());

        // Build unit cell
        cell->SetData(a, b, c, alpha, beta, gamma);
      }

      // Fractional atomic info
      if (strstr(buffer, "x  Element    Atom        Fractional coordinates of atoms  x")) {
        coordsAreFractional = true;
        // Clear old atoms from pmol
        vector<OBAtom*> toDelete;
        FOR_ATOMS_OF_MOL(a, *pmol)
          toDelete.push_back(&*a);
        for (size_t i = 0; i < toDelete.size(); i++)
          pmol->DeleteAtom(toDelete.at(i));

        // Load new atoms from molecule
        ifs.getline(buffer,BUFF_SIZE); // Title line 2
        ifs.getline(buffer,BUFF_SIZE); // ------

        ifs.getline(buffer,BUFF_SIZE); // First entry
        tokenize(vs, buffer);
        int size = vs.size();
        while (size == 7) {
          atomicNum = OBElements::GetAtomicNum(vs[1].c_str());
          x = atof((char*)vs[3].c_str());
          y = atof((char*)vs[4].c_str());
          z = atof((char*)vs[5].c_str());

          // Add atom
          OBAtom *atom = pmol->NewAtom();
          atom->SetAtomicNum(atomicNum);
          vector3 coords (x,y,z);
          atom->SetVector(coords);

          // Reset vars
          ifs.getline(buffer,BUFF_SIZE); // Next entry
          tokenize(vs, buffer);
          size = vs.size();
        }
      }

      // Final enthalpy
      if (strstr(buffer, "Final Enthalpy")) {
        hasEnthalpyData = true;
        tokenize(vs, buffer);
        enthalpy = atof(vs[4].c_str()) * EV_TO_KCAL_PER_MOL;
      }

      // volume
      if (strstr(buffer, "Current cell volume =")) {
        hasVolumeData = true;
        tokenize(vs, buffer);
        volume = atof(vs[4].c_str());
      }

      // pressure
      if (strstr(buffer, " *  Pressure:")) {
        hasPressureData = true;
        tokenize(vs, buffer);
        pressure = atof(vs[2].c_str());
      }
    }

    // Convert coords to cartesian if needed
    if (coordsAreFractional) {
      FOR_ATOMS_OF_MOL(atom, pmol) {
        atom->SetVector(cell->FractionalToCartesian(cell->WrapFractionalCoordinate(atom->GetVector())));
      }
    }

    // Set energy/enthalpy/pv etc
    if (hasEnthalpyData) {
      if (hasVolumeData && hasPressureData) {
        double energy, pv;
        double enthalpy_ev, pv_ev;
        pv = volume * pressure * GPA_A3_TO_KCAL_PER_MOL;
        energy = enthalpy - pv;

        pv_ev = pv / EV_TO_KCAL_PER_MOL;
        enthalpy_ev = enthalpy / EV_TO_KCAL_PER_MOL;

        OBPairData *opd_enthalpy = new OBPairData();
        OBPairData *opd_enthalpy_pv = new OBPairData();
        OBPairData *opd_enthalpy_ev = new OBPairData();
        OBPairData *opd_enthalpy_pv_ev = new OBPairData();

        opd_enthalpy->SetAttribute("Enthalpy (kcal/mol)");
        opd_enthalpy_pv->SetAttribute("Enthalpy PV term (kcal/mol)");
        opd_enthalpy_ev->SetAttribute("Enthalpy (eV)");
        opd_enthalpy_pv_ev->SetAttribute("Enthalpy PV term (eV)");

        snprintf(tag, BUFF_SIZE, "%f", enthalpy);
        opd_enthalpy->SetValue(tag);
        snprintf(tag, BUFF_SIZE, "%f", pv);
        opd_enthalpy_pv->SetValue(tag);
        snprintf(tag, BUFF_SIZE, "%f", enthalpy_ev);
        opd_enthalpy_ev->SetValue(tag);
        snprintf(tag, BUFF_SIZE, "%f", pv_ev);
        opd_enthalpy_pv_ev->SetValue(tag);

        pmol->SetData(opd_enthalpy);
        pmol->SetData(opd_enthalpy_pv);
        pmol->SetData(opd_enthalpy_ev);
        pmol->SetData(opd_enthalpy_pv_ev);
        pmol->SetEnergy(energy);

      }
      else { // special case: H == U @ P=0
        pmol->SetEnergy(enthalpy);
      }
    }

    // set final unit cell
    pmol->SetData(cell);

    pmol->EndModify();

    return true;
  }

} //namespace OpenBabel
