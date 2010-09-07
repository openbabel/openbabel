/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2009 by David C. Lonie for PWscf

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

#define RYDBERG_TO_KCAL_PER_MOL 313.755026
#define RYDBERG_TO_ELECTRON_VOLT 13.60569193
#define BOHR_TO_ANGSTROM .529177

using namespace std;
namespace OpenBabel {
  class PWscfFormat : public OBMoleculeFormat
  {
  public:

    PWscfFormat()
    {
      OBConversion::RegisterFormat("pwscf",this);
    }

    virtual const char* Description()
    {
      return "PWscf format\n"
             "The format used by PWscf, part of Quantum Espresso.\n\n";
    };

    virtual const char* SpecificationURL(){return "http://www.quantum-espresso.org";};

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
  PWscfFormat thePWscfFormat;

  /////////////////////////////////////////////////////////////////

  bool PWscfFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE], tag[BUFF_SIZE];
    double x,y,z;
    double alat = 0;
    vector<string> vs;
    matrix3x3 ortho;
    int atomicNum;
    OBUnitCell *cell = new OBUnitCell();
    pmol->BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE)) {

      if (strstr(buffer, "lattice parameter (a_0)")) {
        tokenize(vs, buffer);
        alat = atof(vs.at(4).c_str());
      }

      // Unit cell info
      if (strstr(buffer, "CELL_PARAMETERS")) {
        // Discover units
        double conv = 1.0;
        tokenize(vs, buffer);

        if (strstr(vs[1].c_str(), "alat")) {
          conv = alat * BOHR_TO_ANGSTROM;
        }
        // Add others if needed

        double v11, v12, v13,
          v21, v22, v23,
          v31, v32, v33;

        ifs.getline(buffer,BUFF_SIZE); // v1
        tokenize(vs, buffer);
        v11 = atof(vs.at(0).c_str()) * conv;
        v12 = atof(vs.at(1).c_str()) * conv;
        v13 = atof(vs.at(2).c_str()) * conv;

        ifs.getline(buffer,BUFF_SIZE); // v2
        tokenize(vs, buffer);
        v21 = atof(vs.at(0).c_str()) * conv;
        v22 = atof(vs.at(1).c_str()) * conv;
        v23 = atof(vs.at(2).c_str()) * conv;

        ifs.getline(buffer,BUFF_SIZE); // v3
        tokenize(vs, buffer);
        v31 = atof(vs.at(0).c_str()) * conv;
        v32 = atof(vs.at(1).c_str()) * conv;
        v33 = atof(vs.at(2).c_str()) * conv;

        // Build unit cell
        cell->SetData(vector3(v11,v12,v13),
                      vector3(v21,v22,v23),
                      vector3(v31,v32,v33));
      }

      // Unit cell info (for non-variable cell calcs)
      if (strstr(buffer, "crystal axes: (cart. coord. in units of a_0)")) {
        double conv = alat * BOHR_TO_ANGSTROM;
        double v11, v12, v13,
          v21, v22, v23,
          v31, v32, v33;

        ifs.getline(buffer,BUFF_SIZE); // v1
        tokenize(vs, buffer);
        v11 = atof(vs.at(3).c_str()) * conv;
        v12 = atof(vs.at(4).c_str()) * conv;
        v13 = atof(vs.at(5).c_str()) * conv;

        ifs.getline(buffer,BUFF_SIZE); // v2
        tokenize(vs, buffer);
        v21 = atof(vs.at(3).c_str()) * conv;
        v22 = atof(vs.at(4).c_str()) * conv;
        v23 = atof(vs.at(5).c_str()) * conv;

        ifs.getline(buffer,BUFF_SIZE); // v3
        tokenize(vs, buffer);
        v31 = atof(vs.at(3).c_str()) * conv;
        v32 = atof(vs.at(4).c_str()) * conv;
        v33 = atof(vs.at(5).c_str()) * conv;

        // Build unit cell
        cell->SetData(vector3(v11,v12,v13),
                      vector3(v21,v22,v23),
                      vector3(v31,v32,v33));
      }

      // Atoms info
      if (strstr(buffer, "ATOMIC_POSITIONS")) {
        // Clear old atoms from pmol
        vector<OBAtom*> toDelete;
        FOR_ATOMS_OF_MOL(a, *pmol)
          toDelete.push_back(&*a);
        for (int i = 0; i < toDelete.size(); i++)
          pmol->DeleteAtom(toDelete.at(i));


        // Discover units
        matrix3x3 conv (1);
        tokenize(vs, buffer);

        if (strstr(vs[1].c_str(), "alat")) {
          conv *= (alat * BOHR_TO_ANGSTROM);
        }
        else if (strstr(vs[1].c_str(), "crystal")) {
          // Set to the zero matrix and test below.
          conv = matrix3x3 (0.0);
        }
        // Add others if needed

        // Load new atoms from molecule
        ifs.getline(buffer,BUFF_SIZE); // First entry
        tokenize(vs, buffer);
        int size = vs.size();
        while (size == 4) {
          atomicNum = etab.GetAtomicNum(vs[0].c_str());
          x = atof((char*)vs[1].c_str());
          y = atof((char*)vs[2].c_str());
          z = atof((char*)vs[3].c_str());
          // Add atom
          OBAtom *atom = pmol->NewAtom();
          atom->SetAtomicNum(atomicNum);
          vector3 coords (x,y,z);
          if (conv.determinant() == 0.0) { // Fractional coords
            atom->SetVector(cell->FractionalToCartesian(coords));
          }
          else {
            atom->SetVector(conv * coords);
          }

          // Reset vars
          ifs.getline(buffer,BUFF_SIZE); // First entry
          tokenize(vs, buffer);
          size = vs.size();
        }
      }

      // Free energy
      if (strstr(buffer, "Final energy =")) {
        tokenize(vs, buffer);
        pmol->SetEnergy(atof(vs[3].c_str()) * RYDBERG_TO_KCAL_PER_MOL);
      }

      // H - PV = U energy
      if (strstr(buffer, "!    total energy              =")) {
        tokenize(vs, buffer);
        pmol->SetEnergy(atof(vs[4].c_str()) * RYDBERG_TO_KCAL_PER_MOL);
      }

      // Enthalphy
      if (strstr(buffer, "Final enthalpy =")) {
        float en, pv;
        tokenize(vs, buffer);

        OBPairData *enthalpy = new OBPairData();
        enthalpy->SetAttribute("Enthalpy (kcal/mol)");
        en = static_cast<float> (atof(vs.at(3).c_str()) * RYDBERG_TO_KCAL_PER_MOL);
        snprintf(tag, BUFF_SIZE, "%f", en);
        enthalpy->SetValue(tag);
        pmol->SetData(enthalpy);

        OBPairData *enthalpy_pv = new OBPairData();
        enthalpy_pv->SetAttribute("Enthalpy PV term (kcal/mol)");
        pv = static_cast<float> (en - pmol->GetEnergy());
        snprintf(tag, BUFF_SIZE, "%f", pv);
        enthalpy_pv->SetValue(tag);
        pmol->SetData(enthalpy_pv);

      }
    }

    // set final unit cell
    pmol->SetData(cell);

    pmol->EndModify();

    return true;
  }

} //namespace OpenBabel
