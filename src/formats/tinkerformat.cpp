/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

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
#include <openbabel/elements.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/generic.h>

#include <openbabel/forcefield.h>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class TinkerFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    TinkerFormat()
    {
      OBConversion::RegisterFormat("txyz",this);
    }

    virtual const char* Description() //required
    {
      return
        "Tinker XYZ format\n"
        "The cartesian XYZ file format used by the molecular mechanics package TINKER.\n"
        "By default, the MM2 atom types are used for writing files but MM3 atom types\n"
        "are provided as an option. Another option provides the ability to take the\n"
        "atom type from the atom class (e.g. as used in SMILES, or set via the API).\n\n"

        "Read Options e.g. -as\n"
        "  s  Generate single bonds only\n\n"
        "Write Options e.g. -xm\n"
        "  m  Write an input file for the CNDO/INDO program.\n"
        "  c  Write atom types using custom atom classes, if available\n"
        "  3  Write atom types for the MM3 forcefield.\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://dasher.wustl.edu/tinker/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  TinkerFormat theTinkerFormat;

  // Sets the MM3 atom type based on the MM2 atom type
  int SetMM3Type(OBAtom *atom);

  /////////////////////////////////////////////////////////////////

  bool TinkerFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int natoms;
    char buffer[BUFF_SIZE];
    vector<string> vs;

    ifs.getline(buffer, BUFF_SIZE);
    tokenize(vs,buffer);
    if (vs.size() < 2)
      return false;
    natoms = atoi(vs[0].c_str());

    // title is 2nd token (usually add tokens for the atom types)
    mol.SetTitle(vs[1]);

    mol.ReserveAtoms(natoms);
    mol.BeginModify();

    string str;
    double x,y,z;
    OBAtom *atom;

    for (int i = 1; i <= natoms; ++i)
    {
        if (!ifs.getline(buffer,BUFF_SIZE))
            return(false);
        tokenize(vs,buffer);
        // e.g. "2  C      2.476285    0.121331   -0.001070     2     1     3    14"
        if (vs.size() < 5)
            return(false);

        atom = mol.NewAtom();
        x = atof((char*)vs[2].c_str());
        y = atof((char*)vs[3].c_str());
        z = atof((char*)vs[4].c_str());
        atom->SetVector(x,y,z); //set coordinates

        //set atomic number
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));

        // add bonding
        if (vs.size() > 6)
          for (unsigned int j = 6; j < vs.size(); ++j)
            mol.AddBond(mol.NumAtoms(), atoi((char *)vs[j].c_str()), 1); // we don't know the bond order

    }
    if (!pConv->IsOption("s",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // clean out remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    mol.EndModify();
    mol.SetTitle(title);
    return(true);
  }

  bool TinkerFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;
    bool mm2Types = false;
    bool mmffTypes = pConv->IsOption("m",OBConversion::OUTOPTIONS) != NULL;
    bool mm3Types = pConv->IsOption("3",OBConversion::OUTOPTIONS) != NULL;
    bool classTypes = pConv->IsOption("c", OBConversion::OUTOPTIONS) != NULL;

    unsigned int i;
    char buffer[BUFF_SIZE];
    OBBond *bond;
    vector<OBBond*>::iterator j;

    // Before we try output of MMFF94 atom types, check if it works
    OBForceField *ff = OpenBabel::OBForceField::FindForceField("MMFF94");
    if (mmffTypes && ff && ff->Setup(mol))
      mmffTypes = ff->GetAtomTypes(mol);
    else
      mmffTypes = false; // either the force field isn't available, or it doesn't work

    if (!mmffTypes && !mm3Types && !classTypes) {
      snprintf(buffer, BUFF_SIZE, "%6d %-20s   MM2 parameters\n",mol.NumAtoms(),mol.GetTitle());
      mm2Types = true;
    }
    else if (mm3Types)
      snprintf(buffer, BUFF_SIZE, "%6d %-20s   MM3 parameters\n",mol.NumAtoms(),mol.GetTitle());
    else if (classTypes)
      snprintf(buffer, BUFF_SIZE, "%6d %-20s   Custom parameters\n",mol.NumAtoms(),mol.GetTitle());
    else
      snprintf(buffer, BUFF_SIZE, "%6d %-20s   MMFF94 parameters\n",mol.NumAtoms(),mol.GetTitle());
    ofs << buffer;

    ttab.SetFromType("INT");

    OBAtom *atom;
    string str,str1;
    int atomType;
    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        str = atom->GetType();
        atomType = 0; // Something is very wrong if this doesn't get set below

        if (mm2Types) {
          ttab.SetToType("MM2");
          ttab.Translate(str1,str);
          atomType = atoi((char*)str1.c_str());
        }
        if (mmffTypes) {
          // Override the MM2 typing
          OBPairData *type = (OpenBabel::OBPairData*)atom->GetData("FFAtomType");
          if (type) {
            str1 = type->GetValue().c_str();
            atomType = atoi((char*)str1.c_str());
          }
        }
        if (mm3Types) {
          // convert to integer for MM3 typing
          atomType = SetMM3Type(atom);
        }
        if (classTypes) {
          // Atom classes are set by the user, so use those
          OBGenericData *data = atom->GetData("Atom Class");
          if (data) {
            OBPairInteger* acdata = dynamic_cast<OBPairInteger*>(data); // Could replace with C-style cast if willing to live dangerously
            if (acdata) {
              int ac = acdata->GetGenericValue();
              if (ac >= 0)
                atomType = ac;
            }
          }
        }

        snprintf(buffer, BUFF_SIZE, "%6d %2s  %12.6f%12.6f%12.6f %5d",
                 i,
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ(),
                 atomType);
        ofs << buffer;

        for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
          {
            snprintf(buffer, BUFF_SIZE, "%6d", (bond->GetNbrAtom(atom))->GetIdx());
            ofs << buffer;
          }

        ofs << endl;
      }

    return(true);
  }

  int SetMM3Type(OBAtom *atom)
  {
    OBAtom *b; // neighbor
    OBBondIterator i, j;
    int countNeighborO, countNeighborS, countNeighborN, countNeighborC;
    countNeighborO = countNeighborS = countNeighborN = countNeighborC = 0;

    // The MM2 typing isn't very good, so we do this ourselves for the most common atom types
    switch (atom->GetAtomicNum()) {
    case 1: // Hydrogen
      b = atom->BeginNbrAtom(j);
      if (b->IsCarboxylOxygen())
        return 24;
      if (b->GetAtomicNum() == OBElements::Sulfur)
        return 44;
      if (b->GetAtomicNum() == OBElements::Nitrogen) {
        if (b->IsAmideNitrogen())
          return 28;
        if (b->GetExplicitDegree() > 3)
          return 48;// ammonium
        return 23; // default amine/imine
      }
      if (b->GetAtomicNum() == OBElements::Carbon && b->GetHyb() == 1)
        return 124; // acetylene

      if (b->GetAtomicNum() == OBElements::Oxygen) {
        if (b->HasAlphaBetaUnsat())
          return 73; // includes non-enol/phenol, but has the right spirit
        return 21; // default alcohol
      }

      return 5; // default H
      break;

    case 2: // Helium
      return 51; break;
    case 3: // Li
      return 163; break;

    case 5: // B
      if (atom->GetExplicitDegree() >= 4)
        return 27; // tetrahedral
      return 26; break;

    case 6: // C
      if (atom->IsInRingSize(3)) { // cyclopropane / cyclopropene
        if (atom->GetHyb() == 3)
          return 22;
        if (atom->GetHyb() == 2) {
          if (atom->CountFreeOxygens() == 1) // propanone
            return 67;
          return 38; // propane
        }
      }
      if (atom->IsInRingSize(4)) { // cyclobutane or cyclobutene
        if (atom->GetHyb() == 3)
          return 56;
        if (atom->GetHyb() == 2) {
          if (atom->CountFreeOxygens() == 1) // butanone
            return 58;
          return 57; // regular cyclobutane
        }
      }

      if (atom->CountBondsOfOrder(2) == 2) { // allene
        if (atom->CountFreeOxygens() == 1) // ketene
          return 106;
        return 68;
      }

      if (atom->GetFormalCharge() == +1)
        return 30;
      if (atom->GetSpinMultiplicity() == 2)
        return 29;

      if (atom->GetHyb() == 3)
        return 1;
      else if (atom->GetHyb() == 2) {
        if (atom->CountFreeOxygens() >= 1)
          return 3;
        return 2;
      }
      else if (atom->GetHyb() == 1)
        return 4;
      break;

    case 7: // N
      // TODO
      if (atom->IsAmideNitrogen())
        return 151;
      if (atom->IsAromatic()) {
        if (atom->GetFormalCharge() == 1)
          return 111;
        if (atom->IsInRingSize(5)) // pyrrole
          return 40;
        if (atom->IsInRingSize(6)) // pyridine
          return 37;
      }

      if (atom->CountFreeOxygens() == 2) // nitro
        return 46;

      if (atom->GetHyb() == 3) {
        if (atom->GetExplicitDegree() > 3)
          return 39; // ammonium
        return 8;
      }
      else if (atom->GetHyb() == 2)
        return 9;
      else if (atom->GetHyb() == 1)
        return 10;
      break;

    case 8: // O
      //TODO
      if (atom->IsPhosphateOxygen())
        return 159;
      if (atom->IsCarboxylOxygen())
        return 75;
      if (atom->IsInRingSize(3))
        return 49; // epoxy

      b = atom->BeginNbrAtom(j);
      if (atom->HasBondOfOrder(2) && b->GetAtomicNum() == OBElements::Carbon) { // O=C
        return 7;
      }

      if (atom->IsAromatic())
        return 41; // furan
      return 6;
      break;

    case 9: // F
      return 11; break;
    case 10: // Ne
      return 52; break;
    case 12: // Mg
      return 59; break;
    case 14: // Si
      return 19; break;

    case 15: // P
      if (atom->CountFreeOxygens() > 0)
        return 153; // phosphate
      if (atom->GetExplicitValence() > 3)
        return 60; // phosphorus V
      return 25; break;

    case 16: // S
      if (atom->IsAromatic())
        return 42; // thiophene
      if (atom->GetFormalCharge() == 1)
        return 16; // sulfonium

    // look at the neighbors
    for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
      {
        switch (b->GetAtomicNum()) {
        case 6:
          if (b->GetHyb() == 2) // S=C
            countNeighborC++; break;
        case 7:
          countNeighborN++; break;
        case 8:
          if (b->GetHvyDegree() == 1)
            countNeighborO++;
          break;
        case 16:
          countNeighborS++; break;
        default:
          continue;
        }
      }

      if (countNeighborO == 1)
        return 17; // sulfoxide
      if (countNeighborO >= 2) {
        if (countNeighborN)
          return 154; // sulfonamide
        return 18; // sulfone or sulfate
      }
      if (countNeighborC)
        return 74; // S=C
      if (countNeighborS == 1)
        return 104; // S-S disulfide
      else if (countNeighborS > 1)
        return 105;

      return 15; break;

    case 17: // Cl
      return 12; break;
    case 18: // Ar
      return 153; break;
    case 20: // Ca
      return 125; break;
    case 26: // Fe
      if (atom->GetFormalCharge() == 2)
        return 61;
      return 62; break;
    case 27: // Co
      if (atom->GetFormalCharge() == 2)
        return 65;
      return 66; break;
    case 28: // Ni
      if (atom->GetFormalCharge() == 2)
        return 63;
      return 64; break;

    case 32: // Ge
      return 31; break;
    case 34: // Se
      return 34; break;
    case 35: // Br
      return 13; break;
    case 36: // Kr
      return 54; break;
    case 38: // Sr
      return 126; break;
    case 50: // Sn
      return 32; break;
    case 52: // Te
      return 35; break;
    case 53: // I
      return 14; break;
    case 54: // Xe
      return 55; break;

    case 56: // Ba
      return 127; break;
    case 57:
      return 128; break;
    case 58:
      return 129; break;
    case 59:
      return 130; break;
    case 60:
      return 131; break;
    case 61:
      return 132; break;
    case 62:
      return 133; break;
    case 63:
      return 134; break;
    case 64:
      return 135; break;
    case 65:
      return 136; break;
    case 66:
      return 137; break;
    case 67:
      return 138; break;
    case 68:
      return 139; break;
    case 69:
      return 140; break;
    case 70:
      return 141; break;
    case 71:
      return 142; break;

    case 82: // Pb
      return 33; break;


    default:
      break;
    }
    return 0;
  }

} //namespace OpenBabel
