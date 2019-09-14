/**********************************************************************
Copyright (C) 2002-2006 by Geoffrey Hutchison
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
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>


using namespace std;
namespace OpenBabel
{

  class ZINDOFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ZINDOFormat()
    {
      OBConversion::RegisterFormat("zin",this);
    }

    virtual const char* Description() //required
    {
      return
        "ZINDO input format\n"
        "The input format for the semiempirical quantum-mechanics program ZINDO.\n"
        "Write Options e.g. -xc\n"
        "  c  Write an input file for the CNDO/INDO program.\n\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE | WRITEONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  ZINDOFormat theZINDOFormat;

  ////////////////////////////////////////////////////////////////

  bool ZINDOFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    int orbitals, valenceE = 0;
    vector<OBAtom*>::iterator i;
    OBAtom *atom;
    bool charged = (pmol->GetTotalCharge() == 0);

    bool cndoStyle = pConv->IsOption("c",OBConversion::OUTOPTIONS);

    if (cndoStyle) {
      ofs << mol.GetTitle() << '\n';
      ofs << "HAMILT= INDO CHARGE=" << pmol->GetTotalCharge() << " RESTART= auto\n";
      ofs << "STOP=CI MAX_CI=50 CI_DUMP=25\n";
      ofs << "MAX_ITS=500\n";
      int spinState = pmol->GetTotalSpinMultiplicity();
      if (spinState != 1) // not a singlet, so give a reasonable default
        ofs << "HIGHSPIN= " << spinState - 1 << "  " << spinState - 1 << '\n';
      ofs << endl;
    }
    else { // old-style ZINDO, so we have to calculate the # of valence electrons
      for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i))
        {
          switch (atom->GetAtomicNum())
            {
            case 1:
              valenceE += 1;
              break;
            case 5:
            case 13:
              valenceE += 3;
              break;
            case 6:
            case 14:
              valenceE += 4;
              break;
            case 7:
            case 15:
            case 33:
              valenceE += 5;
              break;
            case 8:
            case 16:
            case 34:
              valenceE += 6;
              break;
            case 9:
            case 17:
            case 35:
              valenceE += 7;
              break;
            default:
              break;
            }
        }

      orbitals = (valenceE - pmol->GetTotalCharge()) / 2;
      valenceE -= pmol->GetTotalCharge();

      ofs << " $TITLEI" << '\n';
      ofs << '\n';
      ofs << "   " << mol.GetTitle() << '\n';
      ofs << '\n';
      ofs << " $END" << '\n';
      ofs << '\n';
      ofs << " $CONTRL" << '\n';
      ofs << '\n';
      if (charged)
        ofs << " SCFTYP       ROHF   RUNTYP       CI   ENTTYP     COORD" << '\n';
      else
        ofs << " SCFTYP        RHF   RUNTYP       CI   ENTTYP     COORD" << '\n';
      ofs << " UNITS        ANGS   INTTYP        1   IAPX           3" << '\n';
      if (charged)
        {
          ofs << '\n';
          ofs << " NOP = 1 " << '\n';
          ofs << " NDT = 1 " << '\n';
          snprintf(buffer, BUFF_SIZE, " FOP(1) =% 4d% 10.6f",
                   valenceE - 1, 1.0);
          ofs << buffer << '\n';
        }

      snprintf(buffer, BUFF_SIZE, " NAT          %4d   NEL        %4d   MULT           1",
               mol.NumAtoms(),
               valenceE);
      ofs << buffer << '\n';
      ofs << " IPRINT         -1   ITMAX       100" << '\n';
      ofs << '\n';
      ofs << "! ***** BASIS SET AND C. I. SIZE INFORMATION ***** " << '\n';
      ofs << '\n';

      snprintf(buffer, BUFF_SIZE, " DYNAL(1) =     0%5d%5d    0    0 1200%5d",
               mol.NumAtoms() - mol.NumHvyAtoms(),
               mol.NumHvyAtoms(),
               orbitals + 25);
      ofs << buffer << '\n';

      ofs << '\n';
      ofs << " INTFA(1) =   1.000000 1.267000  0.680000  1.000000  1.000000 " << '\n';
      ofs << '\n';
      ofs << "! ***** OUTPUT FILE NAME ***** " << '\n';
      ofs << '\n';
      ofs << "   ONAME =  zindo " << '\n';
      ofs << '\n';
      ofs << " $END" << '\n';
      ofs << '\n';
      ofs << " $DATAIN " << '\n';
      ofs << '\n';
    } // end of old-style ZINDO header


    for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i))
      {
        snprintf(buffer, BUFF_SIZE, "% 10.6f% 10.6f% 10.6f%5d",
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ(),
                 atom->GetAtomicNum());
        ofs << buffer << '\n';
      }

    ofs << '\n' << '\n';

    if (!cndoStyle) { // old-style ZINDO footer for CI specification
      ofs << '\n';
      ofs << " $END " << '\n';
      ofs << '\n';
      ofs << " $CIINPU" << '\n';
      ofs << '\n';
      ofs << "! ***** C. I. SPECIFICATION *****" << '\n';
      ofs << '\n';
      ofs << "    2    1   25    1    0    0    0    1   10    1   10" << '\n';
      ofs << "  -60000.0 0.0000000" << '\n';
      ofs << '\n';

      if (charged)
        snprintf(buffer, BUFF_SIZE, "%5d%5d%5d%5d", 1, orbitals, orbitals, orbitals + 1);
      else
        snprintf(buffer, BUFF_SIZE, "%5d%5d%5d", 1, orbitals, orbitals);
      ofs << buffer << '\n';
      if (charged)
        snprintf(buffer, BUFF_SIZE, "%5d%5d%5d%5d%5d",
                 21, (orbitals - 8), orbitals + 1, orbitals + 1, orbitals + 11);
      else
        snprintf(buffer, BUFF_SIZE, "%5d%5d%5d%5d%5d",
                 21, (orbitals - 9), orbitals, orbitals + 1, orbitals + 10);
      ofs << buffer << "\n\n";
      ofs << " $END \n";
    }

    return(true);
  }

} //namespace OpenBabel
