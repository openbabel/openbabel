/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2014 by Dagmar Lenk

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
#include <openbabel/generic.h>
#include <openbabel/obiter.h>

#include <sstream>
#include <iomanip>
#include <cstdlib>

#define notFound string::npos

using namespace std;
namespace OpenBabel
{

  class EXYZFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    EXYZFormat()
    {
      OBConversion::RegisterFormat("exyz", this, "chemical/x-xyz");
    }

    virtual const char* Description() //required
    {
      return
        "Extended XYZ cartesian coordinates format\n"
"A format used by ORCA-AICCM\n"
"The \"EXYZ\" chemical file format is an extended version of the standard\n"
"\"XYZ\" chemical file format with additional keywords and informations about\n"
"the unit cell and virtual atoms.\n\n"

"* Line one of the file contains the number of atoms in the file.\n"
"* Line two of the file contains a title, comment, filename and/or the\n"
"  following keywords: ``%PBC`` or ``%VIRTUAL`` \n\n"

"Any remaining lines are parsed for atom information until a blank line. These\n"
"lines start with the element symbol, followed by X, Y, and Z coordinates in\n"
"angstroms separated by whitespace and - if ``%VIRTUAL`` is specified - the \n"
"optional word ``VIRTUAL`` to mark virtual atoms. If ``%PBC`` is specified\n"
"a second block will be present containing the 3 vectors for the unit cell\n"
"in angstrom and the offset as shown in the example below::\n\n"

"  4\n"
"  %PBC\n"
"     C        0.00000        1.40272        0.00000\n"
"     H        0.00000        2.49029        0.00000\n"
"     C       -1.21479        0.70136        0.00000\n"
"     H       -2.15666        1.24515        0.00000\n"
"\n"
"  Vector1    2.445200    0.000000    0.000000\n"
"  Vector2    0.000000    1.000000    0.000000\n"
"  Vector3    0.000000    0.000000    1.000000\n"
"  Offset     0.000000    0.000000    0.000000\n\n"

"On **output**, the first line written is the number of atoms in the molecule.\n"
"Line two is the title of the molecule or the filename if no title is defined.\n"
"Remaining lines define the atoms in the file. The first column is the atomic\n"
"symbol (right-aligned on the third character), followed by the XYZ coordinates\n"
"in \"15.5\" format separated by an addition whitespace, in angstroms. This means\n"
"that all coordinates are printed with five decimal places.\n\n"
"The next block starts with a blank line to separate the coordinates from the\n"
"unit cell vectors followed by the vectors of the unit cell marked with the\n"
"keywords Vector1/2/3. The vectors themselves are written in the same format\n"
"as the atom coordinates. The last line contains the keyword Offset and the\n"
"offset of the unit cell. The unit is always angstrom.\n"

        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.org/wiki/EXYZ";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-exyz"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  EXYZFormat theEXYZFormat;

  /////////////////////////////////////////////////////////////////
  bool EXYZFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();
    char buffer[BUFF_SIZE];

    stringstream errorMsg;
    bool unitCell,virtualAtoms;


    unsigned int natoms;	// [ejk] assumed natoms could not be -ve

    if (!ifs)
      return false; // we're attempting to read past the end of the file

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an E-XYZ file: Cannot read the first line.", obWarning);
        return(false);
      }

    if (sscanf(buffer, "%d", &natoms) == 0 || !natoms)
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an E-XYZ file: The first line must contain the number of atoms.", obWarning);
        return(false);
      }

    mol.ReserveAtoms(natoms);

    // The next line contains a title string for the molecule. Use this
    // as the title for the molecule if the line is not
    // empty. Otherwise, use the title given by the calling function.
    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an EXYZ file: Could not read the second line (title/comments/keywords).", obWarning);
        return(false);
      }
    string readTitle(buffer);
    string::size_type location = readTitle.find("Energy");
    if (location != notFound)
      readTitle.erase(location);
    Trim(readTitle);

    location = readTitle.find_first_not_of(" \t\n\r");
    if (readTitle.find_first_not_of(" \t\n\r") != notFound)
      mol.SetTitle(readTitle);
    else
      mol.SetTitle(title);

    string readKeywords(buffer);
    location = readKeywords.find("%PBC");   // file contains unitcell information behind the coords block
    if (readKeywords.find("%PBC") != notFound)
        unitCell = true;
    else
        unitCell = false;

    location = readKeywords.find("%VIRTUAL");   // file contains information about virtual atoms in the column next to the x,y,z values
    if (readKeywords.find("%VIRTUAL") != notFound)
        virtualAtoms = true;
    else
        virtualAtoms = false;

    mol.BeginModify();

    // The next lines contain four items each, separated by white
    // spaces: the atom type, and the coordinates of the atom
    vector<string> vs;
    for (unsigned int i = 1; i <= natoms; i ++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          {
            errorMsg << "Problems reading an XYZ file: "
                     << "Could not read line #" << i+2 << ", file error." << endl
                     << " According to line one, there should be " << natoms
                     << " atoms, and therefore " << natoms+2 << " lines in the file.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
        tokenize(vs,buffer);
        if (vs.size() < 4) // ignore extra columns which some applications add
          {
            errorMsg << "Problems reading an XYZ file: "
                     << "Could not read line #" << i+2 << "." << endl
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "However, OpenBabel found " << vs.size() << " items.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }

        // Atom Type: get the atomic number from the element table, using
        // the first entry in the currently read line. If the entry makes
        // sense, set the atomic number and leave the atomic type open
        // (the type is then later faulted in when atom->GetType() is
        // called). If the entry does not make sense to use, set the atom
        // type manually, assuming that the author of the xyz-file had
        // something "special" in mind.
        OBAtom *atom  = mol.NewAtom();

        int atomicNum = OBElements::GetAtomicNum(vs[0].c_str());
        //set atomic number, or '0' if the atom type is not recognized
        if (atomicNum == 0) {
          // Sometimes people call this an XYZ file, but it's actually Unichem
          // i.e., the first column is the atomic number, not a symbol
          // so we'll first check if we can convert this to an element number
          atomicNum = atoi(vs[0].c_str());
        }

        atom->SetAtomicNum(atomicNum);
        if (atomicNum == 0) // still strange, try using an atom type
          atom->SetType(vs[0]);

        // Read the atom coordinates
        char *endptr;
        double x = strtod((char*)vs[1].c_str(),&endptr);
        if (endptr == (char*)vs[1].c_str())
          {
            errorMsg << "Problems reading an XYZ file: "
                     << "Could not read line #" << i+2 << "." << endl
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "OpenBabel could not interpret item #1 as a number.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
        double y = strtod((char*)vs[2].c_str(),&endptr);
        if (endptr == (char*)vs[2].c_str())
          {
            errorMsg << "Problems reading an XYZ file: "
                     << "Could not read line #" << i+2 << "." << endl
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "OpenBabel could not interpret item #2 as a number.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
        double z = strtod((char*)vs[3].c_str(),&endptr);
        if (endptr == (char*)vs[3].c_str())
          {
            errorMsg << "Problems reading an XYZ file: "
                     << "Could not read line #" << i+2 << "." << endl
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "OpenBabel could not interpret item #3 as a number.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
        atom->SetVector(x,y,z); //set coordinates

        // OK, sometimes there's sym x y z charge -- accepted by Jmol
        if (vs.size() > 5) {
          string::size_type decimal = vs[4].find('.');
          if (decimal !=string::npos) { // period found
            double charge = strtod((char*)vs[4].c_str(),&endptr);
            if (endptr != (char*)vs[4].c_str())
              atom->SetPartialCharge(charge);
          }
        } // attempt to parse charges
      }
    if (!ifs.getline(buffer,BUFF_SIZE)) {           // skip empty line
        errorMsg << "Problems reading an EXYZ file: "
                 << "Expecting unitcell information because of keyword %PBC " << endl
                 << " but nothing more found ";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
        return(false);
    }
    if (unitCell) {
        OBUnitCell *unitCellInfo = new OBUnitCell;
        matrix3x3 unitCellMatrix;
        for (unsigned int i = 1; i <= 3; i ++) {
            if (!ifs.getline(buffer,BUFF_SIZE)) {
                errorMsg << "Problems reading an EXYZ file: "
                         << "Expecting unitcell information because of keyword %PBC " << endl
                         << " but nothing found ";
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
                return(false);
            }
            tokenize(vs,buffer);
            if (vs.size() != 4) { // assume that unit cell vectors are given
                errorMsg << "Problems reading an EXYZ file: "
                         << "OpenBabel found the line '" << buffer << "'" << endl
                         << "According to the specifications (keyword %PBC), this line should contain exactly 4 entries, separated by white space." << endl
                         << "However, OpenBabel found " << vs.size() << " items.";
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
                return(false);
            }
            // Read the vectors of the unit cell
            char *endptr;
            double x = strtod((char*)vs[1].c_str(),&endptr);
            if (endptr == (char*)vs[1].c_str())
            {
                errorMsg << "Problems reading an EXYZ file: "
                         << "Could not read line #" << i+2 << "." << endl
                         << "OpenBabel found the line '" << buffer << "'" << endl
                         << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                         << "OpenBabel could not interpret item #1 as a number.";

                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
                return(false);
            }
            double y = strtod((char*)vs[2].c_str(),&endptr);
            if (endptr == (char*)vs[2].c_str())
            {
                errorMsg << "Problems reading an EXYZ file: "
                         << "Could not read line #" << i+2 << "." << endl
                         << "OpenBabel found the line '" << buffer << "'" << endl
                         << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                         << "OpenBabel could not interpret item #2 as a number.";

                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
                return(false);
            }
            double z = strtod((char*)vs[3].c_str(),&endptr);
            if (endptr == (char*)vs[3].c_str())
            {
                errorMsg << "Problems reading an EXYZ file: "
                         << "Could not read line #" << i+2 << "." << endl
                         << "OpenBabel found the line '" << buffer << "'" << endl
                         << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                         << "OpenBabel could not interpret item #3 as a number.";

                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
                return(false);
            }
            unitCellMatrix.SetRow(i-1, vector3(x,y,z));
        }
        unitCellInfo->SetData(unitCellMatrix);

        if (!ifs.getline(buffer,BUFF_SIZE)) {
            errorMsg << "Problems reading an EXYZ file: "
                     << "Expecting unitcell information because of keyword %PBC " << endl
                     << " looking for center information but nothing found ";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
        }
        tokenize(vs,buffer);
        if (vs.size() != 4) { // assume that unit cell vectors are given
            errorMsg << "Problems reading an EXYZ file: "
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications (keyword %PBC), this line should contain exactly 4 entries, separated by white space." << endl
                     << "However, OpenBabel found " << vs.size() << " items.";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
        }
        // Read the origin (center) vector of the unit cell
        char *endptr;
        double x = strtod((char*)vs[1].c_str(),&endptr);
        if (endptr == (char*)vs[1].c_str())
        {
            errorMsg << "Problems reading an EXYZ file: "
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "OpenBabel could not interpret item #1 as a number.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
        }
        double y = strtod((char*)vs[2].c_str(),&endptr);
        if (endptr == (char*)vs[2].c_str())
        {
            errorMsg << "Problems reading an EXYZ file: "
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "OpenBabel could not interpret item #2 as a number.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
        }
        double z = strtod((char*)vs[3].c_str(),&endptr);
        if (endptr == (char*)vs[3].c_str())
        {
            errorMsg << "Problems reading an EXYZ file: "
                     << "OpenBabel found the line '" << buffer << "'" << endl
                     << "According to the specifications, this line should contain exactly 4 entries, separated by white space." << endl
                     << "OpenBabel could not interpret item #3 as a number.";

            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
        }
        unitCellInfo->SetOffset(vector3(x,y,z));
        mol.SetData(unitCellInfo);
    }


    // clean out any remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool EXYZFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    string keyword = "%PBC";

    if (string(mol.GetTitle()).find("%PBC") != notFound) keyword = ""; // title contains already right keyword

    ofs << mol.NumAtoms() << endl;
    if (fabs(mol.GetEnergy()) > 1.0e-3) // nonzero energy field
        ofs << mol.GetTitle() << " " << mol.GetEnergy() << " "<<  keyword << endl;
    else
        ofs << mol.GetTitle() << " " <<  keyword << endl;

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        ofs << setw(4) << right
            << OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum())
            << setw(15) << setprecision(5) << fixed << showpoint
            << right << atom->GetX() << " " << setw(15) << atom->GetY() << " "
            << setw(15) << atom->GetZ() << endl;
    }
    ofs << endl;
    //
    // now write unitcell information
    //
    if (!mol.HasData(OBGenericDataType::UnitCell)){
        ofs << right << "Vector1"
            << setw(15) << setprecision(5) << fixed << showpoint
            << right << 1. << " " << setw(15) << 0. << " "<< setw(15) << 0. << endl;
        ofs << "Vector2" << setw(15)
            << right << 0. << " " << setw(15) << 1. << " " << setw(15) << 0. << endl;
        ofs << "Vector3" << setw(15)
            << right << 0. << " " << setw(15) << 0. << " "<< setw(15) << 1. << endl;
        ofs << "Offset " << setw(15)
            << right << 0. << " " << setw(15) << 0. << " " << setw(15) << 0. << endl;
    } else {
        OBUnitCell *uC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
        matrix3x3 unitCellMatrix = uC->GetCellMatrix();
        vector3 offset = uC->GetOffset();
        ofs << right << "Vector1"
            << setw(15) << setprecision(5) << fixed << showpoint
            << right << unitCellMatrix.GetRow(0).GetX() << " " << setw(15) << unitCellMatrix.GetRow(0).GetY() << " "
            << setw(15) << unitCellMatrix.GetRow(0).GetZ() << endl;
        ofs << "Vector2" << setw(15)
            << right << unitCellMatrix.GetRow(1).GetX() << " " << setw(15) << unitCellMatrix.GetRow(1).GetY() << " "
            << setw(15) << unitCellMatrix.GetRow(1).GetZ() << endl;
        ofs << "Vector3" << setw(15)
            << right << unitCellMatrix.GetRow(2).GetX() << " " << setw(15) << unitCellMatrix.GetRow(2).GetY() << " "
            << setw(15) << unitCellMatrix.GetRow(2).GetZ() << endl;
        ofs << "Offset " << setw(15)
            << right << offset.GetX() << " " << setw(15) << offset.GetY() << " "
            << setw(15) << offset.GetZ() << endl;
    }
    return(true);
  }

} //namespace OpenBabel
