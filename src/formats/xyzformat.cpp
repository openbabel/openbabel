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
#include <openbabel/obiter.h>

#include <sstream>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class XYZFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    XYZFormat()
    {
      OBConversion::RegisterFormat("xyz", this, "chemical/x-xyz");
    }

    virtual const char* Description() //required
    {
      return
        "XYZ cartesian coordinates format\n"
"A generic coordinate format\n"
"The \"XYZ\" chemical file format is widely supported by many programs, although no\n"
"formal specification has been published. Consequently, Open Babel attempts to be\n"
"extremely flexible in parsing XYZ format files. Similar formats include Tinker XYZ\n"
"and UniChem XYZ which differ slightly in the format of the files. (Notably, UniChem XYZ\n"
"uses the atomic number rather than element symbol for the first column.)\n\n"

"* Line one of the file contains the number of atoms in the file.\n"
"* Line two of the file contains a title, comment, or filename. \n\n"

"Any remaining lines are parsed for atom information. Lines start with the element\n"
"symbol, followed by X, Y, and Z coordinates in angstroms separated by whitespace.\n\n"

"Multiple molecules / frames can be contained within one file.\n\n"

"On **output**, the first line written is the number of atoms in the molecule\n"
"(warning - the number of digits is limited to three for some programs,\n"
"e.g. Maestro). Line two is the title of the molecule or the filename if\n"
"no title is defined. Remaining lines define the atoms in the file. The\n"
"first column is the atomic symbol (right-aligned on the third character),\n"
"followed by the XYZ coordinates in \"10.5\" format, in angstroms. This means\n"
"that all coordinates are printed with five decimal places.\n\n"

"Example::\n\n"

" 12\n"
" benzene example\n"
"   C        0.00000        1.40272        0.00000\n"
"   H        0.00000        2.49029        0.00000\n"
"   C       -1.21479        0.70136        0.00000\n"
"   H       -2.15666        1.24515        0.00000\n"
"   C       -1.21479       -0.70136        0.00000\n"
"   H       -2.15666       -1.24515        0.00000\n"
"   C        0.00000       -1.40272        0.00000\n"
"   H        0.00000       -2.49029        0.00000\n"
"   C        1.21479       -0.70136        0.00000\n"
"   H        2.15666       -1.24515        0.00000\n"
"   C        1.21479        0.70136        0.00000\n"
"   H        2.15666        1.24515        0.00000\n\n"

        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.org/wiki/XYZ";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-xyz"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  XYZFormat theXYZFormat;

  /////////////////////////////////////////////////////////////////
  bool XYZFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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

    unsigned int natoms;	// [ejk] assumed natoms could not be -ve

    if (!ifs)
      return false; // we're attempting to read past the end of the file

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an XYZ file: Cannot read the first line.", obWarning);
        return(false);
      }

    if (sscanf(buffer, "%d", &natoms) == 0 || !natoms)
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an XYZ file: The first line must contain the number of atoms.", obWarning);
        return(false);
      }

    mol.ReserveAtoms(natoms);

    // The next line contains a title string for the molecule. Use this
    // as the title for the molecule if the line is not
    // empty. Otherwise, use the title given by the calling function.
    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an XYZ file: Could not read the second line (title/comments).", obWarning);
        return(false);
      }
    string readTitle(buffer);
    string::size_type location = readTitle.find("Energy");
    if (location != string::npos)
      readTitle.erase(location);
    Trim(readTitle);

    location = readTitle.find_first_not_of(" \t\n\r");
    if (location != string::npos)
      mol.SetTitle(readTitle);
    else
      mol.SetTitle(title);

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

  bool XYZFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    snprintf(buffer, BUFF_SIZE, "%d\n", mol.NumAtoms());
    ofs << buffer;
    if (fabs(mol.GetEnergy()) > 1.0e-3) // nonzero energy field
      snprintf(buffer, BUFF_SIZE, "%s\tEnergy: %15.7f\n",
               mol.GetTitle(), mol.GetEnergy());
    else
      snprintf(buffer, BUFF_SIZE, "%s\n", mol.GetTitle());
    ofs << buffer;

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%-3s%15.5f%15.5f%15.5f\n",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer;
      }

    return(true);
  }

} //namespace OpenBabel
