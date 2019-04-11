/**********************************************************************
Some portions Copyright (C) 2012 by Geoffrey R. Hutchison

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
#include <openbabel/math/spacegroup.h>

#include <sstream>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class POSFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    POSFormat()
    {
      OBConversion::RegisterFormat("pos", this);
    }

    virtual const char* Description() //required
    {
      return
        "POS cartesian coordinates format\n"
"A generic coordinate format\n"
"The \"POS\" file format is a modified version of the \"XYZ\" general format\n\n"

"* Line one of the file contains the number of atoms in the file.\n"
"* Line two of the file contains a title, comment, or filename. \n\n"

"Example::\n\n"

        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://cst-www.nrl.navy.mil/lattice/";}; //optional

    // it's basically the same as XYZ
    virtual const char* GetMIMEType()
    { return "chemical/x-xyz"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return READONEONLY|NOTWRITABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  POSFormat thePOSFormat;

  /////////////////////////////////////////////////////////////////
  bool POSFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* defaultTitle = pConv->GetTitle();
    char buffer[BUFF_SIZE];

    stringstream errorMsg;
    unsigned int natoms;

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
    Trim(readTitle);

    // This should give the location of the spacegroup number
    string::size_type location = readTitle.find("#");
    int spaceGroup = -1;
    if (location != string::npos) {
      // e.g., "#230"
      string spaceGroupNumber = readTitle.substr(location + 1, 4); // +1 to skip #
      string::size_type nonNumber = spaceGroupNumber.find_first_not_of("0123456789");
      if (nonNumber != string::npos)
        spaceGroupNumber.erase(nonNumber);
      // Finally get the space group from the file
      spaceGroup = atoi(spaceGroupNumber.c_str());
    }

    location = readTitle.find_first_not_of(" \t\n\r");
    // Is there non-whitespace
    if (location != string::npos)
      mol.SetTitle(readTitle);
    else
      mol.SetTitle(defaultTitle);

    vector3 v1, v2, v3;
    double x,y,z;
    vector<string> vs;
    OBAtom *atom;
    int atomicNum;
    bool setCellVectors = false;
    // go through remaining lines, particularly looking for cell vectors
    mol.BeginModify();

    while(ifs.peek() != EOF && ifs.good()) {
      ifs.getline(buffer,BUFF_SIZE);
      if (strstr(buffer, "Primitive vectors")) {
        // three lines: a(#) = .. .. ..
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs,buffer);
        if (vs.size() != 5)
          continue;
        v1.SetX(atof(vs[2].c_str()));
        v1.SetY(atof(vs[3].c_str()));
        v1.SetZ(atof(vs[4].c_str()));

        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs,buffer);
        if (vs.size() != 5)
          continue;
        v2.SetX(atof(vs[2].c_str()));
        v2.SetY(atof(vs[3].c_str()));
        v2.SetZ(atof(vs[4].c_str()));

        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs,buffer);
        if (vs.size() != 5)
          continue;
        v3.SetX(atof(vs[2].c_str()));
        v3.SetY(atof(vs[3].c_str()));
        v3.SetZ(atof(vs[4].c_str()));

        setCellVectors = true;
      }
      if (strstr(buffer, "Basis Vectors:")) {
        ifs.getline(buffer,BUFF_SIZE);  // column titles
        ifs.getline(buffer,BUFF_SIZE); // blank
        // real atomic data
        while (ifs.getline(buffer,BUFF_SIZE) && strlen(buffer)) {
          tokenize(vs,buffer);
          if (vs.size() != 7)
            break;

          atom = mol.NewAtom();
          // check to see if first column is number or element symbol
          // (PCModel has files of the form X Y Z symbol)
          atomicNum = OBElements::GetAtomicNum(vs[0].c_str());
          x = atof(vs[4].c_str());
          y = atof(vs[5].c_str());
          z = atof(vs[6].c_str());
          atom->SetVector(x,y,z);
          atom->SetAtomicNum(atomicNum);
        }
      }
    }

    if (setCellVectors) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(v1, v2, v3);
      uc->SetOrigin(fileformatInput);
      mol.SetData(uc);

      // hopefully, we have a space group too
      if (spaceGroup != -1) {
        uc->SetSpaceGroup(spaceGroup);
      }
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    return(true);
  }

} //namespace OpenBabel
