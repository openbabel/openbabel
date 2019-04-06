/**********************************************************************
Copyright (C) 2005-2006 by Geoffrey R. Hutchison
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
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include <cstdlib>

#include <openbabel/math/matrix3x3.h>

using namespace std;
namespace OpenBabel
{

  class FreeFormFractionalFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    FreeFormFractionalFormat()
    {
      OBConversion::RegisterFormat("fract",this);
    }

    virtual const char* Description() //required
    {
      return
"Free Form Fractional format\n"
"General purpose crystallographic format\n"
"The \"free-form\" fractional format attempts to allow for input from a\n"
"range of fractional / crystallography file formats. As such, it has only\n"
"a few restrictions on input:\n\n"

"- Line one of the file contains a title or comment.\n"
"- Line two of the file contains the unit cell parameters separated by\n"
"  whitespace and/or commas (i.e. \"a b c alpha beta gamma\").\n"
"- Any remaining lines are parsed for atom information. Lines start with\n"
"  the element symbol, followed by fractional X, Y, and Z coordinates\n"
"  (in angstroms) separated by whitespace.\n\n"

"Any numeric input (i.e., unit cell parameters, XYZ coordinates) can include\n"
"designations of errors, although this is currently ignored. For example::\n\n"

"  C 1.00067(3) 2.75(2) 3.0678(12)\n\n"

"will be parsed as::\n\n"

"  C 1.00067 2.75 3.0678\n\n"

"When used as an **output** format, The first line written is the title of the\n"
"molecule or the filename if no title is defined. If a molecule has a defined\n"
"unit cell, then the second line will be formatted as::\n\n"

"  a b c alpha beta gamma\n\n"

"where a, b, c are the unit cell vector lengths, and alpha, beta, and gamma are\n"
"the angles between them. These numbers are formatted as \"10.5\", which means that\n"
"5 decimal places will be output for all numbers. In the case where no unit cell\n"
"is defined for the molecule, the vector lengths will be defined as 1.0, and the\n"
"angles to 90.0 degrees.\n\n"

"Remaining lines define the atoms in the file. The first column is the atomic\n"
"symbol, followed by the XYZ coordinates in 10.5 format (in angstroms).\n\n"

"Here is an example file::\n\n"

" ZnO test file\n"
" 3.14 3.24 5.18 90.0 90.0 120.0\n"
" O 0.66667  0.33333  0.3750\n"
" O 0.33333  0.66667  0.8750\n"
" Zn 0.66667  0.33333  0.0000\n"
" Zn 0.33333  0.66667  0.5000\n\n"

"Read Options e.g. -as\n"
"  s  Output single bonds only\n"
"  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.org/wiki/Free_Form_Fractional";}; //optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  FreeFormFractionalFormat theFreeFormFractionalFormat;

/*  const char * TrimErrors(const std::string data)
  {
    string temp = data;
    size_t stdErr = temp.rfind("(");

    if(stdErr!=string::npos)
      temp.erase(stdErr);

    return temp.c_str();
  }
*/
  /////////////////////////////////////////////////////////////////
  bool FreeFormFractionalFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading a free form fractional file: Could not read the first line (title/comments).", obWarning);
        return(false);
      }
    if (strlen(buffer) != 0)
      mol.SetTitle(buffer);
    else
      mol.SetTitle(title);

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading a free form fractional file: Could not read the second line (unit cell parameters a b c alpha beta gamma).",
                              obWarning);
        return(false);
      }
    vector<string> vs;
    tokenize(vs,buffer," \n\t,");
    if (vs.size() != 6)
      return(false);

    //parse cell values
    double A, B, C, Alpha, Beta, Gamma;
    string temp; // used to trim ending (xx) data from strings

    A = atof(vs[0].c_str());
    B = atof(vs[1].c_str());
    C = atof(vs[2].c_str());
    Alpha = atof(vs[3].c_str());
    Beta  = atof(vs[4].c_str());
    Gamma = atof(vs[5].c_str());
    OBUnitCell *uc = new OBUnitCell;
    uc->SetOrigin(fileformatInput);
    uc->SetData(A, B, C, Alpha, Beta, Gamma);
    mol.SetData(uc);

    mol.BeginModify();

    string str;
    double x,y,z;
    vector3 v;
    int atomicNum;
    OBAtom *atom;

    while(ifs.getline(buffer,BUFF_SIZE))
      {
        if (strlen(buffer) == 0 || *buffer == 0x0D) //incl Windows kludge
          // blank line -- consider it the end of this molecule
          break;

        tokenize(vs,buffer);
        if (vs.size() != 4)
          return(false);

        atom = mol.NewAtom();

        // check to see if first column is number or element symbol
        // (PCModel has files of the form X Y Z symbol)
        atomicNum = OBElements::GetAtomicNum(vs[0].c_str());
        if (atomicNum == 0 && (isdigit(vs[0][0]) || ispunct(vs[0][0])))
          {
            x = atof(vs[0].c_str());
            y = atof(vs[1].c_str());
            z = atof(vs[2].c_str());
            atomicNum = OBElements::GetAtomicNum(vs[3].c_str());
          }
        else
          {
            x = atof(vs[1].c_str());
            y = atof(vs[2].c_str());
            z = atof(vs[3].c_str());
          }
        v.Set(x, y, z);
        v = uc->FractionalToCartesian(v);
        atom->SetVector(v);

        atom->SetAtomicNum(atomicNum);
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
    if (!pConv->IsOption("s",OBConversion::INOPTIONS)
        && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool FreeFormFractionalFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    OBUnitCell *uc = NULL;

    ofs << mol.GetTitle() << endl;

    if (!mol.HasData(OBGenericDataType::UnitCell))
      ofs << "   1.00000   1.00000   1.00000  90.00000  90.00000  90.00000\n";
    else
      {
        uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
        snprintf(buffer, BUFF_SIZE,
                 "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f",
                 uc->GetA(), uc->GetB(), uc->GetC(),
                 uc->GetAlpha() , uc->GetBeta(), uc->GetGamma());
        ofs << buffer << "\n";
      }

    vector3 v;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        v = atom->GetVector();
        if (uc != NULL)
          v = uc->CartesianToFractional(v);

        snprintf(buffer, BUFF_SIZE, "%s %10.5f%10.5f%10.5f",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 v.x(),
                 v.y(),
                 v.z());
        ofs << buffer << endl;
      }
    ofs << endl; // add a blank line between molecules
    return(true);
  }

} //namespace OpenBabel
