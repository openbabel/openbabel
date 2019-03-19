/**********************************************************************
Copyright (C) 2003 by Michael Banck <mbanck@gmx.net>
Some portions Copyright (C) 2004 by Chris Morley

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

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>

#include <stdlib.h>

using namespace std;
namespace OpenBabel
{

  class CHTFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CHTFormat()
    {
      OBConversion::RegisterFormat("cht",this);
    }

    virtual const char* Description() //required
    {
      return
        "Chemtool format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://ruby.chemie.uni-freiburg.de/~martin/chemtool/chemtool.html";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  CHTFormat theCHTFormat;

  /////////////////////////////////////////////////////////////////

  bool CHTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    int w, h, x, y; 	// to calculate the geometry
    int bondtype;		// type of bond
    int conv_factor = 50;	// please adjust
    int natoms = 0;		// number of additional (non-carbon) atoms
    OBAtom *atom, *atom1, *atom2;
    OBBond *bond;

    ofs << "Chemtool Version 1.4" << endl;

    // get the geometry
    w = 0;
    h = 0;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i))
      {
        x = (int)(atom->GetX()) * conv_factor;
        y = (int)(atom->GetY()) * conv_factor;
        if (x > w)
          w = x;
        if (y > h)
          h = y;
        if (atom->GetAtomicNum() != 6)
          natoms++;
      }
    ofs << "geometry " << w * 1.1 << " " << h * 1.1 << endl;

    // write out bonds
    ofs << "bonds "<< mol.NumBonds() << endl;
    vector<OBBond*>::iterator j;
    for(bond = mol.BeginBond(j); bond; bond = mol.NextBond(j))
      {
        bondtype = 0;
        atom1 = bond->GetBeginAtom();
        atom2 = bond->GetEndAtom();
        if (bond->GetBondOrder() == 2)
          bondtype = 1;
        if (bond->GetBondOrder() == 3)
          bondtype = 3;
        // @todo: use flag-info, too
        snprintf(buffer, BUFF_SIZE, "%d\t%d\t%d\t%d\t%1d",
                 (int)floor(atom1->GetX() * conv_factor + 0.5),
                 (int)floor(atom1->GetY() * conv_factor + 0.5),
                 (int)floor(atom2->GetX() * conv_factor + 0.5),
                 (int)floor(atom2->GetY() * conv_factor + 0.5),
                 bondtype);
        ofs << buffer << endl;
      }

    // start over, write additional atoms
    ofs << "atoms " << natoms << endl;
    for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i))
      {
        // Carbon does not need to be treated
        if (atom->GetAtomicNum() != 6)
          {
            snprintf(buffer, BUFF_SIZE, "%d\t%d\t%s\t%d",
                     (int)floor(atom->GetX() * conv_factor + 0.5),
                     (int)floor(atom->GetY() * conv_factor + 0.5),
                     OBElements::GetSymbol(atom->GetAtomicNum()),
                     -1 // assume centered Text
                     );
            ofs << buffer << endl;
          }
      }

    // We don't have any splines to write
    ofs << "splines 0" << endl;

    return true;
  }

} //namespace OpenBabel
