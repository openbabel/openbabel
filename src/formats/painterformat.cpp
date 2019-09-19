/**********************************************************************
painterformat.cpp  - Output a set of painter commands for 2D depiction

Copyright (C) 2012 by Noel O'Boyle

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
#include <openbabel/op.h>
#include <openbabel/depict/depict.h>
#include <openbabel/depict/commandpainter.h>

using namespace std;



namespace OpenBabel
{

class PainterFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  PainterFormat()
  {
    OBConversion::RegisterFormat("paint",this);
  }

  virtual const char* Description() //required
  {
    return
    "Painter format\n"
    "Commands used to generate a 2D depiction of a molecule\n\n"

    "This is a utility format that is useful if you want to\n"
    "generate a depiction of a molecule yourself, for example\n"
    "by drawing on a Graphics2D canvas in Java. The format\n"
    "writes out a list of drawing commands as shown\n"
    "in the following example::\n\n"

    "  obabel -:CC(=O)Cl -opaint\n\n"

"  NewCanvas 149.3 140.0\n"
"  SetPenColor 0.0 0.0 0.0 1.0 (rgba)\n"
"  DrawLine 109.3 100.0 to 74.6 80.0\n"
"  SetPenColor 0.0 0.0 0.0 1.0 (rgba)\n"
"  DrawLine 71.6 80.0 to 71.6 53.0\n"
"  DrawLine 77.6 80.0 to 77.6 53.0\n"
"  SetPenColor 0.0 0.0 0.0 1.0 (rgba)\n"
"  DrawLine 74.6 80.0 to 51.3 93.5\n"
"  SetPenColor 0.4 0.4 0.4 1.0 (rgba)\n"
"  SetPenColor 0.4 0.4 0.4 1.0 (rgba)\n"
"  SetPenColor 1.0 0.1 0.1 1.0 (rgba)\n"
"  SetFontSize 16\n"
"  SetFontSize 16\n"
"  SetFontSize 16\n"
"  DrawText 74.6 40.0 \"O\"\n"
"  SetPenColor 0.1 0.9 0.1 1.0 (rgba)\n"
"  SetFontSize 16\n"
"  SetFontSize 16\n"
"  SetFontSize 16\n"
"  SetFontSize 16\n"
"  DrawText 40.0 100.0 \"Cl\"\n\n"

"Note that the origin is considered to be in the top left corner.\n\n"

"The following image was drawn using the information\n"
"in this format as described at\n"
"http://baoilleach.blogspot.co.uk/2012/04/painting-molecules-your-way-introducing.html:\n\n"

".. image:: ../_static/bananamol.png\n\n"

    "Write Options e.g. -xM\n"
    " M Do not include a margin around the depiction\n\n"

    ;
  };


  virtual unsigned int Flags()
  {
      return NOTREADABLE;
  };

  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};
  ////////////////////////////////////////////////////

//Make an instance of the format class
PainterFormat thePainterFormat;

/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////

bool PainterFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream& ofs = *pConv->GetOutStream();

  OBMol workingmol(*pmol); // Copy the molecule

  //*** Coordinate generation ***
  //Generate coordinates only if no existing 2D coordinates
  if(!workingmol.Has2D(true))
  {
    OBOp* pOp = OBOp::FindType("gen2D");
    if(!pOp)
    {
      obErrorLog.ThrowError("PainterFormat", "gen2D not found", obError, onceOnly);
      return false;
    }
    if(!pOp->Do(&workingmol))
    {
      obErrorLog.ThrowError("PainterFormat", string(workingmol.GetTitle()) + "- Coordinate generation unsuccessful", obError);
      return false;
    }
  }
  if(!workingmol.Has2D() && workingmol.NumAtoms()>1)
  {
    string mes("Molecule ");
    mes += workingmol.GetTitle();
    mes += " needs 2D coordinates to display in PNG2format";
    obErrorLog.ThrowError("PainterFormat", mes, obError);
    return false;
  }

  CommandPainter painter(*pConv->GetOutStream());
  OBDepict depictor(&painter);
  if(pConv->IsOption("M"))
    depictor.SetOption(OBDepict::noMargin);
  depictor.DrawMolecule(&workingmol);

  return true; //or false to stop converting
}

} //namespace OpenBabel

