/**********************************************************************
Copyright (C) 2011 by Noel O'Boyle

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
#include <openbabel/op.h>
#include <openbabel/depict/depict.h>
#include <openbabel/depict/cairopainter.h>

using namespace std;



namespace OpenBabel
{

class PNG2Format : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  PNG2Format()
  {
    OBConversion::RegisterFormat("png2",this);
  }

  virtual const char* Description() //required
  {
    return
    "PNG2 format\n"
    "2D depiction of a single molecule using Cairo\n"
    ;
  };


  virtual unsigned int Flags()
  {
      return NOTREADABLE | WRITEBINARY | WRITEONEONLY;
  };

  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};
  ////////////////////////////////////////////////////

//Make an instance of the format class
PNG2Format thePNG2Format;

/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////

bool PNG2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
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
      obErrorLog.ThrowError("PNG2Format", "gen2D not found", obError, onceOnly);
      return false;
    }
    if(!pOp->Do(&workingmol))
    {
      obErrorLog.ThrowError("PNG2Format", string(workingmol.GetTitle()) + "- Coordinate generation unsuccessful", obError);
      return false;
    }
  }
  if(!workingmol.Has2D() && workingmol.NumAtoms()>1)
  {
    string mes("Molecule ");
    mes += workingmol.GetTitle();
    mes += " needs 2D coordinates to display in PNG2format";
    obErrorLog.ThrowError("PNG2Format", mes, obError);
    return false;
  }
  const char* pp = pConv->IsOption("p");
  int size  = pp ? atoi(pp) : 300;
  CairoPainter painter;
  OBDepict depictor(&painter);
  depictor.DrawMolecule(&workingmol);
  painter.WriteImage(ofs, size, size);

  return true; //or false to stop converting
}

} //namespace OpenBabel

