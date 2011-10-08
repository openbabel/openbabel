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
    "2D depiction of a single molecule as a .png file\n\n"

    "The PNG2 format is used 'behind the scenes' by the :ref:`PNG format<PNG_2D_depiction>`\n"
    "if generating image files, and the best way to use it is\n"
    "actually through the PNG format. While it possible to generate\n"
    "a :file:`.png` file directly using the PNG2 format as follows...::\n\n"
    "  obabel -:\"CC(=O)Cl\" -opng2 -O mymol.png\n\n"
    "...it is much better to generate it using the PNG format\n"
    "as this allows you to embed a chemical structure in the\n"
    ":file:`.png` file header which you can later extract::\n\n"
    "  $ obabel -:\"CC(=O)Cl\" -O mymol.png -xO smi\n"
    "  $ obabel mymol.png -osmi\n"
    "  CC(=O)Cl\n\n"

    "The PNG2 format uses the Cairo library to generate the\n"
    ":file:`.png` files.\n"
    "If Cairo was not found when Open Babel was compiled, then\n"
    "this format will be unavailable. However, it will still be possible\n"
    "to use the PNG format to read :file:`.png` files if they contain\n"
    "embedded information.\n\n"

    ".. seealso::\n\n"

    "    :ref:`PNG_2D_depiction`\n\n"

    "Write Options e.g. -xp 500\n"
    " p <pixels> image size, default 300\n\n"

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

