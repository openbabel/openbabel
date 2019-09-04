/**********************************************************************
asciiformat.cpp  - Output a depiction of a molecule as ASCII text

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
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/op.h>
#include <openbabel/depict/depict.h>
#include <openbabel/depict/asciipainter.h>
#include <cstdlib>

using namespace std;



namespace OpenBabel
{

class ASCIIFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  ASCIIFormat()
  {
    OBConversion::RegisterFormat("ascii",this);
  }

  virtual const char* Description() //required
  {
    return
    "ASCII format\n"
    "2D depiction of a single molecule as ASCII text\n\n"

    "This format generates a 2D depiction of a molecule using only ASCII text\n"
    "suitable for a command-line console, or a text file. For example::\n\n"
    "  obabel -:c1ccccc1C(=O)Cl -oascii -xh 20\n"
"  \n"
"         __\n"
"      __/__\\_\n"
"    _/__/    \\__\n"
"  _/_/          \\__\n"
"  |               |\n"
"  |             | |\n"
"  |             | |\n"
"  |             | |\n"
"  |             | |\n"
"  |___            _               Cl\n"
"    \\_\\__       _/ \\_          __\n"
"       \\_\\_  __/     \\__    __/\n"
"         \\__/           \\__/\n"
"                         | |\n"
"                         | |\n"
"                         | |\n"
"                         | |\n"
"  \n"
"                          O\n\n"

    "If the image appears elongated or squat, the aspect ratio should be changed\n"
    "from its default value of 1.5 using the ``-xa <ratio>`` option. To help\n"
    "determine the correct value, use the ``-xs`` option to display a square.\n"

    "Write Options e.g. -xw 69\n"
    
    " w <characters> Image width in characters, default 79\n"
    " h <characters> Image height in characters, default is width/aspect\n"
    " a <ratio> Aspect ratio of character height:width, default is 1.5\n"
    " s         Display a square - this is useful for correcting the aspect ratio\n"
    " t         Write the output molecule index and the title\n"
    " m         Include a margin around the depiction\n\n"
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
ASCIIFormat theASCIIFormat;

/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////

bool ASCIIFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
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
      obErrorLog.ThrowError("ASCIIFormat", "gen2D not found", obError, onceOnly);
      return false;
    }
    if(!pOp->Do(&workingmol))
    {
      obErrorLog.ThrowError("ASCIIFormat", string(workingmol.GetTitle()) + "- Coordinate generation unsuccessful", obError);
      return false;
    }
  }
  if(!workingmol.Has2D() && workingmol.NumAtoms()>1)
  {
    string mes("Molecule ");
    mes += workingmol.GetTitle();
    mes += " needs 2D coordinates to display in ASCIIFormat";
    obErrorLog.ThrowError("ASCIIFormat", mes, obError);
    return false;
  }

  // Default is width is 79, aspect is 1.5, height is width/aspect
  const char* pp = pConv->IsOption("w");
  int width  = pp ? atoi(pp) : 79;
  pp = pConv->IsOption("a");
  double aspect  = pp ? atof(pp) : 1.5;

  pp = pConv->IsOption("h");
  int height  = pp ? atoi(pp) : static_cast<int>(0.5 + width/aspect);

  if (pConv->IsOption("t")) {
    ofs << "#" << pConv->GetOutputIndex() << " " << pmol->GetTitle() << endl;
  }

  ASCIIPainter painter(width, height, aspect);
  OBDepict depictor(&painter);
  // Don't use margin...unless someone specifies it
  if(!pConv->IsOption("m"))
    depictor.SetOption(OBDepict::noMargin);

  if(pConv->IsOption("s")) { // Print the test card -- should be a square
    painter.NewCanvas(100, 100);
    painter.DrawLine(20, 20, 80, 20);
    painter.DrawLine(80, 20, 80, 80);
    painter.DrawLine(80, 80, 20, 80);
    painter.DrawLine(20, 80, 20, 20);
  }
  else
    depictor.DrawMolecule(&workingmol);

  painter.Write(ofs);
  if(pConv->IsOption("s")) {
    ofs << "The above drawing is supposed to show a square. "
        << "If instead you see a squat rectangle, try again with a smaller aspect ratio, e.g.\n   -oascii -xs -xa " << aspect-0.1 << "\n"
        << "If you see a tall rectangle, try again with a larger aspect ratio, e.g.\n   -oascii -xs -xa " << aspect+0.1 << "\n";
  }

  return true;
}

} //namespace OpenBabel

