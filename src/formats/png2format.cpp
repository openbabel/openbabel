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
// Derive directly from OBFormat for objects which are not molecules.
{
public:
  //Register this format type ID in the constructor
  PNG2Format()
  {
    OBConversion::RegisterFormat("png2",this);

  }

  /* The first line of the description should be a brief identifier, <40 chars, because
     it is used in dropdown lists, etc. in some user interfaces. The rest is optional.

     Describe any format specific options here. This text is parsed to provide
     checkboxes, etc for the GUI (for details click the control menu),
     so please try to keep to a similar form.

     Write options are the most common, and the "Write" is optional.
     The option f takes a text parameter, so that it is essential that the option
     is registered in the constructor of the class.
     Finish the options with a blank line as shown, if there are more than one
     group of options, or if there are further comments after them.
  */
  virtual const char* Description() //required
  {
    return
    "XXX format\n"
    "Some comments here, on as many lines as necessay\n"
    "Write Options e.g. -xf3 \n"
    "  f# Number of (fictional) levels\n"
    "  n  Omit (virtual) title\n\n"
    
    "Read Options e.g. -as\n"
    "  s  Consider single bonds only\n"
    ;
  };

  //Optional URL where the file format is specified
  virtual const char* SpecificationURL(){return
     "http://www.mdl.com/downloads/public/ctfile/ctfile.jsp";};

  //Optional
  virtual const char* GetMIMEType()
  { return "chemical/x-xxx"; };


  /* Flags() can return be any of the following combined by |
     or be omitted if none apply
     NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
     READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
  virtual unsigned int Flags()
  {
      return NOTREADABLE | WRITEBINARY | WRITEONEONLY;
  };

   /* This optional function is for formats which can contain more than one
     molecule. It is used to quickly position the input stream after the nth
     molecule without have to convert and discard all the n molecules.
     See obconversion.cpp for details and mdlformat.cpp for an example.*/
  virtual int SkipObjects(int n, OBConversion* pConv)
  {
    return 0;
  };

  ////////////////////////////////////////////////////
  /// Declarations for the "API" interface functions. Definitions are below
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:
  /* Add declarations for any local function or member variables used.
     Generally only a single instance of a format class is used. Keep this in
     mind if you employ member variables. */
};
  ////////////////////////////////////////////////////

//Make an instance of the format class
PNG2Format thePNG2Format;

/////////////////////////////////////////////////////////////////

bool PNG2Format::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL)
      return false;

  istream& ifs = *pConv->GetInStream();

  pmol->BeginModify();

  /** Parse the input stream and use the OpenBabel API to populate the OBMol **/

  // To use an input option
  if(pConv->IsOption("s",OBConversion::INOPTIONS))
  {
    //Code for when -as is specified
  }

  /* If the molecule has other than 3D coordinates for its atoms, it
  is necessary to set the dimension to 0, or 2 */
  int dim;
  dim = 3;
  pmol->SetDimension(dim);

  pmol->EndModify();

  /* For multi-molecule formats, leave the input stream at the start of the
     next molecule, ready for this routine to be called again.

  /* Return true if ok. Returning false means discard the OBMol and stop
     converting, unless the -e option is set. With a multi-molecule inputstream
     this will skip the current molecule and continue with the next, if SkipObjects()
     has been defined. If it has not, and continuation after errors is still required,
     it is necessary to leave the input stream at the beginning of next object when
     returning false;*/
  return true;
}

////////////////////////////////////////////////////////////////

bool PNG2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream& ofs = *pConv->GetOutStream();

  //*** Coordinate generation ***
  //Generate coordinates only if no existing 2D coordinates
  if(!pmol->Has2D(true))
  {
    OBOp* pOp = OBOp::FindType("gen2D");
    if(!pOp)
    {
      obErrorLog.ThrowError("PNG2Format", "gen2D not found", obError, onceOnly);
      return false;
    }
    if(!pOp->Do(pmol))
    {
      obErrorLog.ThrowError("PNG2Format", string(pmol->GetTitle()) + "- Coordinate generation unsuccessful", obError);
      return false;
    }
  }
  if(!pmol->Has2D() && pmol->NumAtoms()>1)
  {
    string mes("Molecule ");
    mes += pmol->GetTitle();
    mes += " needs 2D coordinates to display in PNG2format";
    obErrorLog.ThrowError("PNG2Format", mes, obError);
    return false;
  }

  CairoPainter painter;
  OBDepict depictor(&painter);
  depictor.DrawMolecule(pmol);
  painter.WriteImage(ofs);

  return true; //or false to stop converting
}

} //namespace OpenBabel

