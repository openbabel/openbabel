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
#include <openbabel/mol.h>
#include <openbabel/op.h>
#include <openbabel/depict/depict.h>
#include <openbabel/depict/cairopainter.h>
#include <openbabel/alias.h>

#include <cstdlib>

using namespace std;



namespace OpenBabel
{

class PNG2Format : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  PNG2Format() : _ncols(0), _nrows(0), _nmax(0)
  {
    OBConversion::RegisterFormat("_png2",this);
  }

  virtual const char* Description() //required
  {
    return
    "PNG2 format\n"
    "An internal format to write 2D depiction using Cairo\n"
    "Called from PNGFormat\n";
  };


  virtual unsigned int Flags()
  {
      return NOTREADABLE | WRITEBINARY | DEPICTION2D;
  };
  bool WriteChemObject(OBConversion* pConv);
  bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    int _ncols, _nrows, _nmax;
    vector<OBBase*> _objects;
    CairoPainter _cairopainter; //now local variable in WriteMolecule()
};
  ////////////////////////////////////////////////////

//Make an instance of the format class
PNG2Format thePNG2Format;

/////////////////////////////////////////////////////////////////

bool PNG2Format::WriteChemObject(OBConversion* pConv) // Taken from svgformat.cpp
{
  //Molecules are stored here as pointers to OBBase objects, which are not deleted as usual.
  //When there are no more they are sent to WriteMolecule.
  //This allows their number to be determined whatever their source
  //(they may also have been filtered), so that the table can be properly dimensioned.

  OBBase* pOb = pConv->GetChemObject();

  if(pConv->GetOutputIndex()<=1)
  {
    _objects.clear();
    _nmax=0;
    //_ncols = _nrows = 1;
    pConv->AddOption("pngwritechemobject"); // to show WriteMolecule that this function has been called
    const char* pc = pConv->IsOption("c");
    const char* pr = pConv->IsOption("r");
    if(pr)
      _nrows = atoi(pr);
    if(pc)
      _ncols = atoi(pc);
    if(pr && pc) // both specified: fixes maximum number objects to be output
      _nmax = _nrows * _ncols;

    //explicit max number of objects
    const char* pmax =pConv->IsOption("N");
    if(pmax)
      _nmax = atoi(pmax);
  }

  OBMoleculeFormat::DoOutputOptions(pOb, pConv);

  //save molecule
  _objects.push_back(pOb);

  bool ret=true;
  //Finish if no more input or if the number of molecules has reached the allowed maximum(if specified)
  bool nomore = _nmax && (_objects.size()==_nmax);
  if((pConv->IsLast() || nomore))
  {
    int nmols = _objects.size();
    //Set table properties according to the options and the number of molecules to be output
    if(!(nmols==0 ||                       //ignore this block if there is no input or
        (_nrows && _ncols) ||              //if the user has specified both rows and columns or
        ((!_nrows && !_ncols) && nmols==1))//if neither is specified and there is one output molecule
      )
    {
      if(!_nrows && !_ncols ) //neither specified
      {
        //assign cols/rows in square
        _ncols = (int)ceil(sqrt(((double)nmols)));
      }

      if(_nrows)
        _ncols = (nmols-1) / _nrows + 1; //rounds up
      else if(_ncols)
        _nrows = (nmols-1) / _ncols + 1;
    }

    //output all collected molecules
    int n=0;

    vector<OBBase*>::iterator iter;
    for(iter=_objects.begin(); ret && iter!=_objects.end(); ++iter)
    {
      //need to manually set these to mimic normal conversion
      pConv->SetOutputIndex(++n);
      pConv->SetLast(n==_objects.size());

      ret=WriteMolecule(*iter, pConv);

    }

    //delete all the molecules
    for(iter=_objects.begin();iter!=_objects.end(); ++iter)
      delete *iter;

    _objects.clear();
    _nmax = _ncols = _nrows = 0;
  }
  //OBConversion decrements OutputIndex when returns false because it thinks it is an error
  //So we compensate.
  if(!ret || nomore)
    pConv->SetOutputIndex(pConv->GetOutputIndex()+1);
  return ret && !nomore;
}


////////////////////////////////////////////////////////////////

bool PNG2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream& ofs = *pConv->GetOutStream();

  OBMol workingmol(*pmol); // Copy the molecule

  // CairoPainter cairopainter;

  if (!pConv->IsOption("pngwritechemobject") || (!_nrows && !_ncols))
  { //If WriteMolecule called directly, e.g. from OBConversion::Write()
    _nmax = _nrows = _ncols = 1;
    pConv->SetLast(true);
    pConv->SetOutputIndex(1);
  }

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

  pp = pConv->IsOption("w"); // Default values for width and height are the size
  int width  = pp ? atoi(pp) : size;
  pp = pConv->IsOption("h");
  int height  = pp ? atoi(pp) : size;

  bool transparent=false;
  string background, bondcolor;
  const char* bg = pConv->IsOption("b");
  background = bg ? "black" : "white";
  bondcolor  = bg ? "white" : "black";
  if(bg && (!strcmp(bg, "none") || bg[0]=='0'))
  {
    transparent = true;
    bondcolor = "gray";
  }
  const char* bcol = pConv->IsOption("B");
  if(bcol && *bcol)
    bondcolor = bcol;
  if(bg && *bg)
    background = bg;

  string text;
  if(!pConv->IsOption("d"))
  {    
    text = pmol->GetTitle();
    _cairopainter.SetTitle(text);
  }

  if(pConv->GetOutputIndex()==1) {
    _cairopainter.SetWidth(width);
    _cairopainter.SetHeight(height);
    _cairopainter.SetTableSize(_nrows, _ncols);
  }
  _cairopainter.SetIndex(pConv->GetOutputIndex());

  // Detect if cropping should be done, also remove title in that case...
  if((pConv->GetOutputIndex()==1) && pConv->IsLast() && pConv->IsOption("m")) {
    _cairopainter.SetCropping(true);
    _cairopainter.SetTitle("");
  }

  OBDepict depictor(&_cairopainter);

  // The following options are all taken from svgformat.cpp
  if(!pConv->IsOption("C"))
    depictor.SetOption(OBDepict::drawTermC);
  if(pConv->IsOption("a"))
    depictor.SetOption(OBDepict::drawAllC);

  if(pConv->IsOption("A"))
  {
    AliasData::RevertToAliasForm(workingmol);
    depictor.SetAliasMode();
  }
  _cairopainter.SetBondColor(bondcolor);
  depictor.SetBondColor(bondcolor);
  _cairopainter.SetBackground(background);
  _cairopainter.SetTransparent(transparent);
  if(pConv->IsOption("t"))
    _cairopainter.SetPenWidth(4);
  else
    _cairopainter.SetPenWidth(1);

  //No element-specific atom coloring if requested
  if(pConv->IsOption("u"))
    depictor.SetOption(OBDepict::bwAtoms);
  if(!pConv->IsOption("U"))
    depictor.SetOption(OBDepict::internalColor);
  if(pConv->IsOption("s"))
    depictor.SetOption(OBDepict::asymmetricDoubleBond);

  // Draw it!
  depictor.DrawMolecule(&workingmol);

  if (pConv->IsLast())
  {
    if(!pConv->IsOption("O"))
     //if no embedding just write image
      _cairopainter.WriteImage(ofs);
    else //embedding
    {
      //write image to stringstream, read it into pngformat
      stringstream ss;
      _cairopainter.WriteImage(ss);
      OBConversion conv2(&ss,pConv->GetOutStream());
      conv2.CopyOptions(pConv);
      OBBase Ob; //dummy
      OBFormat* ppng = OBConversion::FindFormat("png");
      if(  !conv2.SetInAndOutFormats(ppng, ppng)
        || !ppng->ReadMolecule(&Ob , &conv2))
      {
        obErrorLog.ThrowError("PNG Format", "Failed to embed molecule(s)",obError);
        _cairopainter.WriteImage(ofs); //just write image without embedding
        return true;
      }
      
      // Write each molecule to png format which will embed
      vector<OBBase*>::iterator iter;
      bool ret=true;
      for(iter=_objects.begin(); ret && iter!=_objects.end(); ++iter)
      {
        conv2.SetLast(iter==_objects.end()-1);
        ret=ppng->WriteMolecule(*iter, &conv2);
      }
      if(ret)
        ofs << ss.rdbuf(); //copy png (with embeds) to normal output stream
    }
  }
  return true; //or false to stop converting
}

} //namespace OpenBabel

