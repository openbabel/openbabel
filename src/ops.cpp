/**********************************************************************
ops.cpp - Some simple plugins for OBMol to implement some general options

Copyright (C) 2007 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>

#ifndef OBERROR
 #define OBERROR
#endif


namespace OpenBabel
{
  /* This works ok but duplicates the -c option*/
  class OpCenter : public OBOp
  {
  public:
    OpCenter(const char* ID) : OBOp(ID, false){};
    virtual const char* Description(){ return "Centers coordinates around (0,0,0)"; }
    virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }

    virtual bool Do(OBBase* pOb, OpMap*, const char*)
    {
      OBMol* pmol = dynamic_cast<OBMol*>(pOb);
      if(pmol)
        pmol->Center();
      return true;
    }
  };

  //////////////////////////////////////////////////////
  OpCenter theOpCenter("center"); //Global instance

  //*************************************************************

  /* This class is another example -- use the FROG code to generate 3D coords.
  class OpFrog : public OBOp
  {
  public:
    OpFrog(const char* ID) : OBOp(ID, false){};
    virtual const char* Description()
    { return "Adds 3D coordinates using the FROG code"; }
    virtual bool WorksWith(OBBase* pOb)const
    { return dynamic_cast<OBMol*>(pOb)!=NULL; }

    virtual bool Do(OBBase* pOb, OpMap*, const char*);
  };

  ///////////////////////////////////////////////////////
  //Global instance with option name used 
  // e.g. babel infile.xxx -O outfile.yyy --3Dfrog
  OpFrog theOpFrog("3Dfrog"); 

  //////////////////////////////////////////////////////
  bool OpFrog::Do(OBBase* pOb, OpMap*, const char*)
  { return true; }
  */

}//namespace

//! \file ops.cpp
//! \brief Base plugin class for operations on molecules
