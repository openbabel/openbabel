/**********************************************************************
ChangeCell.cpp - The option --changecell  Change unutcell.

Copyright(C) 2014 by Okhotnikov Kirill

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
#include <iostream>
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include <openbabel/oberror.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <cstdlib>

namespace OpenBabel
{

class OpChangeCell : public OBOp
{
protected:
  class vc_val
  {
  public:  
    bool mult;
    double value;
    vc_val(): mult(false), value(0.0) {};
  };
  
public:
  OpChangeCell(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return "Change cell size:\n"
                                    "     [keepfract];[*]a;[*]b;[*]c\n"
                                    "     Original cell dimensions can be changed to value a, b or c or multiplied with key '*' " ; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
};

/////////////////////////////////////////////////////////////////
OpChangeCell theOpChangeCell("ChangeCell"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpChangeCell::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  std::vector<std::string> vcr;
  tokenize(vcr, OptionText, ";");

  if( (vcr.size() != 3) && (vcr.size() != 4) )
  {  
    obErrorLog.ThrowError(__FUNCTION__, "Invalid number of arguments!" , obWarning);
    return false;
  }  

  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  
  bool keep_fract = false;
  if( vcr[0][0] == 'k' )
  {  
    keep_fract = true;
    vcr.erase(vcr.begin());
  }
  
  if( vcr.size() != 3 )
  {  
    obErrorLog.ThrowError(__FUNCTION__, "Invalid input. Check first argument!" , obWarning);
    return false;
  }  
  
  std::vector<vc_val> vcvs;
  
  vcvs.resize(3);
  
  for(int i = 0; i < vcvs.size(); i++)
  {
    std::string str = vcr[i];
    Trim(str);
    vcvs[i].mult = false;
    if( str[0] == '*' )
    {
      vcvs[i].mult = true;
      str = str.substr(1);
    }
    vcvs[i].value = atof(str.c_str());
    if(  vcvs[i].value == 0 )
    {  
      obErrorLog.ThrowError(__FUNCTION__, "Wrong value \"" + str +"\"" , obWarning);
      return false;
    }
  }  
  
  OBUnitCell * old_cell;
  if ( ! pmol->HasData(OBGenericDataType::UnitCell) )
    old_cell = NULL;
  else         
    old_cell = (OBUnitCell*)pmol->GetData(OBGenericDataType::UnitCell);
  
  if( old_cell == NULL )
  {
    if(  keep_fract )
    {  
      obErrorLog.ThrowError(__FUNCTION__, "Cannot keep fractional coordinates without unit cell!" , obWarning);
      return false;
    }
    for(int i = 0; i < vcvs.size(); i++)
    {
      if( vcvs[i].mult )
      {
        obErrorLog.ThrowError(__FUNCTION__, "Cannot multiply sizes without unit cell!" , obWarning);
        return false;
      }
    }
  }  
  double a, b, c, alpha, beta, gamma;
  if( old_cell != NULL )
  {
    a = old_cell->GetA(); b = old_cell->GetB(); c = old_cell->GetC(); 
    alpha = old_cell->GetAlpha(); beta = old_cell->GetBeta(); gamma = old_cell->GetGamma();     
  }  
  else
  {
    a = 0.0; b = 0.0; c = 0.0;
    alpha = 90.0; beta = 90.0; gamma = 90.0;
  }
  
  OBUnitCell * new_cell = new OBUnitCell();
  new_cell->SetData( (vcvs[0].mult ? a : 1.0) * vcvs[0].value,  
                     (vcvs[1].mult ? b : 1.0) * vcvs[1].value, 
                     (vcvs[2].mult ? c : 1.0) * vcvs[2].value,
                     alpha, beta, gamma );
  if(old_cell == NULL)
    new_cell->SetSpaceGroup(1);
  else
    new_cell->SetSpaceGroup(old_cell->GetSpaceGroupNumber());
  
  if( keep_fract )
  {
    pmol->BeginModify();
    
    FOR_ATOMS_OF_MOL(a, pmol)
    {
      vector3 old_fract = old_cell->CartesianToFractional(a->GetVector());
      vector3 new_cart = new_cell->FractionalToCartesian(old_fract);
      a->SetVector(new_cart);
    }        
    
    pmol->EndModify();
  }  
  
  if (old_cell != NULL)
    pmol->DeleteData(old_cell);
  
  pmol->SetData(new_cell);
    
  return true;
}
}//namespace


