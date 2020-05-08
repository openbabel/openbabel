/**********************************************************************
smartsfilters.cpp - Some descriptors which analyse molecule using SMARTS

Copyright (C) 2007 by Chris Morley

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
#include <openbabel/oberror.h>
#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include <openbabel/parsmart.h>

using namespace std;
namespace OpenBabel
{

  //**************************************************************
  class SmartsDescriptor : public OBDescriptor
  {
  public:
    //! constructor. Each instance provides an ID, a SMARTS pattern and a description.
    SmartsDescriptor(const char* ID, const char* smarts, const char* descr)
      : OBDescriptor(ID, false), _smarts(smarts), _descr(descr){}

    virtual const char* Description()
    {
      //Adds the SMARTS string to the description
      static string txt;
      txt =  _descr;
      txt += "\n\t SMARTS: ";
      txt += _smarts;
      txt += "\nSmartsDescriptor is definable";
      return txt.c_str();
    }

    double Predict(OBBase* pOb, string* param=nullptr)
    {
      OBMol* pmol = dynamic_cast<OBMol*> (pOb);
      if(!pmol)
        return 0;

      OBSmartsPattern sp;
      if (sp.Init(_smarts) && sp.Match(*pmol))
        return sp.GetUMapList().size();
      else
        return 0.0;
    }

    virtual SmartsDescriptor* MakeInstance(const std::vector<std::string>& textlines)
    {
      return new SmartsDescriptor(textlines[1].c_str(),textlines[2].c_str(),textlines[3].c_str());
    }

  private:
    const char* _smarts;
    const char* _descr;
  };

  //Make global instances

  SmartsDescriptor theHBD("HBD", "[!#6;!H0]","Number of Hydrogen Bond Donors (JoelLib)");
  SmartsDescriptor theHBA1("HBA1", "[$([!#6;+0]);!$([F,Cl,Br,I]);!$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]",
                           "Number of Hydrogen Bond Acceptors 1 (JoelLib)\n"
                           "\t Identification of Biological Activity Profiles Using Substructural\n"
                           "\t Analysis and Genetic Algorithms -- Gillet, Willett and Bradshaw,\n"
                           "\t U. of Sheffield and Glaxo Wellcome.\n"
                           "\t Presented at Random & Rational: Drug Discovery via Rational Design\n"
                           "\t and Combinitorial Chemistry, Strategic Research Institute, Princeton\n"
                           "\t NJ, Sept. 1995" );
  SmartsDescriptor theHBA2("HBA2", "[$([$([#8,#16]);!$(*=N~O);!$(*~N=O);X1,X2]),$([#7;v3;!$([nH]);!$(*(-a)-a)])]",
                           "Number of Hydrogen Bond Acceptors 2 (JoelLib)");
  SmartsDescriptor thenF("nF", "F","Number of Fluorine Atoms");

}//namespace
