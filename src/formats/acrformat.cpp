/**********************************************************************
Copyright (C) 2007 by Daniel Mansfield
Some portions Copyright (C) 2004-2006 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/*
 * File extension module for CaRIne's ASCII Crystal (ACR)
 * By Daniel Mansfield
 * 30th January 2007
 */

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/obmolecformat.h>
#include <stdio.h>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class ACRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID in the constructor
    ACRFormat()
    {
      OBConversion::RegisterFormat("acr", this, "chemical/x-acr");
      //		OBConversion::RegisterOptionParam("f", this, 1);
      //		OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);

    }

    virtual const char* Description() //required
    {
      return
        "ACR format\n"
        "CaRIne ASCII Crystal format (ACR)\n"
        //      "Write Options e.g. -xf3 \n"
        // "  f# Number of (fictional) levels \n"
        //			"  n  Omit (virtual) title\n"
        "Read Options e.g. -as\n"
        "  s  Consider single bonds only\n";
    };

    virtual const char* SpecificationURL()
    {return "http://pros.orange.fr/carine.crystallography/books/31/carine_31_us.pdf";};

    virtual const char* GetMIMEType() { return "chemical/x-acr"; };


	  virtual unsigned int Flags()
	  {
      return READONEONLY | NOTWRITABLE;
	  };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    };

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
  };

  ACRFormat theACRFormat;

  bool ACRFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    istream& ifs = *pConv->GetInStream();

    pmol->BeginModify();

    /** Parse the input stream and use the OpenBabel API to populate the OBMol **/
    char buf[BUFF_SIZE];
    unsigned int atoms, bonds, tmp;
    float scale, dtmp;
    bool atom_input = false, bond_input = false;
    string type;
    //int from, to;
    double X,Y,Z;
    vector<string> vs;

    // read in one at a time
    /* WARNING: Atom id starts from zero in Carine; not so in openbabel.
     * Solution: Let Open Babel to set them. */

    while (true) {
      ifs.getline(buf, BUFF_SIZE);
      if (ifs.eof()) {
        break;
      }

      if (sscanf(buf, "General Scale=%f\n", &dtmp)) {
        scale = dtmp;
        continue;
      } else if (sscanf(buf, "Number of Atoms in Crystal=%d\n", &tmp)) {
        atoms = tmp;
        atom_input = true;

        // read table column names
        ifs.getline(buf, BUFF_SIZE);
        continue;
      } else if (sscanf(buf, "Number of Links in Crystal=%d\n", &tmp)) {
        atom_input = false;
        bond_input = true;
        bonds = tmp;

        // read table column names
        ifs.getline(buf, BUFF_SIZE);
        continue;
      } else if ( '#' == buf[0] || '\r' == buf[0] || '\n' == buf[0] ) {
        // between sections, in both windows and unix.
        continue;
      }
      tokenize(vs, buf, " \t\r\n");

      if (atom_input) {
	if (vs.size() < 9) return false; // timvdm 18/06/2008
        type = vs[1];
        X = atof((char*)vs[6].c_str())/scale;
        Y = atof((char*)vs[7].c_str())/scale;
        Z = atof((char*)vs[8].c_str())/scale;

        OBAtom* a = pmol->NewAtom();
        if (*(type.c_str()) != '*')
          a->SetAtomicNum(OBElements::GetAtomicNum(type.c_str()));
        a->SetVector(X,Y,Z);

      } else if (bond_input) {
	if (vs.size() < 2) return false; // timvdm 18/06/2008
        // add to pmol
        if (!pmol->AddBond(atoi((char*)vs[0].c_str()) + 1, atoi((char*)vs[1].c_str()) + 1,
                           1 /* bond order not specified in Carine, use PerceiveBondOrder later */))
          {
            obErrorLog.ThrowError(__FUNCTION__, "addition of bond between " + vs[0] + " and " + vs[1] + " failed", obError);
            return false;
          }
      }
    }

    /* got sanity? */
    if ( pmol->NumBonds() != bonds ) {
      // then we read a different number of bonds than those promised.
      obErrorLog.ThrowError(__FUNCTION__, "Number of bonds read does not match the number promised", obError);
      return false;
    } else if ( pmol->NumAtoms() != atoms ) {
      obErrorLog.ThrowError(__FUNCTION__, "Number of atoms read does not match the number promised", obError);
      return false;
    }

    pmol->PerceiveBondOrders();

    pmol->EndModify();

    return true;
  }

} //namespace OpenBabel

