/**********************************************************************
titleformat.h - Displays and reads molecule titles

Copyright (C) 2005 Chris Morley

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

using namespace std;
namespace OpenBabel {

class TitleFormat : public OBMoleculeFormat
{
public:
	TitleFormat()
  {
      OBConversion::RegisterFormat("txt",this);
  }

  virtual const char* Description() //required
  {
    return
      "Title format\n"
      "Displays and reads molecule titles\n";
  }
	virtual unsigned int Flags() { return ZEROATOMSOK; }

  /// The "API" interface functions
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};
/////////////////////////////////////////////////////
TitleFormat theTitleFormat;

/////////////////////////////////////////////////////
bool TitleFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
	// Reads titles separated by spaces, tabs or newlines,
	// If option -at set titles can contain spaces.
  OBMol* pmol = pOb->CastAndClear<OBMol>();
	string title;
	istream& ifs = *pConv->GetInStream();
	if(pConv->IsOption("t",OBConversion::INOPTIONS))
	{
		while( ifs && (ifs.peek()!='\t') && (ifs.peek()!='\n') && (ifs.peek()!=EOF))
			title += ifs.get();
		ifs.get(); //delimiter
	}
	else
		ifs >> title;

	pmol->SetTitle(Trim(title));
	return true;
}

bool TitleFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream &ofs = *pConv->GetOutStream();
	ofs << pmol->GetTitle() << endl;
	return true;
}

}//namespace OpenBabel
