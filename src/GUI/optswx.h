/**********************************************************************
optswx.h -  Constructs wxWidgets Option controls from description text

Copyright (C) 2006 by Chris Morley

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
#ifndef OB_OPTSWX_H
#define OB_OPTSWX_H

#include <vector>
#include <openbabel/obconversion.h>

class DynOptionswx
{
public:
	typedef std::vector <std::pair<wxString, wxWindow*> > OMapType;

	DynOptionswx(wxWindow* Par, wxSizer* pSizer)
		: parent(Par), sizer(pSizer){}
	~DynOptionswx();
	void Clear();
	bool Construct(const char* OptionsText, const char* StartText=NULL, int MultiCharFilter=0);
	int SetOptions(OpenBabel::OBConversion& Conv, OpenBabel::OBConversion::Option_type opttyp);
private:
	char* strcasestr(const char* haystack, const char* needle);
	wxWindow* parent;
	wxSizer* sizer;
	OMapType OptionMap;
	std::vector<wxSizer*> Sizers;
};
#endif
