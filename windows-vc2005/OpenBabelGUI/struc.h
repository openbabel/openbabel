/**********************************************************************
struc.h -  Display of chemical structure

Copyright (C) 2006 by Chris Morley

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
#ifndef OB_STRUC_H
#define OB_STRUC_H

#include "stdwx.h"

class OBBase; //forward declaration

class StructureDisplay : public wxWindow
{
public:
	boolDisplay(OBBase* pOb);
	void OnPaint(wxPaintEvent event);
private:
	DECLARE_EVENT_TABLE()

private:
	OBBase* pOb;
	double Scale, XOffset, YOffset;
	
}
