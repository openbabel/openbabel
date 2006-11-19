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
#include "babelconfig.h"
#include "struc.h"
#include "mol.h"

BEGIN_EVENT_TABLE(StructureDisplay,wxWindow)
	EVT_PAINT(StructureDisplay::OnPaint)
END_EVENT_TABLE()

//////////////////////////////////////////////////
bool StructureDisplay::Display(OBBase* pObject)
{
	pOb = pObject;
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(!pmol)
		return false;
	
	//Calculate scaling factor based on range of coordinates
	double xmax = -1E20, ymax = -1E20, xmin = 1E20, ymin = 1E20; 
	FOR_ATOMS_OF_MOL(a,*pmol)
	{
		xmax = wxMax(xmax, a.x());
		ymax = wxMax(ymax, a.x());
		xmin = wxMax(xmin, a.x());
		ymin = wxMax(ymin, a.x());
	}
	int width, height, border=10;
	GetClientSize(&width, &height);
	Scale = wxMax((width-2*border)/(xmax-xmin), (height-2*border)/(ymax/ymin));
	XOffset = wxMax(0, (width-(xmax-xmin)*Scale))/2;
	YOffset = wxMax(0, (height-(ymax-ymin)*Scale))/2;;

}

/////////////////////////////////////////////////
StructureDisplay::OnPaint(wxPaintEvent event)
{
	pOb = pObject;
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(!pmol)
		return false;

	wxPaintDC dc(this);
	wxPen widePen(*wxBLACK, 3);
	wxPen ruboutPen(this->GetBackgroundColour());

		//Draw bonds
	FOR_BONDS_OF_MOL(b,*pmol)
	{	
		OBAtom* a1 = b->GetBeginAtom();	
		OBAtom* a2 = b->GetEndAtom();
		wxCoord x1 = ceil(a1->x())* Scale - XOffset;
		wxCoord y1 = ceil(a1->y())* Scale - YOffset;
		wxCoord x2 = ceil(a2->x())* Scale - XOffset;
		wxCoord y2 = ceil(a2->y())* Scale - YOffset;
		switch(b->GetBondOrder())
			{
			case 1:
				dc.SetPen(*wxBLACK_PEN);
				dc.DrawLine(x1, y1, x2, y2);
				break;
			case 2:
				dc.SetPen(widePen);				
				dc.DrawLine(x1, y1, x2, y2);
				dc.SetPen(ruboutPen);				
				dc.DrawLine(x1, y1, x2, y2);
				break;
			case 3:
				break;
			case 4:
				break;
			}
	}

/*	//Draw atoms
	for(i=0;i<nAtoms;i++)
	{
		if(!ShowHAtoms && *(pAtomArray[i].nam)=='H') continue;
		char Txt[7],Txtn[2]={0x90,0};
		strcpy(Txt,pAtomArray[i].nam);
		if(ShowHAtoms)
		{
			if(pAtomArray[i].Hs>0)
			{
				strcat(Txt,"H");
				int n=pAtomArray[i].Hs;
				if(n>1)
				{
					_itoa(n,Txtn,10);
					strcat(Txt,Txtn);
				}
			}
		}
		else
			if(*Txt=='C' && !isalpha(*(Txt+1)) && *(Txt+1)!='.') *Txt='\0';

		dc.TextOut((int)(pAtomArray[i].x*Scale-XOffset)- tm.tmAveCharWidth/2,
			(int)(pAtomArray[i].y*Scale-YOffset) - tm.tmAscent/2 - 3, 
			Txt,strlen(Txt));
	}
*/
}

