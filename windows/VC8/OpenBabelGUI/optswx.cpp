/**********************************************************************
optswx.cpp -  Constructs wxWidgets Option controls from description text

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
#include "stdwx.h"
#include <sstream>
#include <utility>
#include "optswx.h"

////////////////////////////////////////////////////
DynOptionswx::~DynOptionswx()
{ Clear();}

void DynOptionswx::Clear()
{
	OMapType::iterator itr;
	for(itr=OptionMap.begin();itr!=OptionMap.end();++itr)
	{
		sizer->Detach(itr->second);
		(itr->second)->Destroy();
	}
	OptionMap.clear();
	std::vector<wxSizer*>::iterator itrs;
	for(itrs=Sizers.begin();itrs!=Sizers.end();++itrs)
	{
		sizer->Detach(*itrs);
		delete *itrs;
	}
	Sizers.clear();
}

bool DynOptionswx::Construct(const char* OptionsText, const char* StartText, int MultiCharFilter)
{
	//Looks for "options" (case-insensitive)
	//If StartText is not NULL it must precede "options" on the same line and 
	//not be part of a longer word
	const int ONE=1;
	const int FOUR=4;
	bool NextIsRadio=false;
	char* pNewStr = new char[strlen(OptionsText)+1];
	strcpy(pNewStr,OptionsText); //Make a working copy 
	char* p = pNewStr;
	
	bool OptionsFound=false;
	char* lineend = NULL;
	
	if(StartText)
	{
		do
		{
			p=strcasestr(p,StartText);//locate StartText
			if(!p) break;
			p += strlen(StartText);
		}while(isalpha(*p));//next char is not a letter or number
		if(p)
			lineend=strchr(p,'\n');
	}
	if(p && (p=strcasestr(p,"option")) && (!lineend  || p<lineend))
	{
		OptionsFound=true;

		p = strchr(p,'\n')+1; //options start on next line		
		while(p && *p && *p!='\n') //loop for all options until blank line
		{
			int ProvideEditCtl=0;
			while(*p && !isalnum(*(p++))); //skip space and punctuation
			if(!(*p--)) break;
			wxString oname;
			while(isalnum(*p))
				oname += *p++;

			//Filter on whether option is single or multicharacter
			//MultiCharFilter == 0 Display all options
			//MultiCharFilter == 1 Display only single char options
			//MultiCharFilter == 2 Display only multi char options
			if((MultiCharFilter==1 && oname.size()>1)  
				|| (MultiCharFilter==2 && oname.size()==1))  
			{
				p =strchr(p,'\n');
				if(p) ++p; //to next line
				continue;
			}

			while(isspace(*p++));//skip white space
			while(ispunct(*(--p)) && *p != '-')
			{
				//If first non-white space character after option char is punctuation(except for '-')
				//Provide an edit control rather than a checkbox
				++ProvideEditCtl;
				char endch=0;
				switch(*p++)
				{
				case '\"':
				case '\'':
					endch=*(p-1); break;
				case '<':
					endch='>'; break;
				case '[':
					endch= ']'; break;
				case '(':
					endch=')'; break;
				}
				if(endch)
				{
					char* ptemp = strchr(p, endch);
					if(ptemp)
						p=ptemp+1;
				}
				while(isspace(*p++));//skip white space
				if(!(*p)) break;
			}
			p--; //so not to lose first char of caption

			while(!isalnum(*(++p))); //skip space and punctuation
			if(!(*p)) break;
			char* pCaption = p-1;
			p =strchr(p,'\n');
			if(p)
				*p++ ='\0'; //mark end of this option's text

			if(ProvideEditCtl)
			{
				//First the caption
				char* pdef;
				char* pdefWord=NULL;
				if(pdef=strcasestr(pCaption,"default"))
				{
					//Put the next word in the editbox
					char* tok=strtok(pdef," :-\t</\"\'");
					pdefWord = strtok(NULL," :-</\t;,\"\'>");

					//delete caption after default or after <default etc
					*pdef='\0';
					while(isspace(*(--pdef)));
					if(!strpbrk(pdef,"([<{-;")) pdef++ ;
					*pdef='\0';
				}
				wxStaticText* pEdCaption = new wxStaticText(parent,wxID_STATIC,pCaption);
				OptionMap.push_back(std::make_pair(wxString(),pEdCaption));//string is empty for a caption: not an option

				//Edit boxes for multicharacter options are larger
				const int EDWIDTH = oname.size()>1? 60 : 40;
				wxTextCtrl* pEd;

				//Make a large edit box on the next line if last char is of the caption is :
				bool BigEdit =(pCaption[strlen(pCaption)-1] == ':'); 
				if(BigEdit)
				{
					sizer->Add(pEdCaption,0,wxEXPAND|wxTOP,FOUR);
					while(ProvideEditCtl--)
					{
						pEd = new wxTextCtrl(parent,wxID_ANY,wxEmptyString,
								wxDefaultPosition,wxSize(EDWIDTH,16));
						OptionMap.push_back(std::make_pair(oname,pEd));
						if(ProvideEditCtl)
							oname = ' ' + oname;//editboxes except the first have name preceded by one or more spaces
					}
					sizer->Add(pEd,0,wxEXPAND|wxTOP,ONE);
				}
				else
				{
					wxBoxSizer* pEdSizer = new wxBoxSizer(wxHORIZONTAL);
					while(ProvideEditCtl--)
					{
						pEd = new wxTextCtrl(parent,wxID_ANY,wxEmptyString,
								wxDefaultPosition,wxSize(EDWIDTH,16));
						OptionMap.push_back(std::make_pair(oname,pEd));
						if(ProvideEditCtl)
							oname = ' ' + oname;//editboxes except the first have name preceded by one or more spaces
						pEdSizer->Add(pEd,0);
					}
					pEdSizer->Add(pEdCaption,1,wxLEFT|wxALIGN_CENTER_VERTICAL,FOUR);
					sizer->Add(pEdSizer,0,wxEXPAND|wxTOP,FOUR);
					Sizers.push_back(pEdSizer);
				}
				pEd->AppendText(pdefWord);
			}
			else
			{
				// Checkbox, or radio button ("or" is the last word in the caption && no letters after " or")
				wxControl* pChk;
				
				//First see if should be checked
				//If 'default appears in caption set the checkbox unless it is followed
				//by one of 'no' 'not' or 'none'
				char* pdef=strcasestr(pCaption,"default");
				bool SetChk=(pdef!=NULL);
				if(SetChk)
					strcasecmp(pdef,"no") && !strcasecmp(pdef,"not") && strcasecmp(pdef,"none");

				char* por=strcasestr(pCaption," or");
				bool HasOr = (por && !strpbrk(por+3,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"));

				if(NextIsRadio || HasOr)
				{
					unsigned int style = NextIsRadio ? 0 : wxRB_GROUP; //Style only on first radiobutton
					pChk = new wxRadioButton(parent,wxID_ANY,pCaption,
								wxDefaultPosition, wxDefaultSize, style);
					NextIsRadio = HasOr;
					((wxRadioButton*)pChk)->SetValue(SetChk);
				}
				else
				{
					pChk = new wxCheckBox(parent,wxID_ANY,pCaption);
					((wxCheckBox*)pChk)->SetValue(SetChk);
				}
			
				OptionMap.push_back(std::make_pair(oname,pChk));
				sizer->Add(pChk,0,wxEXPAND|wxTOP,FOUR);
			}
		}
	}
	if(OptionsFound)
	{
		wxStaticLine* pLine = new wxStaticLine(parent);
		sizer->Add(pLine,0,wxTOP|wxEXPAND,FOUR);
		OptionMap.push_back(std::make_pair(wxString(),pLine));//empty string
	}
	delete [] pNewStr;
	return OptionsFound;
}

//////////////////////////////////////////////////////////


int DynOptionswx::SetOptions(OpenBabel::OBConversion& Conv, OpenBabel::OBConversion::Option_type opttyp)
{
	//Now sets options directly in OBConversion
	int count=0;
	OMapType::iterator itr;
	for (itr = OptionMap.begin(); itr != OptionMap.end(); ++itr)
	{
		if(itr->first.empty()) continue; //just a caption or a line
		wxCheckBox* pChk = dynamic_cast<wxCheckBox*> (itr->second);
		if(pChk)
		{	
			if(pChk->IsChecked())
			{
				Conv.AddOption(itr->first, opttyp);
				++count;
			}
		}
		else
		{
			wxRadioButton* pRadio = dynamic_cast<wxRadioButton*> (itr->second);
			if(pRadio)
			{
				if(pRadio->GetValue())
				{
					Conv.AddOption(itr->first, opttyp);
					++count;
				}
			}
			else
			{
				wxTextCtrl* pText = dynamic_cast<wxTextCtrl*> (itr->second);
				if(pText)
				{
					wxString txt = pText->GetValue();
					if(txt.IsEmpty()) continue;
					wxString oname = itr->first;

					//Get the contents of subsequent editboxes 
					OMapType::iterator itr2 = itr;
					while(++itr2!= OptionMap.end())
					{
						if(itr2->first[0]!=' ') //subsequent editboxes have the name preceded by a space
							break;
						txt = txt + ' ' + static_cast<wxTextCtrl*>(itr2->second)->GetValue();
						++itr;
					}
					Conv.AddOption(oname, opttyp, txt);
					++count;
				}
			}
		}
	}
	return count;
}

//////////////////////////////////////////////////////////
char* DynOptionswx::strcasestr(const char* haystack, const char* needle)
{
	//Adapted from http://primates.ximian.com/~fejj/strlib.c
	register unsigned char *h, *n, *hc, *nc;
	size_t needlelen = strlen (needle);
	if (needlelen == 0)
		return (char *) haystack;

	h=(unsigned char *)haystack;
	n=(unsigned char *)needle;
	while (*(h + needlelen - 1))
	{
		if (tolower(*h) == tolower(*n))
		{
			for (hc = h + 1, nc = n + 1; *hc && *nc; hc++, nc++)
				if (tolower(*hc) != tolower(*nc))
					break;
			
			if (!*nc)
				return (char*) h;
		}		
		h++;
	}
	return NULL;
}
