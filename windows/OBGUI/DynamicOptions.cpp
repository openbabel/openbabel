// DynamicOptions.cpp: implementation of the CDynamicOptions class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning (disable : 4786)

#include "stdafx.h"
#include "OBGUI.h"
#include "DynamicOptions.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

CDynamicOptions::CDynamicOptions()
{
	nID=5678;
}

void CDynamicOptions::Construct(const char* OptionTxt, CWnd* pWnd, CRect& CtlRect )
{
/* Reads OptionTxt and starting with the line following each occurrence of "Options"
interprets each line like
  w  use Welsh language (default)
as an option. Constructs a checkbox with the text except for the first character(w in this case)
and inserts it into the window starting at the specified position (window coordinates). Subsequent
checkboxes are displaced downwards by the height of CtlRect. At end, CtlRect is beyond the last chkbox.
When CDynamicOptions::GetOptions() is called it returns a string containing the first characters	
*/

	Clear(); //make sure there are no previous checkboxes

/*	char* txt=strcasestr(OptionTxt,"option"); 
	if(txt==NULL) return;
	txt = strchr(txt,'\n');
	if(txt==NULL) return;
	char* pNewStr = new char[strlen(txt)];
	char* p = pNewStr;
	strcpy(p,++txt); //make a copy of the useful text
*/
	char* pNewStr = new char[strlen(OptionTxt)+1];
	strcpy(pNewStr,OptionTxt); //Make a working copy 
	char* p = pNewStr;

	while(p && (p=strcasestr(p,"option"))) //until end of text; multiple sets of Options found
	{
		bool NextIsRadio=false;

		p = strchr(p,'\n'); //options start on next line
		while(p && *(++p)!='\n') //loop for all options until blank line
		{
			bool ProvideEditCtl=false;
			while(*p && !isalnum(*(p++))); //skip space and punctuation
			if(!(*p)) break;
			char ochar = *(--p); //The option character
			while(isspace(*(++p)));//skip white space
			if(ispunct(*p) && *p != '-')
			{
				//If first non-white space character after option char is punctuation(except for '-')
				//Provide an edit control rather than a checkbox
				ProvideEditCtl=true;
				while(!isspace(*(++p))); //ignore chars until next space
				if(!(*p)) break;
			}
			else
				p--; //so not to lose first char of caption

			while(!isalnum(*(++p))); //skip space and punctuation
			if(!(*p)) break;
			char* pCaption = p-1;
			p =strchr(p,'\n');
			if(p)
				*p++ ='\0'; //mark end of this option's text

			if(ProvideEditCtl)
			{
				const int EDWIDTH = 32;
				const int CHKTEXTINDENT = 19;
				CEdit* pEd = new CEdit;

				//Make a large edit box on the next line if last char is of the caption is :
				bool BigEdit =(pCaption[strlen(pCaption)-1] == ':'); 

				//First the caption
				char* pdef;
				char* pdefWord=NULL;
				if(pdef=strcasestr(pCaption,"default"))
				{
					//Put the next word in the editbox
					char* tok=strtok(pdef," :-\t</\"\'");
					pdefWord = strtok(NULL," :-</\t;.,\"\'>");

					//delete caption after default or after <default etc
					*pdef='\0';
					while(isspace(*(--pdef)));
					if(!strpbrk(pdef,"([<{-:;")) pdef++ ;
					*pdef='\0';
				}

				CStatic* ptxt = new CStatic;
				CRect CapRect(CtlRect);
				CapRect.InflateRect(BigEdit ? -CHKTEXTINDENT : -EDWIDTH,0); 
				VERIFY(ptxt->Create(pCaption, WS_CHILD | WS_VISIBLE | SS_LEFTNOWORDWRAP , 
					CapRect, pWnd, nID++));  
				ptxt->SetFont(pWnd->GetFont(),FALSE); //same font as parent

				//Now the edit
				if(BigEdit)
					CtlRect.OffsetRect(0, CtlRect.Height()-5);//edit is on next line

				CRect EdRect(CtlRect.TopLeft(), 
					CSize(BigEdit ? CtlRect.Width() : EDWIDTH,  CtlRect.Height()-3));

				VERIFY(pEd->Create(WS_CHILD | WS_VISIBLE | WS_TABSTOP | WS_BORDER ,
					EdRect, pWnd, nID++));
				pEd->SetFont(pWnd->GetFont(),FALSE); //same font as parent
				if (pdefWord)
					pEd->SetWindowText(pdefWord);
				
				OptionMap.insert(OMapType::value_type(ochar,pEd));
				OptionMap.insert(OMapType::value_type(' ',ptxt));
			}
			else
			{
				CButton* pChk = new CButton();
				//Decide whether it should be a radio button ("or" is the last word in the caption)
				//or a checkbox
				UINT style;
				char* por=strcasestr(pCaption," or");
				if(por // && no letters after " or"
					&& !strpbrk(por+3,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"))
				{
					//This and the next should be radio buttons
					style = BS_AUTORADIOBUTTON;
					if (!NextIsRadio)
						style |= WS_GROUP; //set only for first of group
					NextIsRadio=true;
					*p='\0';
				}
				else
				{
					if(!NextIsRadio)
						style= BS_AUTOCHECKBOX;
					else
						style=BS_AUTORADIOBUTTON; //no WS_GROUP
					NextIsRadio=false;
				}
				
				int SetChk=0;
				//If 'default appears in caption set the checkbox unless it is followed
				//by one of 'no' 'not' or 'none'
				char* pdef;
				if(pdef=strcasestr(pCaption,"default"))
				{
					SetChk=1;
					char*tok=strtok(pdef," <>()[]-,");
					while(tok=strtok(NULL," <>()[]-,"))
					{
						if(!stricmp(tok,"no") || !strcmp(tok,"not") || !strcmp(tok,"none"))
							SetChk=0;
					}
					
					//delete caption after default or after <default etc
					*pdef='\0';
					while(isspace(*(--pdef)));
					if(!strpbrk(pdef,"([<{-:;")) pdef++ ;
					*pdef='\0';
				}
				VERIFY(pChk->Create(pCaption, style | WS_CHILD | WS_VISIBLE | WS_TABSTOP ,
					CtlRect, pWnd, nID++));
				pChk->SetFont(pWnd->GetFont(),FALSE); //same font as parent
				pChk->SetCheck(SetChk);

				OptionMap.insert(OMapType::value_type(ochar,pChk));
			}
			CtlRect.OffsetRect(0, CtlRect.Height());//move active rect down
		}
	}
	delete [] pNewStr;
}

//////////////////////////////////////////////
CDynamicOptions::~CDynamicOptions()
{
	Clear();
}

/////////////////////////////////////////////
CDynamicOptions::Clear()
{
	//Deletes all the checkboxes
	OMapType::iterator itr;
	for (itr = OptionMap.begin(); itr != OptionMap.end(); itr++)
		delete itr->second;
	OptionMap.clear();
}

////////////////////////////////////////////
void CDynamicOptions::InsertText(const char* txt, CWnd* pWnd, CRect& CtlRect )
{
	//Puts a single line of text (a comment) into the window
	CStatic* ptxt = new CStatic;
	CRect CapRect(CtlRect);
	VERIFY(ptxt->Create(txt, WS_CHILD | WS_VISIBLE, CtlRect, pWnd, nID++));  
	ptxt->SetFont(pWnd->GetFont(),FALSE); //same font as parent
  CString wt;
	ptxt->GetWindowText(wt);
	CtlRect.OffsetRect(0, CtlRect.Height());//move active rect down
	OptionMap.insert(OMapType::value_type(' ',ptxt));
}

////////////////////////////////////////////
const char* CDynamicOptions::GetOptions()
{
	//Returns a string with the option characters of the set checkboxes
	//or option letter "text" for each non-empty editbox
  //e.g  ab"Mon"de"Fri"g

	opts="";
	OMapType::iterator itr;
	for (itr = OptionMap.begin(); itr != OptionMap.end(); itr++)
	{
		if(itr->first==' ') continue; //just a caption
		CButton* pChk = dynamic_cast<CButton*> (itr->second);
		if(pChk)
		{
			if(pChk->GetCheck())
				opts += itr->first;
		}
		else
		{
			CString txt;
			(itr->second)->GetWindowText(txt);
			if(txt.IsEmpty()) continue;
			opts += itr->first;
			opts += '\"';
			opts += txt;
			opts += '\"';
		}						
	}
	return opts.c_str(); //careful not to use this pointer after destruction of CDynamicOptions object
}

//////////////////////////////////////////////////////////
char* CDynamicOptions::strcasestr(const char* haystack, const char* needle)
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