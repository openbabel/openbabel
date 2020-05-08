/**********************************************************************
optswx.cpp -  Constructs wxWidgets Option controls from description text

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
#include <openbabel/babelconfig.h>
#include <stdwx.h>
#include <sstream>
#include <utility>
#include <optswx.h>

/*
The string returned by OBFormat::Description() gives information about the
format and the commandline options that can be used with it. More general
options have similar descriptions in for instance in
OBConversion::TargetClassDescription(). These strings are parsed in
OBGUIFrame::OnChangeFormat and the following routine so that the GUI can
present these options as checkboxes, editboxes, etc. The method was originally
designed to use the existing descriptions as far as possible, so that the
information is extracted without having to add special codes, etc. It is
fairly unobtrusive, but liable to breaking because of alterations which
had more significance than the person making then was expecting.

For format descriptions options for Read and Write are extracted with
separate calls to DynOptionswx::Construct(). The list of options starts with
a line containing "Options". Those for reading have "Read" or "Input"
preceding this on the same line and those for writing have "Write", "Output"
or nothing preceding it on the same line. These are all case insensitive.
The options themselves are one to a line on subsequent lines, ending with a
blank line or the end of the string.

Each option has a name, which can be a single letter. Any punctuation
preceding the name is ignored.

If the first character after the name, other than space or tab, is a
punctuation character, then the option will be displayed as an edit box,
to be used to enter parameters.

--property <attrib> <value> add or replace a property\n

The name of this parameter, e.g. "attrib" and "value" in the example,
is not displyed in the GUI. Additional editboxes are displayed if the first
character after the next whitespace is punctuation. If it is a letter, it is
the start of the caption text which continues to the end of the line (but see
below). The caption is displyed after the edit box(es).Edit boxes for options
with multicharacter names are larger than for those with single character names.

If the last character of the line is ':' then a full-line edit box is
displayed and the caption is displayed before it.

If the text contains the word "default", the subsequent word is written into the
editbox. It will have some punctuation characters stripped from its start and end,
so that
q <day> Starting day, default <Sunday>
 will display as an edit box containing "Sunday" and caption "Starting day".

If the first character after the option name, other than space or tab, is a
letter, then a checkbox is usually displayed.
  c  continuous output: no formatting\n

If the text contains "default", the checkbox is displayed checked and the
caption contains only the text up to this point. If the follow word is
"no", "not or "none" the checkbox is unchecked.

A set of radio buttons are displayed, rather than a checkbox if "or" is the
last word of the caption. Including "default" in the caption checks that button.
For example
1  output CML V1.0  or \n \
2  output CML V2.0 (default)\n \

*/
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
  const int MAXLEADINGSPACES = 3; //not an option line if more
  const int ONE=1;
  const int FOUR=4;
  bool NextIsRadio=false;
  char* pNewStr = new char[strlen(OptionsText)+1];
  strcpy(pNewStr,OptionsText); //Make a working copy
  char* p = pNewStr;

  bool OptionsFound=false;

  std::string qualifiedOptions;
  if(StartText)
    qualifiedOptions = StartText;
  qualifiedOptions += " options";
  p = strcasestr(p, qualifiedOptions.c_str());
  if(p)
  {
    OptionsFound=true;
    p = strchr(p,'\n')+1; //options start on next line
    while(p && p-1 && *p && *p!='\n') //loop for all options until blank line
    {
      int ProvideEditCtl=0;
      bool ProvideExtraCheckbox=false;

      char* plinestart = p;
      while(*p && *p!='\n' && !isalnum(*(p++))); //skip space and punctuation except \n
      if(*p=='\n')
      {
        plinestart = ++p; //reset start of line
        continue;
      }
      //Ignore lines which start with more than MAXLEADINGSPACES whitespace chars(tab is 1 char!)
      if(p - plinestart > MAXLEADINGSPACES)
      {
        p = strchr(p,'\n')+1;
        continue;
      }
      //Ignore lines containing "not displayed in GUI"
      char* lineend = strchr(p, '\n');
      if(p)
      {
        *lineend = '\0';
        char* pw = strstr(p, "not displayed in GUI");
        *lineend = '\n';
        if(pw)
        {
          p = lineend+1;
          continue;
        }
      }

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

      while(*p>=0 && isspace(*p++));//skip white space
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
          endch= ']';
          ProvideExtraCheckbox=true;
          break;
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
        char* pdefWord=nullptr;
        if(pdef=strcasestr(pCaption,"default"))
        {
          //Put the next word in the editbox
          char* tok=strtok(pdef," :-\t</\"\'");
          pdefWord = strtok(nullptr," :-</\t;,\"\'>");

          //delete caption after default or after <default etc
          *pdef='\0';
          while(isspace(*(--pdef)));
          if(!strpbrk(pdef,"([<{-;")) pdef++ ;
          *pdef='\0';
        }
        wxStaticText* pEdCaption = new wxStaticText(parent,wxID_STATIC,wxString(pCaption, wxConvUTF8));
        OptionMap.push_back(std::make_pair(wxString(),pEdCaption));//string is empty for a caption: not an option

        //Edit boxes for multicharacter options are larger
        const int EDWIDTH = 60; //now all are the same size. oname.size()>1? 60 : 40;
        wxTextCtrl* pEd;

        //Make a large edit box on the next line if last char is of the caption is :
        bool BigEdit =(pCaption[strlen(pCaption)-1] == ':');
        if(BigEdit)
        {
          sizer->Add(pEdCaption,0,wxEXPAND|wxTOP,FOUR);
          while(ProvideEditCtl--)
          {
            pEd = new wxTextCtrl(parent,wxID_ANY,wxEmptyString,
                wxDefaultPosition,wxSize(EDWIDTH,18));
            OptionMap.push_back(std::make_pair(oname,pEd));
            if(ProvideEditCtl)
              oname = _T(' ') + oname;//editboxes except the first have name preceded by one or more spaces
          }
          sizer->Add(pEd,0,wxEXPAND|wxTOP,ONE);
        }
        else
        {
          wxBoxSizer* pEdSizer = new wxBoxSizer(wxHORIZONTAL);
          if(ProvideExtraCheckbox)
          {
            wxControl* pChk = new wxCheckBox(parent,wxID_ANY,_T(" "));
            OptionMap.push_back(std::make_pair(oname,pChk));
            pEdSizer->Add(pChk,0,wxALIGN_CENTER_VERTICAL,FOUR);
            oname = _T(' ') + oname;
          }
          while(ProvideEditCtl--)
          {
            pEd = new wxTextCtrl(parent,wxID_ANY,wxEmptyString,
                wxDefaultPosition,wxSize(EDWIDTH,18));
            OptionMap.push_back(std::make_pair(oname,pEd));
            if(ProvideEditCtl)
              oname = _T(' ') + oname;//editboxes except the first have name preceded by one or more spaces
            pEdSizer->Add(pEd,0);
          }
          pEdSizer->Add(pEdCaption,1,wxLEFT|wxALIGN_CENTER_VERTICAL,FOUR);
          sizer->Add(pEdSizer,0,wxEXPAND|wxTOP,FOUR);
          Sizers.push_back(pEdSizer);
        }
        pEd->AppendText(wxString(pdefWord, wxConvUTF8));
      }
      else
      {
        // Checkbox, or radio button ("or" is the last word in the caption && no letters after " or")
        wxControl* pChk;

        //First see if should be checked
        //If 'default appears in caption set the checkbox unless it is followed
        //by one of 'no' 'not' or 'none'
        char* pdef=strcasestr(pCaption,"default");
        bool SetChk=(pdef!=nullptr);
        if(SetChk)
          strcasecmp(pdef,"no") && !strcasecmp(pdef,"not") && strcasecmp(pdef,"none");

        char* por=strcasestr(pCaption," or");
        bool HasOr = (por && !strpbrk(por+3,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"));

        if(NextIsRadio || HasOr)
        {
          unsigned int style = NextIsRadio ? 0 : wxRB_GROUP; //Style only on first radiobutton
          pChk = new wxRadioButton(parent,wxID_ANY,wxString(pCaption, wxConvUTF8),
                wxDefaultPosition, wxDefaultSize, style);
          NextIsRadio = HasOr;
          ((wxRadioButton*)pChk)->SetValue(SetChk);
        }
        else
        {
          pChk = new wxCheckBox(parent,wxID_ANY,wxString(pCaption, wxConvUTF8));
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

    wxString oname = itr->first;
    wxString txt;

    wxCheckBox* pChk = dynamic_cast<wxCheckBox*> (itr->second);
    if(pChk)
    {
      if(!pChk->IsChecked())
      {
        // if a checkbox is not checked, ignore the subsidiary editboxes also
        while(!(++itr)->first.empty()&& itr->first[0]==_T(' '));
        --itr;
        continue;
      }
    }
    else
    {
      wxRadioButton* pRadio = dynamic_cast<wxRadioButton*> (itr->second);
      if(pRadio)
      {
        if(pRadio->GetValue())
          continue;
      }
      else
      {
        wxTextCtrl* pText = dynamic_cast<wxTextCtrl*> (itr->second);
        if(pText)
        {
          txt = pText->GetValue();
          if(txt.IsEmpty()) continue;
          oname = itr->first;
        }
      }
    }

    //Get the contents of subsequent editboxes
    OMapType::iterator itr2 = itr;
    while(++itr2!= OptionMap.end())
    {
      if((itr2->first).empty() || itr2->first[0]!=_T(' ')) //subsequent editboxes have the name preceded by a space
        break;
      txt = txt + _T(' ') + static_cast<wxTextCtrl*>(itr2->second)->GetValue();
      ++itr;
    }
    txt.Trim(true); txt.Trim(false);
    Conv.AddOption(oname.mb_str(), opttyp, txt.mb_str());
    ++count;

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
  return nullptr;
}
