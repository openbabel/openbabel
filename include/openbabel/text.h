/**********************************************************************
text.h - Declaration and Implementation of OBText

Copyright (C) 2008 by Chris Morley

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

#ifndef OB_TEXT_H
#define OB_TEXT_H

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>

namespace OpenBabel
{
  //! \class OBText text.h <openbabel/text.h>
  //! \brief An object containing just text
class OBText : public OBBase
{
private:
  std::string txt;
public:
  //Constructors
  OBText(){}
  OBText(const std::string& text) :txt(text) {}

  ///\return all the text
  std::string GetText()const { return txt; }

  /**\return text from position \param pos up to, but not including,
  the line containing the next occurrence of "OPENBABEL_INSERT".
  \param pos is updated to the start of the next line.
  If "OPENBABEL_INSERT" is not found, and \param ToInsertOnly is false
  the text up to the end of the file is returned and \param pos is set to 0.
  If "OPENBABEL_INSERT" is not found, and \param ToInsertOnly is true
  an empty string is return and \param pos is unchaged.

  Inserting OpenBabel output into boilerplate text.
  Suppose you wanted to insert XML output from OB into into a template XML document
  using the babel interface
    babel template.text inputfile.xxx outputfile.yyy
  The template file could contain a line:
    <!-- OPENBABEL_INSERT here -->
  and still be well-formed XML.
  The template file would be read by TextFormat and passed to the output as an OBText object.
  This could be processed in the output format's WriteChemObject() or WriteMolecule()
  in the following way (see cmlreactformat.cpp)
  <code>
      OBText* ptext = dynamic_cast<OBText*>(pOb);
      if(ptext) {
        string::size_type pos = 0;
        *pConv->GetOutStream() << ptext->GetText(pos); //Output text up to insertion point
        _text = ptext->GetText(pos); //Save text after insertion point to be output at the end
      }
  </code>
  **/
  std::string GetText(std::string::size_type& pos, bool ToInsertOnly=false) const
  {
    std::string::size_type oldpos = pos;
    std::string::size_type newpos = txt.find("OPENBABEL_INSERT", pos);
    if(newpos== std::string::npos)//not found: return rest of txt
    {
      if(ToInsertOnly)
        return("");
      pos = 0;
      return txt.substr(oldpos);
    }

    newpos = txt.rfind('\n', newpos); //to end of previous line
    pos = txt.find("\n", newpos+1)+1; //to past end of line, or 0
    return txt.substr(oldpos, newpos-oldpos);
  }

  void SetText(const std::string& text){ txt = text; }
};

}//namespace
#endif //OB_TEXT_H

//! \file text.h
//! \brief Handle generic non-chemical text processing
