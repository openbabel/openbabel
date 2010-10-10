/**********************************************************************
Copyright (C) 2006 Chris Morley

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
#include <fstream>
#include <string>

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBCOMMON
  #define OBCOMMON
#endif

using namespace std;
namespace OpenBabel
{

///Returns true if character is not one used in an InChI.
bool isnic(char ch)
{
  //This set of characters could be extended
  static std::string nic("\"\'\\@<>!$%&{}[]");
  return ch<0 || nic.find(ch)!=std::string::npos;
};

/// @brief Reads an InChI (possibly split) from an input stream and returns it as unsplit text.
/// The input stream is left after the end of the extracted InChI ready to look for the next one.
std::string GetInChI(std::istream& is);

/*!
This function recovers a normal InChI from an input stream which
contains other arbitrary text. The InChI string can have
extraneous characters inserted, for example because of word wrapping,
provided it follows certain rules.

When this file (getinchi.cpp) is read, 15 InChIs will be extracted, e.g.
 babel -iinchi getinchi.cpp -osmi

Inside an InChI string ignore anything between < and >
This means that an InChI string can be split up by inserting any number of <br /> elements:
InChI=1/C18H25NO6S/c1-14-9-11-15(12-10-14)26(22,23)19(17(21)25-18(2,3)4)13<br />-7-6-8-16(20)24-5/h6,8-12H,7,13H2,1-5H3/b8-6-

Any whitespace after the > is also ignored, so that newline characters can be added:
InChI=1/C29H33NO4Si/c1-5-32-28(31)26-25(34-27(30-26)22-15-9-6-10-16-22)<br />
21-33-35(29(2,3)4,23-17-11-7-12-18-23)24-19-13-8-14-20-24<br />
/h6-20,25-26H,5,21H2,1-4H3/t25-,26-/m0/s1

A second consecutive <...> element ends an unquoted InChI string:
<p>
<small>InChI=1/C47H58N2O10SSi/c1-10-56-43(51)47(36(32-41(50)55-9)30-31-49(44(52)59-45<br />
(3,4)5)60(53,54)37-28-26-34(2)27-29-37)40<br />
(58-42(48-47)35-20-14-11-15-21-35)33-57-61(46(6,7)8,38-22-16-12-17-23-38)39-24-<br />
18-13-19-25-39/h11-29,36,40H,10,30-33H2,1-9H3<br />
/t36-,40-,47-/m0/s1</small>
</p>

  Dmitrii Tchekhovskoi made a proposal for "InChI hyphenation" or "quoted InChI".
  http://sourceforge.net/mailarchive/forum.php?thread_id=10200459&forum_id=45166
  This proposal has not been followed up probably because InChKey was introduced.

  However this function GetInChI() parses quoted InChIs of this form.
  It also extends this proposal, allowing a wider range of corrupted InChIs to be accepted.

The original proposal was essentially:
- When an InChI string is enclosed by " quote characters,
  any whitespace characters it contains (including new lines) are
  ignored.
- Other extraneous strings can also be ignored, but this
  is system dependent.
- The "InChI=" cannot be split.

The extensions are:
- The character that encloses a quoted InChI does not have to be "
  and can be any character that is not used in InChI - a NIC
  [never miss the opportunity for a TLA!]. This means that
  conflicts in systems which have other uses for the quote character
  can be avoided.
  As a special case, '>' is not allowed as a quote character because InChI
  strings in HTML commonly start after <...> elements.
- As well as whitespace characters (which are ignored), a quoted
  InChI can contain an extraneous string which starts and ends with
  a NIC. This allows inserted strings like <br /> to be ignored.
  However, only one such extraneous string is allowed.
- There are no restrictions on splitting "InChI=" by whitespace
  characters, allowing a minimum column width of 1.
  If the splitting were by an extraneous string the minimum column
  width is 2.

The following are some examples of split InChIs.

First two unbroken examples, the first is unquoted
InChI=1/CH4/h1H4 methane
"InChI=1/C4H10O/c1-3-5-4-2/h3-4H2,1-2H3" diethyl ether

Multiple white space splitting
@InChI=1/C15H14O3/c1-11(15(16)17)18-
14-10-6-5-9-13(14)12-7-3-2-4-8-12/h2
-11H,1H3,(H,16,17)@

Split with extraneous text, which starts and ends with a non-InChI character
'InChI=1/C2H6O/c1-2-<br />3/h3H,2H2,1H3'

Table with wrapped InChI column. (View with fixed font.)

'InChI=1/CH4/h1H4'     !flammable!
'InChI=1/C2H2O4/c3-1   !toxic!
(4)2(5)6/h(H,3,4)(H,
5,6)'
'InChI=1/CH4O/c1-2/h   !flammable! !toxic!
2H,1H3'
'InChI=1/H2O/h1H2'
'InChI=1/C10H5ClN2/c   !no information!
11-10-4-2-1-3-9(10)5
-8(6-12)7-13/h1-5H'

Quoted text in emails (but InChI is preserved after one break only).
> "InChI=1/C4H7N3OS/c1-7(8)4-9-5-2-3-6-9/h
> 2-4,8H,1H3/p+1/fC4H8N3OS/h5H/q+1/t9?"
>> "InChI=1/C4H7N3OS/c1-7(8)4-9-5-2-3-6-9/
>> h2-4,8H,1H3/p+1/fC4H8N3OS/h5H/q+1/t9?"

Column width can be 1 if there is no extraneous text other than whitespace.
(When there is an extraneous string with NICs the minimum column width is 2).
'
I
n
C
h
I
=
1
/
C
l
H
/
h
1
H
'
*/

string GetInChI(istream& is)
{
  string prefix("InChI=");
  string result;
  enum statetype {before_inchi, match_inchi, unquoted, quoted};
  statetype state = before_inchi;
  char ch, lastch=0, qch=0;
  size_t split_pos = 0;
  bool inelement=false, afterelement=false;

  while((ch=is.get())!=EOF)
  {
    if(state==before_inchi)
    {
      if(ch>=0 && !isspace(ch))
      {
        if(ch==prefix[0])
        {
          result += ch;
          state = match_inchi;
          qch = lastch;
        }
      }
      lastch = ch;
    }

    else if(ch=='<')
    {
      // Ignore the content of any <...> elements
      // But a second consecutive  <...> element terminates an unquoted InChI
      if(afterelement && state==unquoted)
          return result;
      inelement=true;
    }
    else if(inelement)
    {
      if(afterelement)
      {
        //Now  reading after a <...> inserted in the InChI string
        //Neglect whitespace, but any other character reverts to normal InChI parsing
        if(ch<0 || !isspace(ch))
        {
          is.unget();
          afterelement=false;
          inelement=false;
        }
      }
      else
      {
        if(ch=='>')
          afterelement=true; //look for whitespace after end of element
      }
    }

    else if(ch>=0 && isspace(ch))
    {
      if(state==unquoted)
        return result;
    }

    else if(isnic(ch))
    {
      if(ch==qch && state!=match_inchi)
        return result;
      if(split_pos!=0)
        result.erase(split_pos);
      split_pos = result.size();
    }

    else
    {
      result += ch;
      if(state==match_inchi)
      {
        if(prefix.compare(0,result.size(),result)==0) //true if correct
        {
          if(result.size()==prefix.size())
            state = isnic(qch)&& qch!='>' ? quoted : unquoted;
        }
        else
        {
          is.unget(); //It may be the start of a real "InChI="
          result.erase();
          state = before_inchi;
        }
      }
    }
  }
  return result;
}

} //namespace
