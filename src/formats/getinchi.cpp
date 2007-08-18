/**********************************************************************
Copyright (C) 2006 Chris Morley

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

#include <openbabel/babelconfig.h>
#include <fstream>
#include <string>
using namespace std;
namespace OpenBabel
{

///Returns true if character is not one used in an InChI.
bool inline isnic(char ch)
{
  //This set of characters could be extended
  static string nic("\"\'\\@<>!$%&{}[]");
  return ch<0 || nic.find(ch)!=string::npos;
}

/// @brief Reads an InChI (possibly split) from an input stream and returns it as unsplit text.
/*!
  This function recovers a normal InChI from an input stream which
  contains other arbitary text. The InChI string can have
  extraneous characters inserted, for example because of word wrapping,
  provided it follows certain rules.

  Dmitrii Tchekhovskoi made a proposal for "InChI hyphenation".
  http://sourceforge.net/mailarchive/forum.php?thread_id=10200459&forum_id=45166
  The function here is consistent with this proposal but extends
  it, allowing a wider range of corrupted InChIs to be accepted.

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
- As well as whitespace characters (which are ignored), a quoted
  InChI can contain an extraneous string which starts and ends with
  a NIC. This allows inserted strings like <br /> to be ignored.
  However, only one such extraneous string is allowed.
- There are no restrictions on splitting "InChI=" by whitespace
  characters, allowing a minimum column width of 1.
  If the splitting were by an extraneous string the minimum column
  width is 2.

The following are some examples of split InChIs. 
OpenBabel will find and convert 12 InChIs
in this file, e.g. babel -iinchi getinchi.cpp -osmi

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
string OBAPI GetInChI(istream& is);

string GetInChI(istream& is)
{
  string prefix("InChI=");
  string result;
  enum statetype {before_inchi, match_inchi, unquoted, quoted};
  statetype state = before_inchi;
  char ch, lastch=0, qch=0;
  size_t split_pos = 0;

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
        lastch = ch;
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
            state = isnic(qch) ? quoted : unquoted;
        }
        else
        {
          result.erase();
          state = before_inchi;
        }
      }
    }
  }
  return result;
}

} //namespace
