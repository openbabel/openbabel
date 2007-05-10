/**********************************************************************
descriptor.h - Base class for molecular descriptors
 
Copyright (C) 2007 by Chris Morley
 
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

#ifndef OB_DESCRIPTOR_H
#define OB_DESCRIPTOR_H

#include <string>
#include <sstream>
#include <limits>

#include <openbabel/plugin.h>

namespace OpenBabel
{
class OBBase; //Forward declaration; used only as pointer.

class OBERROR OBDescriptor : public OBPlugin
{
  MAKE_PLUGIN(OBDescriptor)

  public:
    const char* TypeID(){return "descriptors";};

  ///returns the value of a numeric descriptor
  virtual double Predict(OBBase* pOb){return std::numeric_limits<double>::quiet_NaN();}

  ///returns the value of the descriptor and adds it to the object's OBPairData
  double PredictAndSave(OBBase* pOb);

  ///Provides a string value for non-numeric descriptors and returns NaN, or a string representation and returns a numeric value
  virtual double GetStringValue(OBBase* pOb, std::string& svalue);

  ///Parses the filter stream for a relational expression and returns its result when applied to the chemical object
  virtual bool Compare(OBBase* pOb, std::istream& ss, bool noEval);
  
  /// Interprets the --filter option string and returns the combined result of all the comparisons it contains  
  static bool FilterCompare(OBBase* pOb, std::istream& ss, bool noEval);
  
  ///Reads list of descriptor IDs and calls PredictAndSave() for each.
  static void AddProperties(OBBase* pOb, const std::string& DescrList);

  ///Deletes all the OBPairDatas whose attribute names are in the list (if they exist).
  static void DeleteProperties(OBBase* pOb, const std::string& DescrList);

protected:
  ///Read an identifier from the filter string
  static std::string GetIdentifier(std::istream& optionText);

  static double ParsePredicate(std::istream& optionText, char& ch1, char& ch2, std::string& svalue);

  static bool GenericDataCompare(std::string& ID, OBBase* pOb, std::istream& optionText, bool noEval);

    ///Reads a string from the filter stream, optionally preceded by = or !=
  ///Returns false if != operator found, and true otherwise.
  static bool ReadStringFromFilter(std::istream& ss, std::string& result);

  ///Makes a comparison using the operator and a string read from the filter stream with a provided string.
  ///Returns the result of the comparison.
  static bool CompareStringWithFilter(std::istream& optionText, std::string& s, bool noEval);

  //Treats _ as not a punctuation chaarcter
  static bool ispunctU(char ch)
  {
    return ispunct(ch) && ch!='_';
  }
};

template <class T>
static bool DoComparison(char ch1, char ch2, T& val, T& filterval)
{
  switch(ch1)
  {
  case (0):  //no comparison operator is same as =
  case('='):
  return val==filterval; //**needs a better floating point comparison**	
  case('!'):
  return val!=filterval; //**needs a better floating point comparison**	
  case('>'):              
  if(ch2=='=')
    return val>=filterval;
  else
    return val>filterval;
  case('<'):              
  if(ch2=='=')
    return val<=filterval;
  else
    return val<filterval;
  }
  return false;
}

}//namespace
#endif