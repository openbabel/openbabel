/**********************************************************************
descriptor.h - Base class for molecular descriptors

Copyright (C) 2007 by Chris Morley

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

#ifndef OB_DESCRIPTOR_H
#define OB_DESCRIPTOR_H

#include <string>
#include <sstream>
#include <limits>

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>

namespace OpenBabel
{
class OBBase; //Forward declaration; used only as pointer.

// Class introduction in descriptor.cpp
class OBAPI OBDescriptor : public OBPlugin
{
  MAKE_PLUGIN(OBDescriptor)

public:
  const char* TypeID(){return "descriptors";};

  /// \return the value of a numeric descriptor
  virtual double Predict(OBBase* /* pOb */, std::string* /* param */ =nullptr)
  {return std::numeric_limits<double>::quiet_NaN();}

  /// \return the value of the descriptor and adds it to the object's OBPairData
  double PredictAndSave(OBBase* pOb, std::string* param=nullptr);

  ///Provides a string value for non-numeric descriptors and returns NaN, or a string representation and returns a numeric value
  virtual double GetStringValue(OBBase* pOb, std::string& svalue, std::string* param=nullptr);

  ///Parses the filter stream for a relational expression and returns its result when applied to the chemical object
  virtual bool Compare(OBBase* pOb, std::istream& ss, bool noEval, std::string* param=nullptr);

  ///Write information on a plugin class to the string txt.
  ///If the parameter is a descriptor ID, displays the verbose description for that descriptor only
  /// e.g. babel -L descriptors HBA1
  virtual bool Display(std::string& txt, const char* param, const char* ID=nullptr);

  /// Comparison of the values of the descriptor. Used in sorting.
  /// Descriptors may use more complicated ordering than this default (e.g.InChIFilter)
  virtual bool Order(double p1, double p2){ return p1<p2; }
  virtual bool Order(std::string s1, std::string s2){ return s1<s2; }

  /// Interprets the --filter option string and returns the combined result of all the comparisons it contains
  static bool FilterCompare(OBBase* pOb, std::istream& ss, bool noEval);

  ///Reads list of descriptor IDs and calls PredictAndSave() for each.
  static void AddProperties(OBBase* pOb, const std::string& DescrList);

  ///Deletes all the OBPairDatas whose attribute names are in the list (if they exist).
  static void DeleteProperties(OBBase* pOb, const std::string& DescrList);

  ///Reads list of descriptor IDs and OBPairData names and returns a list of values,
  ///each precede by a space or the first character in the list if it is whitespace or punctuation.
  static std::string GetValues(OBBase* pOb, const std::string& DescrList);

  ///Read an identifier and its parameter from the filter string.
  static std::pair<std::string, std::string> GetIdentifier(std::istream& optionText);

protected:

  static double ParsePredicate(std::istream& optionText, char& ch1, char& ch2, std::string& svalue);

  ///Reads a string from the filter stream, optionally preceded by = or !=
  /// \return false if != operator found, and true otherwise.
  static bool ReadStringFromFilter(std::istream& ss, std::string& result);

  ///Makes a comparison using the operator and a string read from the filter stream with a provided string.
  /// \return the result of the comparison and true if NoCompOK==true and there is no comparison operator.
  static bool CompareStringWithFilter(std::istream& optionText, std::string& s, bool noEval, bool NoCompOK=false);

  // Treats _ as not a punctuation character and since 2.3.2 also $ # and %
  static bool ispunctU(char ch)
  {
    return ispunct(ch) && ch!='_' && ch!='$' && ch!='#' && ch!='%';
  }

  /// \return true if s (with or without _ replaced by spaces) is a PairData attribute. On return s is the form which matches.
  static bool MatchPairData(OBBase* pOb, std::string& s);
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

//! \file descriptor.h
//! \brief Base class for molecular descriptors
